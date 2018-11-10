package EhdFlow;

	import java.io.BufferedOutputStream;
	import java.io.FileNotFoundException;
	import java.io.FileOutputStream;
	import java.io.PrintStream;
	import java.util.Scanner;

	public class EhdFlow {

		public static void main(String[] args) {

					Scanner scan = new Scanner(System.in);

					System.out.print("領域の設定");
					System.out.println("X方向の領域XL = 0.006[m]");
					double XL = 0.006;

					System.out.println("Y方向の長さYL = 0.002[m]");
					double YL = 0.002;

					System.out.println("移流時間t = 0.5[s]");
					double t = 0.5;

					System.out.println("タイムステップdt = 0.000001[s]");
					double dt = 0.000001;

					System.out.println("緩和係数a = 0.5");
					double a = 0.5;
	//電位の収束条件
					System.out.println("収束条件J1 = 0.00000001");
					double J1 = 0.00000001;
	//流速の収束条件
					System.out.println("収束条件J2 =0.000001");
					double J2 =0.000001;

					//要素の設定
					System.out.println("要素の大きさdx=dy=h= 0.00005");
					double h = 0.00005;

					//動粘度(液晶)[m2/s]
					double ν = 0.000096;
					//密度(液晶)[kg/m3]
					double ρ = 950;
					//誘電率[F/m]
					double ε=1.00058*9.9*0.00000000001;
					//イオン移動度[m2/Vs]
					double μi = 0.000000001;
					//導電率[A/Vm]
					double σ =  0.000000001;
					//イオン拡散係数
					double Di = 0;
					//要素の設定

					//電極枚数
					int n = 2;
					//電圧[V]
					double V[] = new double[n];
					//電極長さ[m]
					double EL[] = new double[n];
					//ギャップ[m]
					double GL[] = new double[n];

					for(int i=0;i<n;i++){
						System.out.printf("V[%d][V]:",i);
						V[i]=scan.nextDouble();
					}
					for(int i=0;i<n;i++){
						System.out.printf("EL[%d][m]:",i);
						EL[i]=scan.nextDouble();
					}
					for(int i=0;i<n-1;i++){
						System.out.printf("GL[%d][m]:",i);
						GL[i]=scan.nextDouble();
					}

					//要素の大きさ
					System.out.println("要素数");
					int Xn = (int)(XL/h)+2;
					int Yn = (int)(YL/h)+2;

					//0を含めてXn+1,Yn+1
					System.out.printf("Xn=%d\n",Xn+2);
					System.out.printf("Yn=%d\n",Yn+2);

					//装置の大きさ
					int Xnn=0;
					int ds=0;
					int ELn[] = new int[n];
					int GLn[] = new int[n];

					for(int i=0;i<n;i++){
					ELn[i] = (int)(EL[i]/h);
					if(i==n-1)break;
					GLn[i] = (int)(GL[i]/h);
					}

					for(int i=0;i<n;i++){
						ds=ds+ELn[i];
						if(i==n-1) break;
						ds=ds+GLn[i];
					}

					//繰り返し計算(静電ポテンシャル)
					double Φ[][] = new double[Xn][Yn];
					double preΦ[][] = new double[Xn][Yn];
					double q[][] = new double[Xn][Yn];
					double preq[][] = new double[Xn][Yn];
					double preEx[][] = new double[Xn][Yn];
					double preEy[][] = new double[Xn][Yn];
					double term1x[][] = new double[Xn][Yn];
					double term2x[][] = new double[Xn][Yn];
					double term3x[][] = new double[Xn][Yn];
					double term4x[][] = new double[Xn][Yn];
					double term1y[][] = new double[Xn][Yn];
					double term2y[][] = new double[Xn][Yn];
					double term3y[][] = new double[Xn][Yn];
					double term4y[][] = new double[Xn][Yn];
					double jx[][] = new double[Xn][Yn];
					double jy[][] = new double[Xn][Yn];
					double dis1;
					double dis2;

					//EHD
					double Ex[][] = new double[Xn][Yn];
					double Ey[][] = new double[Xn][Yn];
					double absE[][] = new double[Xn][Yn];
					double fx[][] = new double[Xn][Yn];
					double fy[][] = new double[Xn][Yn];

					//静電場の計算
					//電位
					do{
						dis1=0;
						for(int i=0;i<Xn;i++){
							for(int j=0;j<Yn;j++){
						preΦ[i][j]=Φ[i][j];
							}
						}
						//入り口
						for(int j=1;j<Yn-1;j++){
							Φ[0][j]=preΦ[1][j];
						}
						//出口
						for(int j=1;j<Yn-1;j++){
							Φ[Xn-1][j]=preΦ[Xn-2][j];
						}
						//境界条件（上壁面）
						for(int i=0;i<Xn;i++){
							Φ[i][Yn-1]=preΦ[i][Yn-2];
							}
						//境界条件（下壁面）
						for(int i=0;i<(Xn-ds)/2+1;i++){
							Φ[i][0]=preΦ[i][1];
						}
						Xnn=(Xn-ds)/2+1;
						for(int k=0;k<n;k++){
						for(int i=Xnn;i<=Xnn+ELn[k];i++){
							Φ[i][0]=V[k];
							}

						Xnn=Xnn+ELn[k];

						if(k==n-1)break;

						for(int i=Xnn+1;i<Xnn+GLn[k];i++){
							Φ[i][0]=preΦ[i][1];
							}

						Xnn=Xnn+GLn[k];
						}
						for(int i=Xnn+1;i<Xn;i++){
							Φ[i][0]=preΦ[i][1];
							}

						for(int i=1;i<Xn-1;i++){
							for(int j=1;j<Yn-1;j++){
								Φ[i][j] = (Φ[i+1][j]+Φ[i-1][j]+Φ[i][j+1]+Φ[i][j-1])/4;
							if(dis1<Math.abs(Φ[i][j]-preΦ[i][j])/Φ[i][j])dis1=Math.abs(Φ[i][j]-preΦ[i][j])/Φ[i][j];
							}
						}

						System.out.println(dis1);
					}while(J1<dis1);

					//電場のベクトルを求める[V/m]
					for(int i=0;i<Xn;i++){
						for(int j=0;j<Yn;j++){
							if(i==Xn-1) Ex[i][j]=0;
							else Ex[i][j]=-(Φ[i+1][j]-Φ[i][j])/h;
						}
					}
					for(int i=0;i<Xn;i++){
						for(int j=0;j<Yn;j++){
							if(j==Yn-1)Ey[i][j]=0;
							else Ey[i][j]=-(Φ[i][j+1]-Φ[i][j])/h;
						}
					}

					//繰り返し計算
					double u[][] = new double[Xn][Yn];
					double v[][] = new double[Xn][Yn];
					double p[][] = new double[Xn][Yn];
					double preu[][] = new double[Xn][Yn];
					double prev[][] = new double[Xn][Yn];
					double prep[][] = new double[Xn][Yn];
					double sa;
					int tcount=0;
					double dp;
					double D;
					double fux;
					double fvx;
					double visx;
					double visy;
					double fuy;
					double fvy;
					double advu;
					double advv;

					//保存用
					int memot[]=new int[5];
					double memou[][][] = new double[Xn][Yn][5];
					double memov[][][] = new double[Xn][Yn][5];
					double memop[][][] = new double[Xn][Yn][5];
					double memoΦ[][][] = new double[Xn][Yn][5];
					double memoEx[][][] = new double[Xn][Yn][5];
					double memoEy[][][] = new double[Xn][Yn][5];
					double memoabsE[][][] = new double[Xn][Yn][5];
					double memofx[][][] = new double[Xn][Yn][5];
					double memofy[][][] = new double[Xn][Yn][5];
					double memoq[][][] = new double[Xn][Yn][5];
					double memojx[][][] = new double[Xn][Yn][5];
					double memojy[][][] = new double[Xn][Yn][5];
					double memoAx[][][] = new double[Xn][Yn][5];
					double memoAy[][][] = new double[Xn][Yn][5];


					//記録時間
					memot[0]=(int) (1);
					memot[1]=(int) (2);
					memot[2]=(int) ((t-3*dt)/dt);
					memot[3]=(int) ((t-2*dt)/dt);
					memot[4]=(int) ((t-dt)/dt);

					//初期化
					for(int i=0;i<Xn;i++){
						for(int j=0;j<Yn;j++){
							u[i][j]=0;
							v[i][j]=0;
							p[i][j]=0;
						}
					}

					//座標
					double x[]=new double[Xn];
					double y[]=new double[Yn];
					for(int i=0;i<Xn;i++)x[i]=(i-1)*h+h/2;
					for(int j=0;j<Yn;j++)y[j]=(j-1)*h+h/2;

					//dt秒後の繰り返し計算
					do{
						tcount+=1;
						//各値の更新
						for(int i=0;i<Xn;i++){
							for(int j=0;j<Yn;j++){
								preu[i][j] = u[i][j];
								prev[i][j] = v[i][j];
								prep[i][j] = p[i][j];
								preΦ[i][j] = Φ[i][j];
								preq[i][j] = q[i][j];
								preEx[i][j] = Ex[i][j];
								preEy[i][j] = Ey[i][j];
							}
						}
						//電荷密度
						//x方向の電流
						for(int i=0;i<Xn-1;i++){
							for(int j=0;j<Yn-1;j++){
								if(i==Xn-1){
									jx[i][j]=0;
								}
								else {
									term1x[i][j]=(preq[i][j]+preq[i+1][j])/2*μi*preEx[i][j];
									term2x[i][j]=(preq[i][j]+preq[i+1][j])/2*preu[i][j];
									//term3x[i][j]=Di*(preq[i+1][j]-preq[i][j])/h;
									term4x[i][j]=σ*preEx[i][j];
									jx[i][j]=term1x[i][j]+term2x[i][j]/**-term3x[i][j]**/+term4x[i][j];
								}
							}
						}
						//y方向の電流
						for(int i=0;i<Xn;i++){
							for(int j=0;j<Yn;j++){
								if(j==Yn-1){
									jy[i][j]=0;
								}
								else{
									term1y[i][j]=(preq[i][j]+preq[i][j+1])/2*μi*preEy[i][j];
									term2y[i][j]=(preq[i][j]+preq[i][j+1])/2*prev[i][j];
									//term3y[i][j]=Di*(preq[i][j+1]-preq[i][j])/h;
									term4y[i][j]=σ*preEy[i][j];
									jy[i][j]=term1y[i][j]+term2y[i][j]/**-term3y[i][j]**/+term4y[i][j];
								}
							}
						}

						//電荷密度の更新
						for(int i=1;i<Xn-1;i++){
							for(int j=1;j<Yn-1;j++){
								q[i][j]=preq[i][j]-dt/h*(jx[i][j]-jx[i-1][j]+jy[i][j]-jy[i][j-1]);
							}
						}
						//境界条件
						//入り口
						for(int j=1;j<Yn-1;j++){
							q[0][j]=preq[1][j];
						}
						//出口
						for(int j=1;j<Yn-1;j++){
							q[Xn-1][j]=preq[Xn-2][j];
						}
						//境界条件（上壁面）
						for(int i=0;i<Xn;i++){
							q[i][Yn-1]=preq[i][Yn-2];
						}
						//境界条件（下壁面）
						for(int i=0;i<(Xn-ds)/2;i++){
							q[i][0]=preq[i][1];
							}
						Xnn=(Xn-ds)/2;
						for(int k=0;k<n;k++){
							for(int i=Xnn;i<=Xnn+ELn[k];i++){
								q[i][0]=ε*preEy[i][0];
							}
							Xnn=Xnn+ELn[k];
							if(k==n-1)break;
							for(int i=Xnn+1;i<Xnn+GLn[k];i++){
								q[i][0]=preq[i][1];
							}
							Xnn=Xnn+GLn[k];
						}
						for(int i=Xnn+1;i<Xn;i++){
							q[i][0]=preq[i][1];
						}
						do{
							dis2=0;
							for(int i=0;i<Xn;i++){
								for(int j=0;j<Yn;j++){
									preΦ[i][j]=Φ[i][j];
								}
							}
							//入り口
							for(int j=1;j<Yn-1;j++){
								Φ[0][j]=preΦ[1][j];
							}
							//出口
							for(int j=1;j<Yn-1;j++){
								Φ[Xn-1][j]=preΦ[Xn-2][j];
							}
							//境界条件（上壁面）
							for(int i=0;i<Xn;i++){
								Φ[i][Yn-1]=preΦ[i][Yn-2];
								}
							//境界条件（下壁面）
							for(int i=0;i<(Xn-ds)/2;i++){
								Φ[i][0]=preΦ[i][1];
								}
							Xnn=(Xn-ds)/2;
							for(int k=0;k<n;k++){
							for(int i=Xnn;i<=Xnn+ELn[k];i++){
								Φ[i][0]=V[k];
								}
							Xnn=Xnn+ELn[k];
							if(k==n-1)break;
							for(int i=Xnn+1;i<Xnn+GLn[k];i++){
								Φ[i][0]=preΦ[i][1];
								}
							Xnn=Xnn+GLn[k];
							}
							for(int i=Xnn+1;i<Xn;i++){
								Φ[i][0]=preΦ[i][1];
								}

							for(int i=0;i<Xn;i++){
								for(int j=0;j<Yn;j++){
							preΦ[i][j]=Φ[i][j];
								}
							}

							for(int i=1;i<Xn-1;i++){
								for(int j=1;j<Yn-1;j++){
									Φ[i][j]=(Φ[i+1][j]+Φ[i-1][j]+Φ[i][j+1]+Φ[i][j-1]+h*h*q[i][j]/ε)/4;
								if(dis2<Math.abs(Φ[i][j]-preΦ[i][j])) dis2=Math.abs(Φ[i][j]-preΦ[i][j]);
								}
							}
							System.out.println(dis2);
						}while(J1<dis2);

						//電場のベクトルを求める[V/m]
						for(int i=0;i<Xn;i++){
							for(int j=0;j<Yn;j++){
								if(i==Xn-1) Ex[i][j]=-(Φ[i][j]-Φ[i-1][j])/h;
								else Ex[i][j]=-(Φ[i+1][j]-Φ[i][j])/h;
							}
						}
						for(int i=0;i<Xn;i++){
							for(int j=0;j<Yn;j++){
								if(j==Yn-1)Ey[i][j]=-(Φ[i][j]-Φ[i][j-1])/h;
								else Ey[i][j]=-(Φ[i][j+1]-Φ[i][j])/h;
							}
						}
						//電場の大きさ
						for(int i=0;i<Xn;i++){
							for(int j=0;j<Yn;j++){
								if(i==0&&j==0)absE[i][j]=Math.sqrt((Ex[i][j]*Ex[i][j]+Ey[i][j]*Ey[i][j]));
								else if(i==Xn-1&&j==Yn-1)absE[i][j]=Math.sqrt((Ex[i-1][j]*Ex[i-1][j]+Ey[i][j-1]*Ey[i][j-1]));
								else if(j==0)absE[i][j]=Math.sqrt((Ex[i][j]+Ex[i-1][j])*(Ex[i][j]+Ex[i-1][j])+Ey[i][j]*Ey[i][j]);
								else if(i==0)absE[i][j]=Math.sqrt(Ex[i][j]*Ex[i][j]+(Ey[i][j]+Ey[i][j-1])*(Ey[i][j]+Ey[i][j-1]));
								else if(j==Yn-1)absE[i][j]=Math.sqrt((Ex[i][j]+Ex[i-1][j])*(Ex[i][j]+Ex[i-1][j])+Ey[i][j-1]*Ey[i][j-1]);
								else if(i==Xn-1)absE[i][j]=Math.sqrt(Ex[i-1][j]*Ex[i-1][j]+(Ey[i][j]+Ey[i][j-1])*(Ey[i][j]+Ey[i][j-1]));
								else absE[i][j]=Math.sqrt((Ex[i][j]+Ex[i-1][j])*(Ex[i][j]+Ex[i-1][j])+(Ey[i][j]+Ey[i][j-1])*(Ey[i][j]+Ey[i][j-1]));
							}
						}
						//力[N]
						for(int i=1;i<Xn-1;i++){
							for(int j=1;j<Yn-1;j++){
								//fx[i][j]=(ε*(absE[i][j]+absE[i][j])/2*(absE[i+1][j]-absE[i][j])/h)/ρ+(preq[i][j]+preq[i][j+1])*Ex[i][j]/(2*ρ);
								fx[i][j]=(ε*(absE[i][j]+absE[i+1][j])/2*(absE[i+1][j]-absE[i][j])/h)/ρ;
								//fx[i][j]=(preq[i][j]+preq[i][j+1])*Ex[i][j]/(2*ρ);
								if(j==Yn-2) fy[i][j]=0;
								//else fy[i][j]=(ε*(absE[i][j]+absE[i][j+1])/2*(absE[i][j+1]-absE[i][j])/h)/ρ+q[i][j]*(Ey[i][j]+Ey[i][j-1])/(2*ρ);
								else fy[i][j]=(ε*(absE[i][j]+absE[i][j+1])/2*(absE[i][j+1]-absE[i][j])/h)/ρ;
								//else fy[i][j]=(preq[i][j]+preq[i][j+1])*Ey[i][j]/(2*ρ);
							}
						}

						//流れの計算
						//入り口
						for(int j=1;j<Yn-1;j++){
							u[0][j]=preu[1][j];
							v[0][j]=0;
						}

						//出口
						for(int j=1;j<Yn-1;j++){
							u[Xn-1][j]=preu[Xn-2][j];
							if(j==Yn-2) v[Xn-1][j]=0;
							else v[Xn-1][j]=prev[Xn-2][j];
						}

						//境界条件（壁面）
						for(int i=0;i<Xn;i++){
							//底面
							u[i][0]=-preu[i][1];
							v[i][0]=0;
							//上面
							u[i][Yn-1]=-preu[i][Yn-2];
							v[i][Yn-2]=0;
							}

						for(int i=0;i<Xn;i++){
							for(int j=0;j<Yn;j++){
								preu[i][j] = u[i][j];
								prev[i][j] = v[i][j];
								prep[i][j] = p[i][j];
							}
						}

						//領域内
						for(int i=1;i<Xn-1;i++){
							for(int j=1;j<Yn-1;j++){
								advv=(prev[i][j]+prev[i][j-1]+prev[i+1][j]+prev[i+1][j-1])/4;
								//各項の計算
								fux=(preu[i][j]*(preu[i+1][j]-preu[i-1][j])-Math.abs(preu[i][j])*(preu[i+1][j]+preu[i-1][j]-2*preu[i][j]))/(2*h);
								fvx=(advv*(preu[i][j+1]+preu[i][j-1])-Math.abs(advv)*(preu[i][j+1]+preu[i][j-1]-2*preu[i][j]))/(2*h);
								visx=ν*(preu[i+1][j]+preu[i-1][j]+preu[i][j+1]+preu[i][j-1]-4*preu[i][j])/h/h;
								//uの仮値の決定
								u[i][j]=preu[i][j]+dt*(visx-fux-fvx-(prep[i+1][j]-prep[i][j])/(h*ρ)+fx[i][j]);
							}
						}

						for(int i=1;i<Xn-1;i++){
							for(int j=1;j<Yn-1;j++){
								if(j==Yn-2)v[i][j]=0;
								else{
								advu=(preu[i][j]+preu[i-1][j]+preu[i][j+1]+preu[i-1][j+1])/4;
								//各項の計算
								fuy=(advu*(prev[i+1][j]+prev[i-1][j])-Math.abs(advu)*(prev[i+1][j]+prev[i-1][j]-2*prev[i][j]))/(2*h);
								fvy=(prev[i][j]*(prev[i][j+1]-prev[i][j-1])-Math.abs(prev[i][j])*(prev[i][j+1]+prev[i][j-1]-2*prev[i][j]))/(2*h);
								visy=ν*(prev[i+1][j]+prev[i-1][j]+prev[i][j+1]+prev[i][j-1]-4*prev[i][j])/h/h;
								//vの仮値の計算
								v[i][j]=prev[i][j]+dt*(visy-fuy-fvy-(prep[i][j+1]-prep[i][j])/(h*ρ)+fy[i][j]);
								}
							}
						}

						do{
							sa=0;
								for(int i=1;i<Xn-1;i++){
									for(int j=1;j<Yn-1;j++){
										D=(u[i][j]-u[i-1][j]+v[i][j]-v[i][j-1])/h;
										if(Math.abs(D)>sa) sa=Math.abs(D);
										dp=-a*D/(2*dt/(h*h));
										p[i][j]=p[i][j]+dp;
										if(j==Yn-2){
											if(i==1){
												u[i][j]=u[i][j]+dt*dp/h;
												u[i-1][j]=u[1][j];
												v[i][j]=0;
												v[i][j-1]=v[i][j-1]-dt*dp/h;
											}
											if(i==Xn-1){
												u[i][j]=u[i-1][j];
												u[i-1][j]=u[i-1][j]-dt*dp/h;
												v[i][j]=0;
												v[i][j-1]=v[i-1][j-1];
											}
											else {
												u[i][j]=u[i][j]+dt*dp/h;
												u[i-1][j]=u[i-1][j]-dt*dp/h;
												v[i][j]=0;
												v[i][j-1]=v[i][j-1]-dt*dp/h;
												}
											}
										if(j==1){
											if(i==1){
												u[i][j]=u[i][j]+dt*dp/h;
												u[i-1][j]=u[1][j];
												v[i][j]=v[i][j]+dt*dp/h;
												v[i][j-1]=0;
												}
											if(i==Xn-1){
												u[i][j]=u[i-1][j];
												u[i-1][j]=u[i-1][j]-dt*dp/h;
												v[i][j]=v[i-1][j];
												v[i][j-1]=0;
											}
											else{
												u[i][j]=u[i][j]+dt*dp/h;
												u[i-1][j]=u[i-1][j]-dt*dp/h;
												v[i][j]=v[i][j]+dt*dp/h;
												v[i][j-1]=0;
												}
											}
										else {
											if(i==1){
												u[i][j]=u[i][j]+dt*dp/h;
												u[i-1][j]=u[1][j];
												v[i][j]=v[i][j]+dt*dp/h;
												v[i][j-1]=v[i][j-1]-dt*dp/h;
												}
											if(i==Xn-1){
												u[i][j]=u[i-1][j];
												u[i-1][j]=u[i-1][j]-dt*dp/h;
												v[i][j]=v[i-1][j];
												v[i][j-1]=v[i-1][j-1];
											}
											else{
												u[i][j]=u[i][j]+dt*dp/h;
												u[i-1][j]=u[i-1][j]-dt*dp/h;
												v[i][j]=v[i][j]+dt*dp/h;
												v[i][j-1]=v[i][j-1]-dt*dp/h;
												}
										}
										}
									}
								 System.out.printf("Dmax=%1.8f,t=%1.8f\n",sa,tcount*dt);
									}while(sa>J2);

						//データの保存
						for(int rt=0;rt<5;rt++){
						if(tcount==memot[rt]){
							for(int i=0;i<Xn;i++){
								for(int j=0;j<Yn;j++){
									memou[i][j][rt]=u[i][j];
									memov[i][j][rt]=v[i][j];
									memop[i][j][rt]=p[i][j];
									memoΦ[i][j][rt]=Φ[i][j];
									memoabsE[i][j][rt]=absE[i][j];
									memoEx[i][j][rt]=Ex[i][j];
									memoEy[i][j][rt]=Ey[i][j];
									memofx[i][j][rt]=fx[i][j];
									memofy[i][j][rt]=fy[i][j];
									memoq[i][j][rt]=q[i][j];
									memojx[i][j][rt]=jx[i][j];
									memojy[i][j][rt]=jy[i][j];
									memoAx[i][j][rt]=jx[i][j]*h;
									memoAy[i][j][rt]=jy[i][j]*h;
								}
							}
						}
						}

					}while(tcount*dt<=t);
			//ファイルへの出力
					try {
						FileOutputStream fosU[] = new FileOutputStream[5];
						FileOutputStream fosP[] = new FileOutputStream[5];
						FileOutputStream fosΦ[] = new FileOutputStream[5];
						FileOutputStream fosE[] = new FileOutputStream[5];
						FileOutputStream fosabsE[] = new FileOutputStream[5];
						FileOutputStream fosF[] = new FileOutputStream[5];
						FileOutputStream fosQ[] = new FileOutputStream[5];
						FileOutputStream fosJ[] = new FileOutputStream[5];
						FileOutputStream fosA[] = new FileOutputStream[5];
						FileOutputStream fosJt;

						BufferedOutputStream bosU[] = new BufferedOutputStream[5];
						BufferedOutputStream bosP[] = new BufferedOutputStream[5];
						BufferedOutputStream bosΦ[]= new BufferedOutputStream[5];
						BufferedOutputStream bosE[] = new BufferedOutputStream[5];
						BufferedOutputStream bosabsE[] = new BufferedOutputStream[5];
						BufferedOutputStream bosF[] = new BufferedOutputStream[5];
						BufferedOutputStream bosQ[] = new BufferedOutputStream[5];
						BufferedOutputStream bosJ[] = new BufferedOutputStream[5];
						BufferedOutputStream bosA[] = new BufferedOutputStream[5];
						BufferedOutputStream bosJt;

						PrintStream outU[] = new PrintStream[5];
						PrintStream outP[] = new PrintStream[5];
						PrintStream outΦ[] = new PrintStream[5];
						PrintStream outE[] = new PrintStream[5];
						PrintStream outabsE[] = new PrintStream[5];
						PrintStream outF[] = new PrintStream[5];
						PrintStream outQ[] = new PrintStream[5];
						PrintStream outJ[] = new PrintStream[5];
						PrintStream outA[] = new PrintStream[5];
						PrintStream outJt;

						for(int rt=0;rt<5;rt++){
							fosU[rt] = new FileOutputStream("x,y,U(t="+(dt*memot[rt])+").csv");
							fosP[rt] = new FileOutputStream("x,y,p(t="+(dt*memot[rt])+").csv");
							fosΦ[rt] = new FileOutputStream("x,y,Φ(t="+(dt*memot[rt])+").csv");
							fosE[rt] = new FileOutputStream("x,y,E(t="+(dt*memot[rt])+").csv");
							fosabsE[rt] = new FileOutputStream("x,y,absE(t="+(dt*memot[rt])+").csv");
							fosF[rt] = new FileOutputStream("x,y,F(t="+(dt*memot[rt])+").csv");
							fosQ[rt] = new FileOutputStream("x,y,Q(t="+(dt*memot[rt])+").csv");
							fosJ[rt] = new FileOutputStream("x,y,J(t="+(dt*memot[rt])+").csv");
							fosA[rt] = new FileOutputStream("x,y,A(t="+(dt*memot[rt])+").csv");
							fosJt = new FileOutputStream("x,y,Jt(t="+(dt*memot[rt])+").csv");

							bosU[rt] = new BufferedOutputStream(fosU[rt]);
							bosP[rt] = new BufferedOutputStream(fosP[rt]);
							bosΦ[rt] = new BufferedOutputStream(fosΦ[rt]);
							bosE[rt] = new BufferedOutputStream(fosE[rt]);
							bosabsE[rt] = new BufferedOutputStream(fosabsE[rt]);
							bosF[rt] = new BufferedOutputStream(fosF[rt]);
							bosQ[rt] = new BufferedOutputStream(fosQ[rt]);
							bosJ[rt] = new BufferedOutputStream(fosJ[rt]);
							bosA[rt] = new BufferedOutputStream(fosA[rt]);
							bosJt = new BufferedOutputStream(fosJt);

							outU[rt] = new PrintStream(bosU[rt]);
							outP[rt] = new PrintStream(bosP[rt]);
							outΦ[rt] = new PrintStream(bosΦ[rt]);
							outE[rt] = new PrintStream(bosE[rt]);
							outabsE[rt] = new PrintStream(bosabsE[rt]);
							outF[rt] = new PrintStream(bosF[rt]);
							outQ[rt] = new PrintStream(bosQ[rt]);
							outJ[rt] = new PrintStream(bosJ[rt]);
							outA[rt] = new PrintStream(bosA[rt]);
							outJt = new PrintStream(bosJt);

							outU[rt].printf("\n\n\n");
							outP[rt].printf("\n\n\n");
							outΦ[rt].printf("\n\n\n");
							outE[rt].printf("\n\n\n");
							outabsE[rt].printf("\n\n\n");
							outF[rt].printf("\n\n\n");
							outQ[rt].printf("\n\n\n");
							outJ[rt].printf("\n\n\n");
							outA[rt].printf("\n\n\n");
							outJt.printf("\n\n\n");

							for(int i=0;i<Xn;i++){
							for(int j=0;j<Yn;j++){
								outU[rt].printf("%3.8f,%3.8f,0,%3.8f,%3.8f,0\n",x[i],y[j],memou[i][j][rt],memov[i][j][rt]);
								outP[rt].printf("%3.8f,%3.8f,0,%3.8f\n",x[i],y[j],memop[i][j][rt]);
								outΦ[rt].printf("%3.8f,%3.8f,0,%3.8f\n",x[i],y[j],memoΦ[i][j][rt]);
								outE[rt].printf("%3.8f,%3.8f,0,%3.8f,%3.8f,0\n",x[i],y[j],memoEx[i][j][rt],memoEy[i][j][rt]);
								outabsE[rt].printf("%3.8f,%3.8f,0,%3.8f\n",x[i],y[j],memoabsE[i][j][rt]);
								outF[rt].printf("%3.8f,%3.8f,0,%3.8f,%3.8f,0\n",x[i],y[j],memofx[i][j][rt],memofy[i][j][rt]);
								outQ[rt].printf("%3.8f,%3.8f,0,%3.8f\n",x[i],y[j],memoq[i][j][rt]);
								outJ[rt].printf("%3.8f,%3.8f,0,%3.8f,%3.8f,0\n",x[i],y[j],memojx[i][j][rt],memojy[i][j][rt]);
								outA[rt].printf("%3.8f,%3.8f,0,%3.8f,%3.8f,0\n",x[i],y[j],memoAx[i][j][rt],memoAy[i][j][rt]);
								}
						}

						outU[rt].close();
						outP[rt].close();
						outΦ[rt].close();
						outE[rt].close();
						outabsE[rt].close();
						outF[rt].close();
						outQ[rt].close();
						outJ[rt].close();
						outA[rt].close();
						outJt.close();
						}
					}catch (FileNotFoundException e) {
						e.printStackTrace();
					}
		}
	}
