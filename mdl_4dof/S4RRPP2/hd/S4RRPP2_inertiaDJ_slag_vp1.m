% Calculate time derivative of joint inertia matrix for
% S4RRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:52
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4RRPP2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP2_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP2_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP2_inertiaDJ_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP2_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP2_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPP2_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:52:25
% EndTime: 2018-11-14 13:52:25
% DurationCPUTime: 0.20s
% Computational Cost: add. (701->75), mult. (576->96), div. (0->0), fcn. (328->4), ass. (0->47)
t61 = rSges(5,1) + pkin(3);
t47 = sin(qJ(1));
t60 = t47 * pkin(1);
t59 = -rSges(4,1) - pkin(2);
t46 = qJ(1) + qJ(2);
t42 = sin(t46);
t58 = t42 * rSges(3,2);
t43 = cos(t46);
t45 = qJD(1) + qJD(2);
t57 = t43 * t45;
t32 = t43 * qJ(3);
t56 = qJD(3) * t42 + t45 * t32;
t55 = t43 * pkin(2) + t42 * qJ(3);
t54 = pkin(1) * qJD(1);
t53 = -pkin(2) - t61;
t52 = t47 * t54;
t48 = cos(qJ(1));
t51 = t48 * t54;
t50 = t59 * t42;
t24 = t43 * rSges(3,1) - t58;
t16 = t43 * rSges(4,1) + t42 * rSges(4,3) + t55;
t49 = t53 * t42;
t20 = -rSges(3,1) * t57 + t45 * t58;
t12 = t42 * rSges(5,2) + t61 * t43 + t55;
t23 = -t42 * rSges(3,1) - t43 * rSges(3,2);
t15 = t43 * rSges(4,3) + t32 + t50;
t19 = t23 * t45;
t11 = t43 * rSges(5,2) + t32 + t49;
t7 = rSges(4,3) * t57 + t45 * t50 + t56;
t3 = rSges(5,2) * t57 + t45 * t49 + t56;
t30 = qJD(3) * t43;
t8 = t30 + (t59 * t43 + (-rSges(4,3) - qJ(3)) * t42) * t45;
t4 = t30 + ((-rSges(5,2) - qJ(3)) * t42 + t53 * t43) * t45;
t44 = t48 * pkin(1);
t22 = t24 + t44;
t21 = t23 - t60;
t18 = t20 - t51;
t17 = t19 - t52;
t14 = t44 + t16;
t13 = t15 - t60;
t10 = t44 + t12;
t9 = t11 - t60;
t6 = t8 - t51;
t5 = t7 - t52;
t2 = t4 - t51;
t1 = t3 - t52;
t25 = [0.2e1 * m(3) * (t22 * t17 + t21 * t18) + 0.2e1 * m(4) * (t13 * t6 + t14 * t5) + 0.2e1 * m(5) * (t10 * t1 + t9 * t2); m(3) * (t24 * t17 + t23 * t18 + t19 * t22 + t20 * t21) + m(4) * (t8 * t13 + t7 * t14 + t15 * t6 + t16 * t5) + m(5) * (t12 * t1 + t3 * t10 + t11 * t2 + t4 * t9); 0.2e1 * m(3) * (t24 * t19 + t23 * t20) + 0.2e1 * m(4) * (t15 * t8 + t16 * t7) + 0.2e1 * m(5) * (t11 * t4 + t12 * t3); m(4) * ((t13 * t45 - t5) * t43 + (t14 * t45 + t6) * t42) + m(5) * ((t45 * t9 - t1) * t43 + (t10 * t45 + t2) * t42); m(4) * ((t15 * t45 - t7) * t43 + (t16 * t45 + t8) * t42) + m(5) * ((t11 * t45 - t3) * t43 + (t12 * t45 + t4) * t42); 0; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t25(1) t25(2) t25(4) t25(7); t25(2) t25(3) t25(5) t25(8); t25(4) t25(5) t25(6) t25(9); t25(7) t25(8) t25(9) t25(10);];
Mq  = res;
