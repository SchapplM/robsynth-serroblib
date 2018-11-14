% Calculate time derivative of joint inertia matrix for
% S4RPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2018-11-14 13:48
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4RPPR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR2_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR2_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR2_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR2_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:47:26
% EndTime: 2018-11-14 13:47:26
% DurationCPUTime: 0.18s
% Computational Cost: add. (439->73), mult. (646->105), div. (0->0), fcn. (518->6), ass. (0->48)
t38 = pkin(6) + qJ(4);
t31 = sin(t38);
t42 = cos(qJ(1));
t32 = cos(t38);
t41 = sin(qJ(1));
t55 = t41 * t32;
t20 = -t42 * t31 + t55;
t58 = -pkin(1) - pkin(2);
t57 = -rSges(3,1) - pkin(1);
t40 = cos(pkin(6));
t29 = t40 * pkin(3) + pkin(2);
t56 = -pkin(1) - t29;
t39 = sin(pkin(6));
t54 = t41 * t39;
t52 = t42 * t39;
t49 = qJD(1) * t42;
t51 = qJ(2) * t49 + qJD(2) * t41;
t35 = t41 * qJ(2);
t50 = t42 * pkin(1) + t35;
t48 = t58 * t41;
t47 = t41 * t31 + t42 * t32;
t11 = (qJD(1) - qJD(4)) * t47;
t12 = -qJD(1) * t55 + t20 * qJD(4) + t31 * t49;
t4 = -t12 * rSges(5,1) - t11 * rSges(5,2);
t3 = t11 * rSges(5,1) - t12 * rSges(5,2);
t13 = t20 * rSges(5,1) - rSges(5,2) * t47;
t14 = -rSges(5,1) * t47 - t20 * rSges(5,2);
t46 = -t41 * t40 + t52;
t45 = t42 * t40 + t54;
t44 = t42 * rSges(3,3) + t57 * t41;
t43 = pkin(3) * t52 + t56 * t41;
t36 = t42 * qJ(2);
t34 = qJD(2) * t42;
t22 = t46 * qJD(1);
t21 = t45 * qJD(1);
t18 = t42 * rSges(3,1) + t41 * rSges(3,3) + t50;
t17 = t36 + t44;
t16 = t34 + (t57 * t42 + (-rSges(3,3) - qJ(2)) * t41) * qJD(1);
t15 = t44 * qJD(1) + t51;
t10 = rSges(4,1) * t45 - rSges(4,2) * t46 + t42 * pkin(2) + t50;
t9 = rSges(4,1) * t46 + rSges(4,2) * t45 + t36 + t48;
t8 = -t21 * rSges(4,1) + t22 * rSges(4,2) + t34 + (t58 * t42 - t35) * qJD(1);
t7 = t22 * rSges(4,1) + t21 * rSges(4,2) + qJD(1) * t48 + t51;
t6 = pkin(3) * t54 + t42 * t29 - t14 + t50;
t5 = -t13 + t36 + t43;
t2 = t34 + (t56 * t42 + (-pkin(3) * t39 - qJ(2)) * t41) * qJD(1) - t3;
t1 = t43 * qJD(1) - t4 + t51;
t19 = [0.2e1 * m(3) * (t18 * t15 + t17 * t16) + 0.2e1 * m(4) * (t10 * t7 + t9 * t8) + 0.2e1 * m(5) * (t6 * t1 + t5 * t2); m(3) * (-t42 * t15 + t41 * t16 + (t17 * t42 + t18 * t41) * qJD(1)) + m(4) * (t41 * t8 - t42 * t7 + (t10 * t41 + t42 * t9) * qJD(1)) + m(5) * (-t42 * t1 + t41 * t2 + (t41 * t6 + t42 * t5) * qJD(1)); 0; 0; 0; 0; m(5) * (t14 * t1 + t13 * t2 + t3 * t5 + t4 * t6); m(5) * (t3 * t41 - t4 * t42 + (t13 * t42 + t14 * t41) * qJD(1)); 0; 0.2e1 * m(5) * (t13 * t3 + t14 * t4);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t19(1) t19(2) t19(4) t19(7); t19(2) t19(3) t19(5) t19(8); t19(4) t19(5) t19(6) t19(9); t19(7) t19(8) t19(9) t19(10);];
Mq  = res;
