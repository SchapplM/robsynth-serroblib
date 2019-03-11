% Calculate time derivative of joint inertia matrix for
% S4RPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRP2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP2_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_inertiaDJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP2_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP2_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP2_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:30:42
% EndTime: 2019-03-08 18:30:42
% DurationCPUTime: 0.28s
% Computational Cost: add. (487->87), mult. (1068->124), div. (0->0), fcn. (868->4), ass. (0->54)
t39 = cos(qJ(3));
t40 = cos(qJ(1));
t37 = sin(qJ(3));
t38 = sin(qJ(1));
t59 = t38 * t37;
t45 = t40 * t39 + t59;
t43 = t45 * qJD(3);
t15 = t45 * qJD(1) - t43;
t54 = qJD(1) * t38;
t57 = t40 * t37;
t46 = -t38 * t39 + t57;
t53 = qJD(1) * t40;
t64 = t46 * qJD(3) - t37 * t53;
t16 = -t39 * t54 - t64;
t65 = -t16 * rSges(5,1) - t15 * rSges(5,2) + t64 * pkin(3);
t63 = -pkin(1) - pkin(2);
t62 = -rSges(3,1) - pkin(1);
t30 = t39 * pkin(3) + pkin(2);
t61 = -pkin(1) - t30;
t60 = -pkin(2) + t30;
t56 = qJ(2) * t53 + qJD(2) * t38;
t33 = t38 * qJ(2);
t55 = t40 * pkin(1) + t33;
t52 = t63 * t38;
t50 = t61 * t38;
t49 = qJD(1) * t60;
t8 = -t16 * rSges(4,1) - t15 * rSges(4,2);
t7 = t15 * rSges(4,1) - t16 * rSges(4,2);
t17 = -rSges(4,1) * t46 - rSges(4,2) * t45;
t18 = -rSges(4,1) * t45 + rSges(4,2) * t46;
t47 = t15 * rSges(5,1) - t16 * rSges(5,2);
t44 = t40 * rSges(3,3) + t62 * t38;
t42 = -rSges(5,1) * t46 - rSges(5,2) * t45 - pkin(3) * t57;
t41 = -rSges(5,1) * t45 + rSges(5,2) * t46 - pkin(3) * t59 - t40 * t30;
t35 = t40 * pkin(2);
t34 = t40 * qJ(2);
t32 = qJD(2) * t40;
t22 = t40 * rSges(3,1) + t38 * rSges(3,3) + t55;
t21 = t34 + t44;
t20 = t32 + (t62 * t40 + (-rSges(3,3) - qJ(2)) * t38) * qJD(1);
t19 = t44 * qJD(1) + t56;
t14 = -t18 + t35 + t55;
t13 = t34 + t52 - t17;
t12 = t35 + t41;
t11 = t60 * t38 + t42;
t10 = -t41 + t55;
t9 = t34 + t50 - t42;
t6 = t32 + (t63 * t40 - t33) * qJD(1) - t7;
t5 = qJD(1) * t52 + t56 - t8;
t4 = t40 * t49 + (t37 * t54 - t43) * pkin(3) + t47;
t3 = t38 * t49 + t65;
t2 = t32 + pkin(3) * t43 + (t61 * t40 + (-pkin(3) * t37 - qJ(2)) * t38) * qJD(1) - t47;
t1 = qJD(1) * t50 + t56 - t65;
t23 = [0.2e1 * m(3) * (t22 * t19 + t21 * t20) + 0.2e1 * m(4) * (t13 * t6 + t14 * t5) + 0.2e1 * m(5) * (t10 * t1 + t9 * t2); m(3) * (-t40 * t19 + t38 * t20 + (t21 * t40 + t22 * t38) * qJD(1)) + m(4) * (t38 * t6 - t40 * t5 + (t13 * t40 + t14 * t38) * qJD(1)) + m(5) * (-t40 * t1 + t38 * t2 + (t10 * t38 + t40 * t9) * qJD(1)); 0; m(4) * (t7 * t13 + t8 * t14 + t17 * t6 + t18 * t5) + m(5) * (t12 * t1 + t3 * t10 + t11 * t2 + t4 * t9); m(4) * (t7 * t38 - t8 * t40 + (t17 * t40 + t18 * t38) * qJD(1)) + m(5) * (-t3 * t40 + t4 * t38 + (t11 * t40 + t12 * t38) * qJD(1)); 0.2e1 * m(4) * (t17 * t7 + t18 * t8) + 0.2e1 * m(5) * (t11 * t4 + t12 * t3); 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t23(1) t23(2) t23(4) t23(7); t23(2) t23(3) t23(5) t23(8); t23(4) t23(5) t23(6) t23(9); t23(7) t23(8) t23(9) t23(10);];
Mq  = res;
