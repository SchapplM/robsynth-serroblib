% Calculate time derivative of joint inertia matrix for
% S4RPPR7
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPPR7_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR7_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR7_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR7_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR7_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:35
% EndTime: 2019-12-31 16:41:37
% DurationCPUTime: 0.92s
% Computational Cost: add. (999->149), mult. (1506->226), div. (0->0), fcn. (1108->6), ass. (0->90)
t54 = pkin(6) + qJ(4);
t45 = sin(t54);
t46 = cos(t54);
t78 = rSges(5,1) * t45 + rSges(5,2) * t46;
t61 = cos(qJ(1));
t60 = sin(qJ(1));
t98 = Icges(5,4) * t45;
t70 = Icges(5,2) * t46 + t98;
t23 = Icges(5,6) * t61 + t60 * t70;
t97 = Icges(5,4) * t46;
t72 = Icges(5,1) * t45 + t97;
t25 = Icges(5,5) * t61 + t60 * t72;
t75 = t23 * t46 + t25 * t45;
t124 = t61 * t75;
t123 = t61 * t78;
t51 = t61 * qJ(2);
t57 = sin(pkin(6));
t115 = pkin(3) * t57;
t63 = t78 + t115;
t59 = -pkin(5) - qJ(3);
t86 = -rSges(5,3) - pkin(1) + t59;
t80 = t86 * t60;
t15 = t61 * t63 + t51 + t80;
t52 = t61 * rSges(5,3);
t27 = t78 * t60 + t52;
t99 = t61 * pkin(1) + t60 * qJ(2);
t16 = t115 * t60 - t59 * t61 + t27 + t99;
t122 = -t15 * t60 + t16 * t61;
t109 = rSges(5,2) * t45;
t91 = qJD(1) * t61;
t100 = qJ(2) * t91 + qJD(2) * t60;
t83 = qJD(3) * t61 + t100;
t110 = rSges(5,1) * t46;
t88 = qJD(4) * t60;
t84 = t88 * t110 + t78 * t91;
t5 = -t88 * t109 + (t61 * t115 + t80) * qJD(1) + t83 + t84;
t38 = -t109 + t110;
t49 = qJD(2) * t61;
t81 = -qJD(3) * t60 + t49;
t87 = qJD(4) * t61;
t6 = t38 * t87 + (t86 * t61 + (-qJ(2) - t63) * t60) * qJD(1) + t81;
t121 = -t5 * t61 + t60 * t6;
t68 = Icges(5,5) * t45 + Icges(5,6) * t46;
t120 = -Icges(5,3) * t60 + t61 * t68;
t119 = -Icges(5,6) * t60 + t61 * t70;
t118 = -Icges(5,5) * t60 + t61 * t72;
t117 = 2 * m(5);
t55 = t60 ^ 2;
t56 = t61 ^ 2;
t116 = m(5) * t38;
t112 = rSges(3,2) - pkin(1);
t105 = t23 * t45;
t104 = t119 * t45;
t103 = t25 * t46;
t102 = t118 * t46;
t101 = t60 * rSges(5,3);
t93 = rSges(4,3) + qJ(3);
t21 = Icges(5,3) * t61 + t60 * t68;
t92 = qJD(1) * t21;
t90 = qJD(4) * t45;
t89 = qJD(4) * t46;
t85 = -pkin(1) - t93;
t34 = t78 * qJD(4);
t82 = (t55 + t56) * t34;
t79 = rSges(4,1) * t57 + rSges(4,2) * cos(pkin(6));
t74 = -t118 * t45 - t119 * t46;
t73 = Icges(5,1) * t46 - t98;
t71 = -Icges(5,2) * t45 + t97;
t69 = Icges(5,5) * t46 - Icges(5,6) * t45;
t67 = rSges(3,3) * t61 + t112 * t60;
t66 = t74 * t60;
t65 = qJD(4) * t73;
t64 = qJD(4) * t71;
t62 = t60 * t85 + t61 * t79;
t30 = -rSges(3,2) * t61 + t60 * rSges(3,3) + t99;
t29 = t51 + t67;
t28 = t101 - t123;
t20 = t49 + (t112 * t61 + (-rSges(3,3) - qJ(2)) * t60) * qJD(1);
t19 = qJD(1) * t67 + t100;
t18 = t60 * t79 + t61 * t93 + t99;
t17 = t51 + t62;
t10 = t120 * qJD(1) + t69 * t88;
t9 = -t69 * t87 + t92;
t8 = (t85 * t61 + (-qJ(2) - t79) * t60) * qJD(1) + t81;
t7 = qJD(1) * t62 + t83;
t4 = -t120 * t60 - t61 * t74;
t3 = t60 * t21 - t124;
t2 = -t120 * t61 + t66;
t1 = t21 * t61 + t75 * t60;
t11 = [0.2e1 * m(3) * (t19 * t30 + t20 * t29) + 0.2e1 * m(4) * (t17 * t8 + t18 * t7) - t45 * t65 - t72 * t89 - t46 * t64 + t70 * t90 + (t15 * t6 + t16 * t5) * t117; m(3) * (-t19 * t61 + t60 * t20 + (t29 * t61 + t30 * t60) * qJD(1)) + m(4) * (t60 * t8 - t61 * t7 + (t17 * t61 + t18 * t60) * qJD(1)) + m(5) * ((t15 * t61 + t16 * t60) * qJD(1) + t121); 0; m(4) * (t60 * t7 + t61 * t8 + (-t17 * t60 + t18 * t61) * qJD(1)) + m(5) * (t122 * qJD(1) + t60 * t5 + t6 * t61); 0; 0; (-qJD(4) * t75 - (t119 * qJD(1) + t60 * t64) * t45 + (t118 * qJD(1) + t60 * t65) * t46) * t61 / 0.2e1 + (-qJD(4) * t74 - (qJD(1) * t23 - t71 * t87) * t45 + (qJD(1) * t25 - t73 * t87) * t46) * t60 / 0.2e1 + m(5) * (t121 * t38 + t122 * t34) - (t56 / 0.2e1 + t55 / 0.2e1) * t68 * qJD(4) + ((t105 / 0.2e1 - t103 / 0.2e1 + t16 * t116) * t60 + (t104 / 0.2e1 - t102 / 0.2e1 + t15 * t116) * t61) * qJD(1); -m(5) * t82; 0; ((-t60 * t27 + t28 * t61) * (-t60 * t84 + (t109 * t55 - t38 * t56) * qJD(4) + ((-t27 + t52) * t61 + (t101 - t28 + t123) * t60) * qJD(1)) - t38 * t82) * t117 - qJD(1) * t60 * (t1 * t61 + t2 * t60) + t61 * ((t61 * t10 + (t2 + t124) * qJD(1)) * t61 + (-t1 * qJD(1) + (-t118 * t89 + t119 * t90) * t60 + (t9 + (t103 - t105) * qJD(4) + (-t21 + t74) * qJD(1)) * t61) * t60) + (t3 * t61 + t4 * t60) * t91 + t60 * ((t60 * t9 + (-t3 + t66) * qJD(1)) * t60 + (t4 * qJD(1) + (t23 * t90 - t25 * t89 + t92) * t61 + (t10 + (t102 - t104) * qJD(4) + t75 * qJD(1)) * t60) * t61);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t11(1), t11(2), t11(4), t11(7); t11(2), t11(3), t11(5), t11(8); t11(4), t11(5), t11(6), t11(9); t11(7), t11(8), t11(9), t11(10);];
Mq = res;
