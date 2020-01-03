% Calculate time derivative of joint inertia matrix for
% S5RPPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPPR6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR6_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR6_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR6_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR6_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:37
% EndTime: 2019-12-31 17:47:40
% DurationCPUTime: 1.08s
% Computational Cost: add. (2594->280), mult. (7210->402), div. (0->0), fcn. (7398->8), ass. (0->131)
t156 = 2 * m(6);
t155 = rSges(6,3) + pkin(6);
t154 = -pkin(2) - qJ(4);
t111 = sin(qJ(1));
t113 = cos(qJ(1));
t141 = qJD(1) * t113;
t153 = qJ(2) * t141 + qJD(2) * t111;
t106 = sin(pkin(8));
t152 = t106 * t111;
t107 = sin(pkin(7));
t151 = t107 * t113;
t108 = cos(pkin(8));
t150 = t108 * t113;
t109 = cos(pkin(7));
t149 = t109 * t108;
t110 = sin(qJ(5));
t148 = t109 * t110;
t112 = cos(qJ(5));
t147 = t109 * t112;
t146 = t109 * t113;
t145 = t111 * t108;
t144 = t113 * pkin(1) + t111 * qJ(2);
t102 = t113 * qJ(2);
t143 = t113 * pkin(3) + t102;
t142 = qJD(1) * t111;
t140 = qJD(3) * t107;
t139 = qJD(4) * t109;
t136 = t109 * t142;
t87 = t106 * t151 + t145;
t71 = t110 * t146 + t87 * t112;
t89 = t107 * t152 - t150;
t82 = t89 * qJD(1);
t49 = -t71 * qJD(5) + t82 * t110 - t112 * t136;
t70 = -t87 * t110 + t112 * t146;
t50 = t70 * qJD(5) - t110 * t136 - t82 * t112;
t88 = t106 * t113 + t107 * t145;
t81 = t88 * qJD(1);
t26 = t50 * rSges(6,1) + t49 * rSges(6,2) + t81 * rSges(6,3);
t137 = t107 * t150;
t86 = -t137 + t152;
t39 = t71 * rSges(6,1) + t70 * rSges(6,2) + t86 * rSges(6,3);
t138 = t113 * t140 + t153;
t135 = t109 * t141;
t134 = (-pkin(3) - qJ(2)) * t111;
t133 = -qJ(3) * t107 - pkin(1);
t132 = pkin(2) * t146 + qJ(3) * t151 + t144;
t124 = t106 * t147 - t107 * t110;
t84 = t106 * t148 + t107 * t112;
t57 = -Icges(6,4) * t124 + Icges(6,2) * t84 + Icges(6,6) * t149;
t58 = -Icges(6,1) * t124 + Icges(6,4) * t84 + Icges(6,5) * t149;
t77 = t84 * qJD(5);
t78 = t124 * qJD(5);
t62 = Icges(6,5) * t77 + Icges(6,6) * t78;
t63 = Icges(6,4) * t77 + Icges(6,2) * t78;
t64 = Icges(6,1) * t77 + Icges(6,4) * t78;
t131 = -t124 * t64 + t62 * t149 + t78 * t57 + t77 * t58 + t84 * t63;
t130 = pkin(3) * t141 + t113 * t139 + t138;
t73 = t111 * t148 + t112 * t89;
t80 = t87 * qJD(1);
t47 = -t73 * qJD(5) - t80 * t110 + t112 * t135;
t72 = -t110 * t89 + t111 * t147;
t48 = t72 * qJD(5) + t110 * t135 + t80 * t112;
t129 = -t48 * rSges(6,1) - t47 * rSges(6,2);
t128 = -rSges(6,1) * t73 - rSges(6,2) * t72;
t127 = rSges(3,1) * t109 - rSges(3,2) * t107;
t126 = t111 * pkin(3) + qJ(4) * t146 + t132;
t125 = -pkin(1) - t127;
t123 = t154 * t109 + t133;
t122 = (-rSges(5,3) + t154) * t109 + t133;
t121 = -pkin(1) + (rSges(4,2) - pkin(2)) * t109 + (-rSges(4,3) - qJ(3)) * t107;
t100 = qJD(2) * t113;
t120 = t100 + (-t139 - t140) * t111;
t119 = t123 * t111;
t118 = t122 * t111;
t117 = rSges(3,3) * t113 + t125 * t111;
t79 = -qJD(1) * t137 + t106 * t142;
t25 = t79 * rSges(6,3) - t129;
t59 = -rSges(6,1) * t124 + rSges(6,2) * t84 + rSges(6,3) * t149;
t65 = rSges(6,1) * t77 + rSges(6,2) * t78;
t11 = -t25 * t149 + t59 * t79 - t65 * t88;
t12 = t26 * t149 - t59 * t81 - t65 * t86;
t29 = t39 * t149 - t59 * t86;
t40 = -rSges(6,3) * t88 - t128;
t30 = -t40 * t149 - t59 * t88;
t116 = t11 * t113 + t111 * t12 + (-t111 * t30 + t113 * t29) * qJD(1);
t115 = rSges(4,1) * t113 + t121 * t111;
t15 = -pkin(4) * t82 + pkin(6) * t81 + qJD(1) * t119 + t130 + t26;
t16 = -t80 * pkin(4) - t155 * t79 + (t123 * t113 + t134) * qJD(1) + t120 + t129;
t27 = -pkin(4) * t89 + t155 * t88 + t119 + t128 + t143;
t28 = pkin(4) * t87 + pkin(6) * t86 + t126 + t39;
t31 = -rSges(5,1) * t82 - rSges(5,2) * t81 + qJD(1) * t118 + t130;
t32 = -t80 * rSges(5,1) + t79 * rSges(5,2) + (t122 * t113 + t134) * qJD(1) + t120;
t41 = -rSges(5,1) * t89 - rSges(5,2) * t88 + t118 + t143;
t42 = t87 * rSges(5,1) - t86 * rSges(5,2) + rSges(5,3) * t146 + t126;
t114 = m(6) * (t111 * t15 + t113 * t16 + t28 * t141 - t27 * t142) / 0.2e1 + m(5) * (t111 * t31 + t113 * t32 + t42 * t141 - t41 * t142) / 0.2e1;
t75 = t111 * rSges(3,3) + t127 * t113 + t144;
t74 = t102 + t117;
t67 = t100 + ((-rSges(3,3) - qJ(2)) * t111 + t125 * t113) * qJD(1);
t66 = t117 * qJD(1) + t153;
t61 = t111 * rSges(4,1) + (-rSges(4,2) * t109 + rSges(4,3) * t107) * t113 + t132;
t60 = t102 + t115;
t56 = -Icges(6,5) * t124 + Icges(6,6) * t84 + Icges(6,3) * t149;
t46 = -t111 * t140 + t100 + ((-rSges(4,1) - qJ(2)) * t111 + t121 * t113) * qJD(1);
t45 = t115 * qJD(1) + t138;
t38 = Icges(6,1) * t73 + Icges(6,4) * t72 - Icges(6,5) * t88;
t37 = Icges(6,1) * t71 + Icges(6,4) * t70 + Icges(6,5) * t86;
t36 = Icges(6,4) * t73 + Icges(6,2) * t72 - Icges(6,6) * t88;
t35 = Icges(6,4) * t71 + Icges(6,2) * t70 + Icges(6,6) * t86;
t34 = Icges(6,5) * t73 + Icges(6,6) * t72 - Icges(6,3) * t88;
t33 = Icges(6,5) * t71 + Icges(6,6) * t70 + Icges(6,3) * t86;
t24 = Icges(6,1) * t50 + Icges(6,4) * t49 + Icges(6,5) * t81;
t23 = Icges(6,1) * t48 + Icges(6,4) * t47 + Icges(6,5) * t79;
t22 = Icges(6,4) * t50 + Icges(6,2) * t49 + Icges(6,6) * t81;
t21 = Icges(6,4) * t48 + Icges(6,2) * t47 + Icges(6,6) * t79;
t20 = Icges(6,5) * t50 + Icges(6,6) * t49 + Icges(6,3) * t81;
t19 = Icges(6,5) * t48 + Icges(6,6) * t47 + Icges(6,3) * t79;
t18 = -t56 * t88 + t57 * t72 + t58 * t73;
t17 = t56 * t86 + t57 * t70 + t58 * t71;
t14 = -t124 * t38 + t34 * t149 + t36 * t84;
t13 = -t124 * t37 + t33 * t149 + t35 * t84;
t10 = t131 * t149;
t9 = -t34 * t88 + t36 * t72 + t38 * t73;
t8 = -t33 * t88 + t35 * t72 + t37 * t73;
t7 = t34 * t86 + t36 * t70 + t38 * t71;
t6 = t33 * t86 + t35 * t70 + t37 * t71;
t5 = t49 * t57 + t50 * t58 + t56 * t81 + t62 * t86 + t63 * t70 + t64 * t71;
t4 = t47 * t57 + t48 * t58 + t56 * t79 - t62 * t88 + t63 * t72 + t64 * t73;
t3 = t25 * t86 + t26 * t88 - t39 * t79 + t40 * t81;
t2 = -t124 * t23 + t19 * t149 + t21 * t84 + t36 * t78 + t38 * t77;
t1 = -t124 * t24 + t20 * t149 + t22 * t84 + t35 * t78 + t37 * t77;
t43 = [(t15 * t28 + t16 * t27) * t156 + 0.2e1 * m(4) * (t45 * t61 + t46 * t60) + 0.2e1 * m(5) * (t31 * t42 + t32 * t41) + 0.2e1 * m(3) * (t66 * t75 + t67 * t74) + t131; m(6) * (t111 * t16 - t113 * t15 + (t111 * t28 + t113 * t27) * qJD(1)) + m(4) * (t111 * t46 - t113 * t45 + (t111 * t61 + t113 * t60) * qJD(1)) + m(5) * (t111 * t32 - t113 * t31 + (t111 * t42 + t113 * t41) * qJD(1)) + m(3) * (t111 * t67 - t113 * t66 + (t111 * t75 + t113 * t74) * qJD(1)); 0; 0.2e1 * (m(4) * (t111 * t45 + t113 * t46 + t61 * t141 - t60 * t142) / 0.2e1 + t114) * t107; 0; 0; 0.2e1 * t114 * t109; 0; 0; 0; t10 + m(6) * (t11 * t27 + t12 * t28 + t15 * t29 + t16 * t30) + (-t2 / 0.2e1 - t4 / 0.2e1) * t88 + (t1 / 0.2e1 + t5 / 0.2e1) * t86 + (t13 / 0.2e1 + t17 / 0.2e1) * t81 + (t14 / 0.2e1 + t18 / 0.2e1) * t79; m(6) * (t11 * t111 - t113 * t12 + (t111 * t29 + t113 * t30) * qJD(1)); m(6) * (t116 * t107 - t3 * t109); m(6) * (t3 * t107 + t116 * t109); (t30 * t11 + t29 * t12 + (t39 * t88 + t40 * t86) * t3) * t156 + t81 * (t17 * t149 + t6 * t86 - t7 * t88) + t86 * ((t20 * t86 + t22 * t70 + t24 * t71 + t33 * t81 + t49 * t35 + t50 * t37) * t86 + t6 * t81 - (t19 * t86 + t21 * t70 + t23 * t71 + t34 * t81 + t49 * t36 + t50 * t38) * t88 + t7 * t79 + t5 * t149) + t79 * (t18 * t149 + t8 * t86 - t88 * t9) - t88 * ((-t20 * t88 + t22 * t72 + t24 * t73 + t33 * t79 + t47 * t35 + t48 * t37) * t86 + t8 * t81 - (-t19 * t88 + t21 * t72 + t23 * t73 + t34 * t79 + t47 * t36 + t48 * t38) * t88 + t9 * t79 + t4 * t149) + (t1 * t86 + t13 * t81 + t14 * t79 - t2 * t88 + t10) * t149;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t43(1), t43(2), t43(4), t43(7), t43(11); t43(2), t43(3), t43(5), t43(8), t43(12); t43(4), t43(5), t43(6), t43(9), t43(13); t43(7), t43(8), t43(9), t43(10), t43(14); t43(11), t43(12), t43(13), t43(14), t43(15);];
Mq = res;
