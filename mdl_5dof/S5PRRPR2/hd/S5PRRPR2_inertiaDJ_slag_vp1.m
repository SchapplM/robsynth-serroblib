% Calculate time derivative of joint inertia matrix for
% S5PRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR2_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR2_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR2_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR2_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:18
% EndTime: 2019-12-05 16:17:23
% DurationCPUTime: 1.48s
% Computational Cost: add. (5477->248), mult. (4646->369), div. (0->0), fcn. (4218->8), ass. (0->124)
t116 = sin(qJ(5));
t117 = cos(qJ(5));
t114 = sin(pkin(9));
t115 = cos(pkin(9));
t145 = Icges(6,4) * t117;
t79 = -Icges(6,6) * t115 + (-Icges(6,2) * t116 + t145) * t114;
t146 = Icges(6,4) * t116;
t80 = -Icges(6,5) * t115 + (Icges(6,1) * t117 - t146) * t114;
t158 = -t116 * t80 - t117 * t79;
t157 = 2 * m(4);
t156 = 2 * m(5);
t155 = 2 * m(6);
t112 = pkin(8) + qJ(2);
t109 = sin(t112);
t154 = pkin(2) * t109;
t113 = qJD(2) + qJD(3);
t111 = qJ(3) + t112;
t107 = sin(t111);
t108 = cos(t111);
t138 = t115 * t116;
t121 = t107 * t138 + t108 * t117;
t137 = t115 * t117;
t77 = t107 * t116 + t108 * t137;
t57 = -t77 * qJD(5) + t121 * t113;
t75 = t107 * t137 - t108 * t116;
t76 = t107 * t117 - t108 * t138;
t58 = t76 * qJD(5) - t75 * t113;
t153 = t58 * rSges(6,1) + t57 * rSges(6,2);
t140 = t113 * t108;
t152 = qJ(4) * t140 + qJD(4) * t107;
t151 = pkin(2) * qJD(2);
t135 = qJD(5) * t114;
t83 = (-Icges(6,2) * t117 - t146) * t135;
t149 = t116 * t83;
t82 = (-Icges(6,5) * t116 - Icges(6,6) * t117) * t135;
t84 = (-Icges(6,1) * t116 - t145) * t135;
t127 = t114 * t117 * t84 - t115 * t82;
t26 = (t158 * qJD(5) - t149) * t114 + t127;
t147 = t26 * t115;
t100 = t107 * qJ(4);
t144 = t107 * t114;
t143 = t108 * t114;
t142 = t108 * t115;
t141 = t113 * t107;
t139 = t113 * t114;
t136 = t108 * pkin(3) + t100;
t54 = t77 * rSges(6,1) + t76 * rSges(6,2) + rSges(6,3) * t143;
t134 = t109 * t151;
t110 = cos(t112);
t133 = t110 * t151;
t132 = t107 * t139;
t131 = t108 * t139;
t128 = -rSges(5,1) * t115 - pkin(3);
t87 = t108 * rSges(4,1) - rSges(4,2) * t107;
t71 = -rSges(4,1) * t140 + rSges(4,2) * t141;
t59 = -t75 * qJD(5) + t76 * t113;
t60 = -t121 * qJD(5) + t77 * t113;
t124 = -rSges(6,1) * t60 - t59 * rSges(6,2);
t123 = -rSges(6,1) * t75 + rSges(6,2) * t121;
t122 = t128 * t107;
t86 = -rSges(4,1) * t107 - rSges(4,2) * t108;
t40 = pkin(4) * t142 + pkin(7) * t143 + t136 + t54;
t70 = t86 * t113;
t120 = -pkin(4) * t115 - pkin(3) + (-rSges(6,3) - pkin(7)) * t114;
t64 = rSges(5,1) * t142 - rSges(5,2) * t143 + t107 * rSges(5,3) + t136;
t101 = t108 * qJ(4);
t63 = rSges(5,2) * t144 + t108 * rSges(5,3) + t101 + t122;
t119 = t120 * t107;
t47 = Icges(6,5) * t75 - Icges(6,6) * t121 + Icges(6,3) * t144;
t49 = Icges(6,4) * t75 - Icges(6,2) * t121 + Icges(6,6) * t144;
t51 = Icges(6,1) * t75 - Icges(6,4) * t121 + Icges(6,5) * t144;
t16 = -t115 * t47 + (-t116 * t49 + t117 * t51) * t114;
t48 = Icges(6,5) * t77 + Icges(6,6) * t76 + Icges(6,3) * t143;
t50 = Icges(6,4) * t77 + Icges(6,2) * t76 + Icges(6,6) * t143;
t52 = Icges(6,1) * t77 + Icges(6,4) * t76 + Icges(6,5) * t143;
t17 = -t115 * t48 + (-t116 * t50 + t117 * t52) * t114;
t28 = Icges(6,5) * t60 + Icges(6,6) * t59 + Icges(6,3) * t131;
t30 = Icges(6,4) * t60 + Icges(6,2) * t59 + Icges(6,6) * t131;
t32 = Icges(6,1) * t60 + Icges(6,4) * t59 + Icges(6,5) * t131;
t3 = -t115 * t28 + (-t116 * t30 + t117 * t32 + (-t116 * t51 - t117 * t49) * qJD(5)) * t114;
t78 = -Icges(6,3) * t115 + (Icges(6,5) * t117 - Icges(6,6) * t116) * t114;
t35 = -t121 * t79 + t78 * t144 + t75 * t80;
t36 = t78 * t143 + t76 * t79 + t77 * t80;
t27 = Icges(6,5) * t58 + Icges(6,6) * t57 - Icges(6,3) * t132;
t29 = Icges(6,4) * t58 + Icges(6,2) * t57 - Icges(6,6) * t132;
t31 = Icges(6,1) * t58 + Icges(6,4) * t57 - Icges(6,5) * t132;
t4 = -t115 * t27 + (-t116 * t29 + t117 * t31 + (-t116 * t52 - t117 * t50) * qJD(5)) * t114;
t8 = t57 * t79 + t58 * t80 + t76 * t83 + t77 * t84 + (t108 * t82 - t78 * t141) * t114;
t9 = t59 * t79 + t60 * t80 - t121 * t83 + t75 * t84 + (t107 * t82 + t78 * t140) * t114;
t118 = -t147 + (t3 + t9) * t144 / 0.2e1 - (t36 + t17) * t132 / 0.2e1 + (t4 + t8 + (t16 + t35) * t113) * t143 / 0.2e1;
t45 = rSges(5,2) * t132 + rSges(5,3) * t140 + t113 * t122 + t152;
t99 = qJD(4) * t108;
t46 = rSges(5,2) * t131 + t99 + (t128 * t108 + (-rSges(5,3) - qJ(4)) * t107) * t113;
t22 = t113 * t119 + t152 + t153;
t39 = t101 + t119 + t123;
t23 = t99 + (t120 * t108 - t100) * t113 + t124;
t106 = pkin(2) * t110;
t85 = (-rSges(6,1) * t116 - rSges(6,2) * t117) * t135;
t81 = -t115 * rSges(6,3) + (rSges(6,1) * t117 - rSges(6,2) * t116) * t114;
t73 = t106 + t87;
t72 = t86 - t154;
t66 = t71 - t133;
t65 = t70 - t134;
t62 = t106 + t64;
t61 = t63 - t154;
t53 = rSges(6,3) * t144 - t123;
t44 = t46 - t133;
t43 = t45 - t134;
t42 = -t115 * t54 - t81 * t143;
t41 = t115 * t53 + t81 * t144;
t38 = t106 + t40;
t37 = t39 - t154;
t34 = rSges(6,3) * t131 - t124;
t33 = -rSges(6,3) * t132 + t153;
t21 = t115 * t34 + (t107 * t85 + t81 * t140) * t114;
t20 = -t115 * t33 + (-t108 * t85 + t81 * t141) * t114;
t19 = t23 - t133;
t18 = t22 - t134;
t13 = t48 * t143 + t50 * t76 + t52 * t77;
t12 = t47 * t143 + t49 * t76 + t51 * t77;
t11 = -t121 * t50 + t48 * t144 + t52 * t75;
t10 = -t121 * t49 + t47 * t144 + t51 * t75;
t5 = ((-t113 * t54 + t34) * t108 + (-t113 * t53 - t33) * t107) * t114;
t1 = [0; 0; -t114 * t149 + (t18 * t38 + t19 * t37) * t155 + (t43 * t62 + t44 * t61) * t156 + (t65 * t73 + t66 * t72) * t157 + t127 + t158 * t135; 0; m(6) * (t18 * t40 + t19 * t39 + t22 * t38 + t23 * t37) + m(5) * (t43 * t64 + t44 * t63 + t45 * t62 + t46 * t61) + m(4) * (t65 * t87 + t66 * t86 + t70 * t73 + t71 * t72) + t26; (t22 * t40 + t23 * t39) * t155 + (t45 * t64 + t46 * t63) * t156 + (t70 * t87 + t71 * t86) * t157 + t26; 0; m(6) * ((t113 * t37 - t18) * t108 + (t113 * t38 + t19) * t107) + m(5) * ((t113 * t61 - t43) * t108 + (t113 * t62 + t44) * t107); m(6) * ((t113 * t39 - t22) * t108 + (t113 * t40 + t23) * t107) + m(5) * ((t113 * t63 - t45) * t108 + (t113 * t64 + t46) * t107); 0; m(6) * t5; m(6) * (t18 * t42 + t19 * t41 + t20 * t38 + t21 * t37) + t118; m(6) * (t20 * t40 + t21 * t39 + t22 * t42 + t23 * t41) + t118; m(6) * ((t113 * t41 - t20) * t108 + (t113 * t42 + t21) * t107); (t42 * t20 + t41 * t21 + (-t107 * t54 + t108 * t53) * t5 * t114) * t155 - (-t115 * t36 + (t107 * t12 + t108 * t13) * t114) * t132 + (-t8 * t115 + ((t76 * t29 + t77 * t31 + t57 * t50 + t58 * t52) * t108 - t13 * t141 + (t76 * t30 + t77 * t32 + t57 * t49 + t58 * t51) * t107 + t12 * t140 + ((t108 * t27 - t48 * t141) * t108 + (t108 * t28 - t47 * t141) * t107) * t114) * t114) * t143 + (-t115 * t35 + (t10 * t107 + t108 * t11) * t114) * t131 + (-t9 * t115 + ((-t121 * t29 + t75 * t31 + t59 * t50 + t60 * t52) * t108 - t11 * t141 + (-t121 * t30 + t75 * t32 + t59 * t49 + t60 * t51) * t107 + t10 * t140 + ((t107 * t27 + t48 * t140) * t108 + (t107 * t28 + t47 * t140) * t107) * t114) * t114) * t144 - t115 * (-t147 + ((t113 * t16 + t4) * t108 + (-t113 * t17 + t3) * t107) * t114);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
