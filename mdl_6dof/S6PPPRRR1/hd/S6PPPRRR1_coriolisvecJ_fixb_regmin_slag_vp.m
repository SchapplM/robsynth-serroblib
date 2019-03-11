% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
% 
% Output:
% tauc_reg [6x20]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PPPRRR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:40:34
% EndTime: 2019-03-08 18:40:38
% DurationCPUTime: 1.44s
% Computational Cost: add. (1959->208), mult. (5358->357), div. (0->0), fcn. (5240->16), ass. (0->135)
t72 = cos(pkin(13));
t74 = cos(pkin(7));
t135 = t72 * t74;
t67 = sin(pkin(13));
t71 = cos(pkin(14));
t142 = t67 * t71;
t75 = cos(pkin(6));
t60 = qJD(1) * t75 + qJD(2);
t66 = sin(pkin(14));
t69 = sin(pkin(7));
t70 = sin(pkin(6));
t36 = t69 * t66 * t60 + (t66 * t135 + t142) * t70 * qJD(1);
t78 = sin(qJ(4));
t81 = cos(qJ(4));
t138 = t69 * t71;
t87 = (t71 * t135 - t66 * t67) * t70;
t35 = qJD(1) * t87 + t60 * t138;
t73 = cos(pkin(8));
t149 = t35 * t73;
t114 = t69 * t70 * t72;
t48 = -qJD(1) * t114 + t60 * t74 + qJD(3);
t68 = sin(pkin(8));
t97 = t48 * t68 + t149;
t85 = t78 * t36 - t97 * t81;
t76 = sin(qJ(6));
t121 = qJD(6) * t76;
t77 = sin(qJ(5));
t110 = t77 * t121;
t79 = cos(qJ(6));
t119 = t79 * qJD(5);
t80 = cos(qJ(5));
t157 = t80 * t119 - t110;
t127 = qJD(4) * t78;
t113 = t68 * t127;
t125 = qJD(4) * t81;
t17 = t48 * t113 + t36 * t125 + t127 * t149;
t21 = t81 * t36 + t97 * t78;
t155 = qJD(4) * t21 - t17;
t126 = qJD(4) * t80;
t61 = -qJD(6) + t126;
t154 = -qJD(6) - t61;
t136 = t71 * t73;
t140 = t68 * t81;
t153 = (t81 * t136 - t66 * t78) * t69 + t74 * t140;
t137 = t69 * t75;
t40 = t70 * t142 + (t70 * t135 + t137) * t66;
t39 = t71 * t137 + t87;
t50 = t74 * t75 - t114;
t96 = t39 * t73 + t50 * t68;
t152 = -t40 * t78 + t96 * t81;
t117 = qJD(5) * qJD(6);
t120 = qJD(6) * t79;
t122 = qJD(5) * t80;
t47 = (t77 * t120 + t76 * t122) * qJD(4) + t76 * t117;
t19 = qJD(4) * pkin(10) + t21;
t28 = -t35 * t68 + t48 * t73;
t10 = t19 * t80 + t28 * t77;
t16 = t85 * qJD(4);
t4 = t10 * qJD(5) - t77 * t16;
t151 = t4 * t76;
t150 = t4 * t79;
t46 = t157 * qJD(4) + t79 * t117;
t147 = t46 * t76;
t128 = qJD(4) * t77;
t53 = t76 * t128 - t119;
t146 = t53 * t61;
t124 = qJD(5) * t76;
t55 = t79 * t128 + t124;
t145 = t55 * t61;
t144 = t61 * t76;
t143 = t61 * t79;
t141 = t68 * t78;
t83 = qJD(4) ^ 2;
t139 = t68 * t83;
t134 = t76 * t80;
t133 = t79 * t80;
t103 = pkin(5) * t77 - pkin(11) * t80;
t57 = t103 * qJD(5);
t132 = t21 - t57;
t64 = t77 ^ 2;
t131 = -t80 ^ 2 + t64;
t130 = qJD(4) * pkin(4);
t123 = qJD(5) * t77;
t118 = qJD(4) * qJD(5);
t115 = t78 * t139;
t112 = t68 * t125;
t109 = t61 * t120;
t107 = t77 * t118;
t18 = t85 - t130;
t106 = -qJD(4) * t18 + t16;
t105 = t80 * t112;
t104 = t77 * t112;
t58 = -pkin(5) * t80 - pkin(11) * t77 - pkin(4);
t15 = t58 * qJD(4) + t85;
t8 = qJD(5) * pkin(11) + t10;
t1 = t15 * t79 - t76 * t8;
t2 = t15 * t76 + t79 * t8;
t25 = t40 * t81 + t96 * t78;
t29 = -t39 * t68 + t50 * t73;
t12 = t25 * t80 + t29 * t77;
t102 = t12 * t79 - t152 * t76;
t101 = -t12 * t76 - t152 * t79;
t100 = t19 * t77 - t28 * t80;
t11 = t25 * t77 - t29 * t80;
t43 = t74 * t141 + (t78 * t136 + t66 * t81) * t69;
t49 = -t68 * t138 + t73 * t74;
t31 = t43 * t80 + t49 * t77;
t99 = -t153 * t76 + t31 * t79;
t98 = -t153 * t79 - t31 * t76;
t30 = t43 * t77 - t49 * t80;
t95 = qJD(4) * t64 - t61 * t80;
t52 = t80 * t141 + t73 * t77;
t94 = -t79 * t140 - t52 * t76;
t93 = t76 * t140 - t52 * t79;
t51 = t77 * t141 - t73 * t80;
t82 = qJD(5) ^ 2;
t91 = pkin(10) * t82 - t155;
t90 = qJD(5) * (t18 - t85 - t130);
t3 = -t100 * qJD(5) - t80 * t16;
t7 = -qJD(5) * pkin(5) + t100;
t84 = qJD(5) * t7 + qJD(6) * t15 + t61 * t85 + t3;
t56 = t103 * qJD(4);
t45 = t52 * qJD(5) + t104;
t44 = -t51 * qJD(5) + t105;
t38 = t43 * qJD(4);
t37 = t153 * qJD(4);
t27 = t31 * qJD(5) + t77 * t37;
t26 = -t30 * qJD(5) + t80 * t37;
t23 = t25 * qJD(4);
t22 = t152 * qJD(4);
t14 = qJD(4) * t57 + t17;
t13 = t79 * t14;
t6 = -t11 * qJD(5) + t22 * t80;
t5 = t12 * qJD(5) + t22 * t77;
t9 = [0, 0, 0, 0, -t23 * qJD(4), -t22 * qJD(4), 0, 0, 0, 0, 0, -t5 * qJD(5) + (-t123 * t152 - t23 * t80) * qJD(4), -t6 * qJD(5) + (-t122 * t152 + t23 * t77) * qJD(4), 0, 0, 0, 0, 0 -(-t102 * qJD(6) + t23 * t79 - t6 * t76) * t61 + t101 * t107 + t5 * t53 + t11 * t47 (t101 * qJD(6) + t23 * t76 + t6 * t79) * t61 - t102 * t107 + t5 * t55 + t11 * t46; 0, 0, 0, 0, -t38 * qJD(4), -t37 * qJD(4), 0, 0, 0, 0, 0, -t27 * qJD(5) + (-t123 * t153 - t38 * t80) * qJD(4), -t26 * qJD(5) + (-t122 * t153 + t38 * t77) * qJD(4), 0, 0, 0, 0, 0 -(-t99 * qJD(6) - t76 * t26 + t79 * t38) * t61 + t98 * t107 + t27 * t53 + t30 * t47 (t98 * qJD(6) + t79 * t26 + t76 * t38) * t61 - t99 * t107 + t27 * t55 + t30 * t46; 0, 0, 0, 0, -t115, -t81 * t139, 0, 0, 0, 0, 0, -t80 * t115 + (-t45 - t104) * qJD(5), t77 * t115 + (-t44 - t105) * qJD(5), 0, 0, 0, 0, 0 -(t93 * qJD(6) + t79 * t113 - t76 * t44) * t61 + t94 * t107 + t45 * t53 + t51 * t47 (t94 * qJD(6) + t76 * t113 + t79 * t44) * t61 + t93 * t107 + t45 * t55 + t51 * t46; 0, 0, 0, 0, t155, 0, 0.2e1 * t80 * t107, -0.2e1 * t131 * t118, t82 * t80, -t82 * t77, 0, t77 * t90 - t91 * t80, t91 * t77 + t80 * t90, t46 * t79 * t77 + t157 * t55 (-t53 * t79 - t55 * t76) * t122 + (-t147 - t47 * t79 + (t53 * t76 - t55 * t79) * qJD(6)) * t77, t61 * t110 - t80 * t46 + (t77 * t55 + t95 * t79) * qJD(5), t77 * t109 + t47 * t80 + (-t77 * t53 - t95 * t76) * qJD(5) (-t61 - t126) * t123 (t58 * t121 + t132 * t79) * t61 + (t8 * t120 - t13 + (qJD(5) * t53 + t109) * pkin(10) + t84 * t76) * t80 + (t7 * t120 + pkin(10) * t47 + t85 * t53 + t151 + (-pkin(10) * t144 + (-pkin(10) * t134 + t58 * t79) * qJD(4) + t1) * qJD(5)) * t77 (t58 * t120 - t132 * t76) * t61 + (qJD(5) * pkin(10) * t55 + (t14 + (-pkin(10) * t61 - t8) * qJD(6)) * t76 + t84 * t79) * t80 + (-t7 * t121 + pkin(10) * t46 + t85 * t55 + t150 + (-pkin(10) * t143 - (pkin(10) * t133 + t58 * t76) * qJD(4) - t2) * qJD(5)) * t77; 0, 0, 0, 0, 0, 0, -t77 * t83 * t80, t131 * t83, 0, 0, 0, t106 * t77, t106 * t80, -t55 * t143 + t147 (t46 + t146) * t79 + (-t47 + t145) * t76, -t109 + (t61 * t133 + (-t55 + t124) * t77) * qJD(4), t61 * t121 + (-t61 * t134 + (t53 + t119) * t77) * qJD(4), t61 * t128, -pkin(5) * t47 - t150 + (t100 * t76 + t56 * t79) * t61 - t10 * t53 + (pkin(11) * t143 + t7 * t76) * qJD(6) + (-t1 * t77 + (-pkin(11) * t123 - t7 * t80) * t76) * qJD(4), -pkin(5) * t46 + t151 - (-t100 * t79 + t56 * t76) * t61 - t10 * t55 + (-pkin(11) * t144 + t7 * t79) * qJD(6) + (-t7 * t133 + (-pkin(11) * t119 + t2) * t77) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53 * t55, -t53 ^ 2 + t55 ^ 2, t46 - t146, -t145 - t47, t107, t154 * t2 - t76 * t3 - t55 * t7 + t13, t154 * t1 - t76 * t14 - t79 * t3 + t53 * t7;];
tauc_reg  = t9;
