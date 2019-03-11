% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPPRR5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:37:42
% EndTime: 2019-03-09 01:37:46
% DurationCPUTime: 1.64s
% Computational Cost: add. (2663->255), mult. (4516->385), div. (0->0), fcn. (2454->6), ass. (0->142)
t70 = sin(pkin(9));
t71 = cos(pkin(9));
t54 = t71 * qJD(2) + t70 * qJD(3);
t75 = sin(qJ(5));
t77 = cos(qJ(5));
t97 = pkin(5) * t75 - pkin(8) * t77;
t32 = t97 * qJD(5) - t54;
t26 = t32 * qJD(1);
t64 = qJD(1) * qJ(2) + qJD(3);
t57 = pkin(3) * qJD(1) + t64;
t73 = -pkin(1) - qJ(3);
t58 = t73 * qJD(1) + qJD(2);
t28 = t70 * t57 + t71 * t58;
t25 = qJD(1) * pkin(7) + t28;
t16 = qJD(4) * t77 - t25 * t75;
t55 = qJD(2) * t70 - qJD(3) * t71;
t46 = t55 * qJD(1);
t7 = t16 * qJD(5) + t77 * t46;
t74 = sin(qJ(6));
t76 = cos(qJ(6));
t17 = qJD(4) * t75 + t25 * t77;
t14 = qJD(5) * pkin(8) + t17;
t27 = t71 * t57 - t70 * t58;
t89 = -pkin(5) * t77 - pkin(8) * t75 - pkin(4);
t15 = t89 * qJD(1) - t27;
t92 = t14 * t74 - t15 * t76;
t1 = -t92 * qJD(6) + t74 * t26 + t76 * t7;
t132 = qJD(1) * t77;
t59 = -qJD(6) + t132;
t171 = -t92 * t59 + t1;
t127 = qJD(6) * t74;
t115 = t75 * t127;
t125 = t76 * qJD(5);
t116 = t77 * t125;
t170 = t115 - t116;
t45 = qJD(1) * t54;
t169 = 0.2e1 * t45;
t6 = t14 * t76 + t15 * t74;
t2 = -qJD(6) * t6 + t76 * t26 - t74 * t7;
t168 = t6 * t59 - t2;
t91 = t16 * t75 - t17 * t77;
t167 = t91 * qJD(1) - t45;
t68 = t75 ^ 2;
t90 = qJD(1) * t68 - t59 * t77;
t166 = -t59 * t115 - t90 * t125;
t123 = qJD(5) * qJD(6);
t126 = qJD(6) * t76;
t114 = t75 * t126;
t129 = qJD(5) * t77;
t84 = t74 * t129 + t114;
t30 = qJD(1) * t84 + t74 * t123;
t165 = 0.2e1 * qJD(1);
t8 = t17 * qJD(5) + t75 * t46;
t162 = t74 * t8;
t161 = t76 * t8;
t160 = t8 * t75;
t13 = -qJD(5) * pkin(5) - t16;
t159 = t13 * t74;
t158 = t13 * t76;
t133 = qJD(1) * t75;
t47 = t74 * t133 - t125;
t157 = t47 * t59;
t156 = t47 * t75;
t131 = qJD(5) * t74;
t49 = t76 * t133 + t131;
t155 = t49 * t47;
t154 = t49 * t59;
t153 = t59 * t74;
t152 = t59 * t76;
t79 = qJD(1) ^ 2;
t151 = t70 * t79;
t150 = t71 * t79;
t149 = t74 * t77;
t148 = t75 * t76;
t147 = t76 * t77;
t146 = t77 * t30;
t78 = qJD(5) ^ 2;
t145 = t78 * t75;
t144 = t78 * t77;
t143 = -t47 * t116 - t30 * t148;
t117 = t75 * t125;
t39 = -t70 * t149 - t71 * t76;
t42 = t71 * t147 + t70 * t74;
t142 = t42 * qJD(1) + qJD(6) * t39 - t70 * t117;
t41 = t70 * t147 - t71 * t74;
t88 = t71 * t149 - t70 * t76;
t141 = -t41 * qJD(1) - t88 * qJD(6) - t71 * t117;
t130 = qJD(5) * t75;
t118 = t74 * t130;
t140 = -t88 * qJD(1) - qJD(6) * t41 + t70 * t118;
t139 = -t39 * qJD(1) - t42 * qJD(6) + t71 * t118;
t72 = pkin(3) + qJ(2);
t138 = t70 * t72 + t71 * t73;
t69 = t77 ^ 2;
t137 = t68 - t69;
t136 = t68 + t69;
t135 = t78 + t79;
t24 = -qJD(1) * pkin(4) - t27;
t134 = qJD(1) * t24;
t128 = qJD(6) * t47;
t124 = qJD(1) * qJD(5);
t122 = t75 * t79 * t77;
t121 = t49 * t129;
t120 = t70 * t129;
t119 = t71 * t129;
t113 = t59 * t126;
t66 = qJD(2) * t165;
t112 = 0.2e1 * t124;
t111 = t135 * t77;
t110 = t135 * t75;
t108 = t75 * t124;
t29 = t170 * qJD(1) - t76 * t123;
t107 = t49 * t130 + t29 * t77;
t106 = -t70 * t73 + t71 * t72;
t105 = -t46 - t134;
t104 = -t29 + t128;
t102 = t75 * t113;
t101 = t49 * t114;
t100 = t70 * t112;
t99 = -0.2e1 * t71 * t124;
t98 = t77 * t108;
t44 = pkin(7) + t138;
t96 = -qJD(6) * t44 * t77 + t32;
t95 = -t6 * t74 + t76 * t92;
t94 = -t6 * t76 - t74 * t92;
t87 = -t44 * t78 + t169;
t43 = -pkin(4) - t106;
t86 = qJD(5) * (qJD(1) * t43 + t24 - t55);
t85 = t90 * t74;
t31 = -t106 + t89;
t83 = -qJD(6) * t31 + t44 * t130 - t55 * t77;
t82 = t95 * qJD(6) + t1 * t76 - t2 * t74;
t81 = t7 * t77 + t160 + (-t16 * t77 - t17 * t75) * qJD(5);
t80 = t81 + t134;
t53 = t97 * qJD(1);
t12 = t44 * t147 + t31 * t74;
t11 = -t44 * t149 + t31 * t76;
t10 = t16 * t76 + t53 * t74;
t9 = -t16 * t74 + t53 * t76;
t4 = t83 * t74 + t96 * t76;
t3 = t96 * t74 - t83 * t76;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, qJ(2) * t66, 0, 0, 0, 0, 0, 0, t66, 0, qJD(3) * t165, t64 * qJD(2) - t58 * qJD(3) + (qJ(2) * qJD(2) - qJD(3) * t73) * qJD(1), 0, 0, 0, 0, 0, 0, t169, -0.2e1 * t46, 0, t45 * t106 + t46 * t138 + t27 * t54 + t28 * t55, 0.2e1 * t98, -t137 * t112, t144, -0.2e1 * t98, -t145, 0, t75 * t86 + t77 * t87, -t75 * t87 + t77 * t86, t136 * t46 + t81, -t24 * t54 - t43 * t45 + t81 * t44 - t91 * t55, -t29 * t148 - t170 * t49, -t101 + (-t121 + (t29 + t128) * t75) * t74 + t143, t107 - t166, t30 * t74 * t75 + t47 * t84, t102 + t146 + (-t85 - t156) * qJD(5) (-t59 - t132) * t130, -t4 * t59 + (-t2 + (t44 * t47 + t159) * qJD(5)) * t77 + (t13 * t126 + t30 * t44 + t47 * t55 + t162 + (qJD(1) * t11 - t92) * qJD(5)) * t75, t3 * t59 + (t1 + (t44 * t49 + t158) * qJD(5)) * t77 + (-t13 * t127 - t29 * t44 + t49 * t55 + t161 + (-qJD(1) * t12 - t6) * qJD(5)) * t75, t11 * t29 - t12 * t30 - t3 * t47 - t4 * t49 + t95 * t129 + (qJD(6) * t94 - t1 * t74 - t2 * t76) * t75, t44 * t160 + t1 * t12 + t11 * t2 + t3 * t6 - t4 * t92 + (t44 * t129 + t55 * t75) * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, -t79 * qJ(2), 0, 0, 0, 0, 0, 0, -t79, 0, 0 (-qJD(3) - t64) * qJD(1), 0, 0, 0, 0, 0, 0, -t150, t151, 0, -t45 * t70 + t46 * t71 + (-t27 * t71 - t28 * t70) * qJD(1), 0, 0, 0, 0, 0, 0, t75 * t100 - t71 * t111, t77 * t100 + t71 * t110, -t136 * t151, t167 * t70 + t80 * t71, 0, 0, 0, 0, 0, 0, t47 * t119 - t139 * t59 + (t30 * t71 + (-qJD(5) * t88 - t47 * t70) * qJD(1)) * t75, t49 * t119 + t141 * t59 + (-t29 * t71 + (-qJD(5) * t42 - t49 * t70) * qJD(1)) * t75, -t139 * t49 - t141 * t47 - t29 * t88 - t42 * t30, t71 * t160 + t1 * t42 - t2 * t88 + t141 * t6 - t139 * t92 + (-t70 * t133 + t119) * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79 (qJD(2) + t58) * qJD(1), 0, 0, 0, 0, 0, 0, -t151, -t150, 0, t45 * t71 + t46 * t70 + (-t27 * t70 + t28 * t71) * qJD(1), 0, 0, 0, 0, 0, 0, -t70 * t111 + t75 * t99, t70 * t110 + t77 * t99, t136 * t150, -t167 * t71 + t80 * t70, 0, 0, 0, 0, 0, 0, t47 * t120 - t140 * t59 + (t30 * t70 + (qJD(5) * t39 + t47 * t71) * qJD(1)) * t75, t49 * t120 + t142 * t59 + (-t29 * t70 + (-qJD(5) * t41 + t49 * t71) * qJD(1)) * t75, -t140 * t49 - t142 * t47 + t39 * t29 - t41 * t30, t70 * t160 + t1 * t41 + t2 * t39 + t142 * t6 - t140 * t92 + (t71 * t133 + t120) * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t145, -t144, 0, -t91 * qJD(5) + t7 * t75 - t8 * t77, 0, 0, 0, 0, 0, 0, t102 - t146 + (-t85 + t156) * qJD(5), t107 + t166, t101 + (t104 * t75 + t121) * t74 + t143 (-qJD(5) * t94 - t8) * t77 + (qJD(5) * t13 + t82) * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122, t137 * t79, 0, t122, 0, 0, t105 * t75, t105 * t77, 0, 0, -t49 * t152 - t29 * t74 (-t29 + t157) * t76 + (-t30 + t154) * t74, -t113 + (t59 * t147 + (-t49 + t131) * t75) * qJD(1), -t47 * t153 - t30 * t76, t59 * t127 + (-t59 * t149 + (t47 + t125) * t75) * qJD(1), t59 * t133, -pkin(5) * t30 - t17 * t47 + t59 * t9 - t161 + (pkin(8) * t152 + t159) * qJD(6) + (t92 * t75 + (-pkin(8) * t130 - t13 * t77) * t74) * qJD(1), pkin(5) * t29 - t10 * t59 - t17 * t49 + t162 + (-pkin(8) * t153 + t158) * qJD(6) + (-t13 * t147 + (-pkin(8) * t125 + t6) * t75) * qJD(1), t10 * t47 + t49 * t9 + ((qJD(6) * t49 - t30) * pkin(8) + t171) * t76 + (pkin(8) * t104 + t168) * t74, -pkin(5) * t8 + pkin(8) * t82 - t10 * t6 - t13 * t17 + t9 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, -t47 ^ 2 + t49 ^ 2, -t29 - t157, -t155, -t154 - t30, t108, -t13 * t49 - t168, t13 * t47 - t171, 0, 0;];
tauc_reg  = t5;
