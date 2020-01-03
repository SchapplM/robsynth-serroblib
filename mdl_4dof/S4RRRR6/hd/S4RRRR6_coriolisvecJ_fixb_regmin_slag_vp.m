% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% 
% Output:
% tauc_reg [4x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRRR6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:30:51
% EndTime: 2019-12-31 17:30:57
% DurationCPUTime: 1.91s
% Computational Cost: add. (1880->265), mult. (5306->412), div. (0->0), fcn. (4025->8), ass. (0->133)
t94 = cos(pkin(4));
t152 = qJD(1) * t94;
t124 = qJD(2) + t152;
t93 = sin(pkin(4));
t153 = qJD(1) * t93;
t97 = sin(qJ(2));
t138 = t97 * t153;
t96 = sin(qJ(3));
t99 = cos(qJ(3));
t177 = t99 * t124 - t96 * t138;
t48 = qJD(4) - t177;
t100 = cos(qJ(2));
t143 = qJD(1) * t100;
t135 = t93 * t143;
t84 = -qJD(3) + t135;
t156 = t100 * t93;
t141 = pkin(6) * t156;
t172 = pkin(1) * t97;
t62 = t141 + (pkin(7) + t172) * t94;
t63 = (-pkin(2) * t100 - pkin(7) * t97 - pkin(1)) * t93;
t176 = t99 * t62 + t96 * t63;
t175 = t100 * t97;
t53 = t96 * t124 + t99 * t138;
t140 = pkin(1) * t152;
t67 = pkin(6) * t135 + t97 * t140;
t40 = t124 * pkin(7) + t67;
t47 = qJD(1) * t63;
t20 = t40 * t99 + t47 * t96;
t107 = (pkin(2) * t97 - pkin(7) * t100) * t93;
t66 = qJD(2) * t107;
t58 = qJD(1) * t66;
t171 = pkin(1) * t100;
t163 = t93 * t97;
t88 = pkin(6) * t163;
t68 = (t94 * t171 - t88) * qJD(2);
t59 = qJD(1) * t68;
t102 = -t20 * qJD(3) + t99 * t58 - t96 * t59;
t142 = qJD(1) * qJD(2);
t132 = t93 * t142;
t122 = t97 * t132;
t6 = -pkin(3) * t122 - t102;
t174 = (pkin(3) * t53 + t48 * pkin(8)) * t48 + t6;
t121 = t100 * t132;
t33 = qJD(3) * t53 + t96 * t121;
t95 = sin(qJ(4));
t98 = cos(qJ(4));
t31 = t53 * t98 - t84 * t95;
t32 = t177 * qJD(3) + t99 * t121;
t12 = t31 * qJD(4) - t98 * t122 + t95 * t32;
t173 = -qJD(3) * t176 + t99 * t66 - t96 * t68;
t101 = qJD(1) ^ 2;
t144 = qJD(4) * t98;
t145 = qJD(4) * t95;
t11 = t95 * t122 - t84 * t144 - t53 * t145 + t98 * t32;
t170 = t11 * t95;
t29 = t53 * t95 + t98 * t84;
t169 = t29 * t48;
t168 = t31 * t48;
t167 = t177 * t84;
t166 = t53 * t84;
t165 = t84 * t96;
t164 = t84 * t99;
t162 = t95 * t33;
t161 = t98 * t33;
t160 = t67 + t84 * (pkin(3) * t96 - pkin(8) * t99);
t64 = -pkin(6) * t138 + t100 * t140;
t65 = qJD(1) * t107;
t158 = t99 * t64 + t96 * t65;
t157 = -t100 ^ 2 + t97 ^ 2;
t155 = t100 * t99;
t90 = t93 ^ 2;
t154 = t101 * t90;
t151 = qJD(2) * t97;
t150 = qJD(2) * t99;
t149 = qJD(3) * t95;
t148 = qJD(3) * t96;
t147 = qJD(3) * t98;
t146 = qJD(3) * t99;
t139 = t95 * t156;
t137 = t93 * t151;
t136 = t93 * t94 * t101;
t134 = qJD(2) * t156;
t133 = t90 * t142;
t128 = t48 * t98;
t83 = -pkin(3) * t99 - pkin(8) * t96 - pkin(2);
t127 = pkin(8) * t138 - qJD(4) * t83 + t158;
t125 = 0.2e1 * t133;
t123 = qJD(2) + 0.2e1 * t152;
t43 = (t98 * t155 + t95 * t97) * t153;
t119 = t98 * t146 - t43;
t118 = -0.2e1 * pkin(1) * t133;
t69 = (t94 * t172 + t141) * qJD(2);
t60 = qJD(1) * t69;
t10 = t33 * pkin(3) - t32 * pkin(8) + t60;
t106 = t47 * t146 - t40 * t148 + t96 * t58 + t99 * t59;
t5 = pkin(8) * t122 + t106;
t117 = t95 * t10 + t98 * t5;
t39 = -t124 * pkin(2) - t64;
t16 = -pkin(3) * t177 - t53 * pkin(8) + t39;
t18 = -pkin(8) * t84 + t20;
t3 = t16 * t98 - t18 * t95;
t4 = t16 * t95 + t18 * t98;
t61 = t88 + (-pkin(2) - t171) * t94;
t70 = t96 * t163 - t94 * t99;
t71 = t99 * t163 + t94 * t96;
t21 = pkin(3) * t70 - pkin(8) * t71 + t61;
t23 = -pkin(8) * t156 + t176;
t116 = t21 * t98 - t23 * t95;
t115 = t21 * t95 + t23 * t98;
t19 = -t40 * t96 + t47 * t99;
t113 = -t62 * t96 + t63 * t99;
t112 = -t64 * t96 + t65 * t99;
t36 = t98 * t156 + t71 * t95;
t109 = -t48 * t144 - t162;
t108 = -t48 * t145 + t161;
t105 = t63 * t146 - t62 * t148 + t96 * t66 + t99 * t68;
t104 = pkin(1) * (-t94 * t142 + t154);
t17 = pkin(3) * t84 - t19;
t103 = -pkin(8) * t33 + (t17 + t19) * t48;
t2 = -t4 * qJD(4) + t98 * t10 - t95 * t5;
t42 = t95 * t99 * t135 - t98 * t138;
t37 = t71 * t98 - t139;
t35 = -t70 * qJD(3) + t99 * t134;
t34 = t71 * qJD(3) + t96 * t134;
t24 = -pkin(3) * t138 - t112;
t22 = pkin(3) * t156 - t113;
t15 = -t36 * qJD(4) + t95 * t137 + t35 * t98;
t14 = -qJD(4) * t139 - t98 * t137 + t71 * t144 + t35 * t95;
t13 = t34 * pkin(3) - t35 * pkin(8) + t69;
t8 = -pkin(3) * t137 - t173;
t7 = pkin(8) * t137 + t105;
t1 = t3 * qJD(4) + t117;
t9 = [0, 0, 0, t125 * t175, -t157 * t125, t123 * t134, -t123 * t137, 0, t97 * t118 - t69 * t124 - t60 * t94, t100 * t118 - t68 * t124 - t59 * t94, t32 * t71 + t35 * t53, t177 * t35 - t32 * t70 - t33 * t71 - t34 * t53, -t35 * t84 + (-t100 * t32 + (qJD(1) * t71 + t53) * t151) * t93, t34 * t84 + (t100 * t33 + (-qJD(1) * t70 + t177) * t151) * t93, (-t90 * t143 - t84 * t93) * t151, -t173 * t84 - t69 * t177 + t61 * t33 + t60 * t70 + t39 * t34 + (-t102 * t100 + (qJD(1) * t113 + t19) * t151) * t93, t105 * t84 + t69 * t53 + t61 * t32 + t60 * t71 + t39 * t35 + (t106 * t100 + (-qJD(1) * t176 - t20) * t151) * t93, t11 * t37 + t15 * t31, -t11 * t36 - t12 * t37 - t14 * t31 - t15 * t29, t11 * t70 + t15 * t48 + t31 * t34 + t33 * t37, -t12 * t70 - t14 * t48 - t29 * t34 - t33 * t36, t33 * t70 + t34 * t48, (-qJD(4) * t115 + t98 * t13 - t95 * t7) * t48 + t116 * t33 + t2 * t70 + t3 * t34 + t8 * t29 + t22 * t12 + t6 * t36 + t17 * t14, -(qJD(4) * t116 + t95 * t13 + t98 * t7) * t48 - t115 * t33 - t1 * t70 - t4 * t34 + t8 * t31 + t22 * t11 + t6 * t37 + t17 * t15; 0, 0, 0, -t154 * t175, t157 * t154, -t100 * t136, t97 * t136, 0, -pkin(6) * t121 + t97 * t104 + t67 * t124, pkin(6) * t122 + t100 * t104 + t64 * t124, -t53 * t164 + t32 * t96, (t32 - t167) * t99 + (-t33 + t166) * t96, -t84 * t146 + (t84 * t155 + (t96 * qJD(2) - t53) * t97) * t153, t84 * t148 + (-t100 * t165 + (-t177 + t150) * t97) * t153, t84 * t138, -pkin(2) * t33 - t60 * t99 + t112 * t84 + t67 * t177 + (pkin(7) * t164 + t39 * t96) * qJD(3) + (-t19 * t97 + (-pkin(7) * t151 - t100 * t39) * t96) * t153, -pkin(2) * t32 + t60 * t96 - t158 * t84 - t67 * t53 + (-pkin(7) * t165 + t39 * t99) * qJD(3) + (-t39 * t155 + (-pkin(7) * t150 + t20) * t97) * t153, t11 * t98 * t96 + (-t96 * t145 + t119) * t31, t29 * t43 + t31 * t42 + (-t29 * t98 - t31 * t95) * t146 + (-t170 - t12 * t98 + (t29 * t95 - t31 * t98) * qJD(4)) * t96, -t11 * t99 + t119 * t48 + (-t31 * t84 + t108) * t96, t12 * t99 + (-t95 * t146 + t42) * t48 + (t29 * t84 + t109) * t96, -t48 * t165 - t33 * t99, t83 * t161 - t17 * t42 - t24 * t29 + (t127 * t95 - t160 * t98) * t48 + (t17 * t149 - t2 + (qJD(3) * t29 + t109) * pkin(7)) * t99 + (t17 * t144 + t6 * t95 - t84 * t3 + (t48 * t149 + t12) * pkin(7)) * t96, -t83 * t162 - t17 * t43 - t24 * t31 + (t127 * t98 + t160 * t95) * t48 + (t17 * t147 + t1 + (qJD(3) * t31 - t108) * pkin(7)) * t99 + (-t17 * t145 + t6 * t98 + t84 * t4 + (t48 * t147 + t11) * pkin(7)) * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53 * t177, -t177 ^ 2 + t53 ^ 2, t32 + t167, -t166 - t33, t122, -t20 * t84 - t39 * t53 + t102, -t177 * t39 - t19 * t84 - t106, t128 * t31 + t170, (t11 - t169) * t98 + (-t12 - t168) * t95, t128 * t48 - t31 * t53 + t162, -t48 ^ 2 * t95 + t29 * t53 + t161, -t48 * t53, -pkin(3) * t12 + t103 * t95 - t174 * t98 - t20 * t29 - t3 * t53, -pkin(3) * t11 + t103 * t98 + t174 * t95 - t20 * t31 + t4 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31 * t29, -t29 ^ 2 + t31 ^ 2, t11 + t169, -t12 + t168, t33, -t17 * t31 + t4 * t48 + t2, t17 * t29 - t117 + (-qJD(4) + t48) * t3;];
tauc_reg = t9;
