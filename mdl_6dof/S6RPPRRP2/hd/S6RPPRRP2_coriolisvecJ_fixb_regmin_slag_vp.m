% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRRP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:01:17
% EndTime: 2019-03-09 02:01:22
% DurationCPUTime: 1.70s
% Computational Cost: add. (3519->273), mult. (8382->354), div. (0->0), fcn. (6065->8), ass. (0->138)
t113 = sin(qJ(4));
t115 = cos(qJ(4));
t110 = cos(pkin(10));
t104 = t110 * qJD(2);
t108 = sin(pkin(10));
t100 = sin(pkin(9)) * pkin(1) + qJ(3);
t96 = t100 * qJD(1);
t67 = t104 + (-pkin(7) * qJD(1) - t96) * t108;
t158 = qJD(1) * t110;
t78 = t108 * qJD(2) + t110 * t96;
t68 = pkin(7) * t158 + t78;
t32 = t113 * t67 + t115 * t68;
t202 = t32 * qJD(4);
t112 = sin(qJ(5));
t155 = t112 * qJD(4);
t114 = cos(qJ(5));
t156 = qJD(5) * t114;
t93 = t115 * t108 + t113 * t110;
t194 = t93 * qJD(1);
t163 = t113 * t108;
t92 = -t115 * t110 + t163;
t195 = t92 * qJD(1);
t41 = (qJD(5) - t195) * t155 + t194 * t156;
t201 = t93 * qJD(3);
t85 = -qJD(1) * t163 + t115 * t158;
t82 = qJD(5) - t85;
t200 = t194 * qJD(4);
t29 = qJD(4) * pkin(8) + t32;
t94 = -cos(pkin(9)) * pkin(1) - t110 * pkin(3) - pkin(2);
t83 = t94 * qJD(1) + qJD(3);
t39 = -t85 * pkin(4) - pkin(8) * t194 + t83;
t10 = t112 * t39 + t114 * t29;
t8 = t82 * qJ(6) + t10;
t188 = t8 * t82;
t157 = qJD(5) * t112;
t122 = t92 * qJD(3);
t196 = -t113 * t68 + t115 * t67;
t20 = -qJD(1) * t122 + qJD(4) * t196;
t87 = t92 * qJD(4);
t120 = qJD(1) * t87;
t88 = t93 * qJD(4);
t79 = qJD(1) * t88;
t47 = t79 * pkin(4) + pkin(8) * t120;
t147 = t112 * t20 - t114 * t47 + t29 * t156 + t39 * t157;
t189 = t79 * pkin(5);
t2 = t147 - t189;
t198 = -t2 + t188;
t154 = t114 * qJD(4);
t69 = t112 * t194 - t154;
t197 = -t92 * t41 - t88 * t69;
t40 = -qJD(5) * t154 + t114 * t120 + t157 * t194;
t71 = t114 * t194 + t155;
t175 = -t92 * t40 + t71 * t88;
t180 = pkin(7) + t100;
t89 = t180 * t108;
t90 = t180 * t110;
t54 = t113 * t90 + t115 * t89;
t28 = -qJD(4) * pkin(4) - t196;
t13 = t69 * pkin(5) - t71 * qJ(6) + t28;
t190 = pkin(8) * t79;
t193 = t13 * t82 - t190;
t192 = t71 ^ 2;
t191 = t82 ^ 2;
t187 = t10 * t82;
t186 = t69 * t85;
t185 = t71 * t69;
t148 = t71 * t82;
t184 = t71 * t194;
t183 = t71 * t87;
t182 = t194 * t69;
t181 = t93 * t79;
t134 = pkin(5) * t112 - qJ(6) * t114;
t179 = t112 * qJD(6) - t82 * t134 + t32;
t170 = t114 * t93;
t171 = t114 * t87;
t178 = -t41 * t170 + t69 * t171;
t58 = pkin(4) * t194 - t85 * pkin(8);
t177 = t112 * t58 + t114 * t196;
t176 = -t112 * t41 - t69 * t156;
t51 = t92 * pkin(4) - t93 * pkin(8) + t94;
t55 = -t113 * t89 + t115 * t90;
t174 = t112 * t51 + t114 * t55;
t173 = pkin(8) * qJD(5);
t74 = t112 * t79;
t172 = t112 * t82;
t76 = t114 * t79;
t167 = t79 * qJ(6);
t9 = -t112 * t29 + t114 * t39;
t166 = qJD(6) - t9;
t165 = qJD(5) * t69;
t164 = qJD(5) * t93;
t160 = t87 * qJD(4);
t159 = t108 ^ 2 + t110 ^ 2;
t152 = t82 * t173;
t151 = t93 * t157;
t149 = t93 * t156;
t21 = qJD(1) * t201 + t202;
t146 = -t40 + t165;
t144 = qJD(1) * t159;
t143 = t82 * t151;
t142 = t71 * t149;
t4 = t41 * pkin(5) + t40 * qJ(6) - t71 * qJD(6) + t21;
t141 = -t4 - t152;
t140 = -t13 * t87 + t4 * t93;
t126 = t112 * t47 + t114 * t20 + t39 * t156 - t29 * t157;
t1 = t82 * qJD(6) + t126 + t167;
t7 = -t82 * pkin(5) + t166;
t139 = t82 * t7 + t1;
t138 = -t112 * t8 + t114 * t7;
t137 = t112 * t7 + t114 * t8;
t136 = t21 * t93 - t28 * t87;
t135 = t82 * t87 - t181;
t133 = t108 * (-t108 * t96 + t104) - t110 * t78;
t131 = t74 + (-t114 * t85 + t156) * t82;
t130 = -t82 * t157 + t85 * t172 + t76;
t129 = t13 * t71 + t147;
t128 = t82 * t28 - t190;
t127 = t79 * t170 - t82 * t171 - t143;
t35 = -qJD(4) * t54 - t122;
t59 = t88 * pkin(4) + t87 * pkin(8);
t125 = t112 * t59 + t114 * t35 + t51 * t156 - t55 * t157;
t119 = t135 * t112 - t82 * t149;
t118 = t138 * qJD(5) + t1 * t114 + t2 * t112;
t117 = t119 - t197;
t36 = t55 * qJD(4) + t201;
t95 = -t114 * pkin(5) - t112 * qJ(6) - pkin(4);
t81 = t88 * qJD(4);
t37 = t71 * pkin(5) + t69 * qJ(6);
t22 = t134 * t93 + t54;
t17 = t69 * t82 - t40;
t15 = -t92 * pkin(5) + t112 * t55 - t114 * t51;
t14 = t92 * qJ(6) + t174;
t12 = -pkin(5) * t194 + t112 * t196 - t114 * t58;
t11 = qJ(6) * t194 + t177;
t6 = (-pkin(5) * t87 + qJ(6) * t164) * t112 + (qJ(6) * t87 + (pkin(5) * qJD(5) - qJD(6)) * t93) * t114 + t36;
t5 = -t88 * pkin(5) + t174 * qJD(5) + t112 * t35 - t114 * t59;
t3 = t88 * qJ(6) + t92 * qJD(6) + t125;
t16 = [0, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t144 (t100 * t144 - t133) * qJD(3), -t93 * t120 - t194 * t87, t92 * t120 - t194 * t88 - t87 * t85 - t181, -t160, -t81, 0, -t36 * qJD(4) + t94 * t79 + t83 * t88, -t83 * t87 + (-t195 * t94 - t35) * qJD(4), -t71 * t151 + (-t40 * t93 - t183) * t114, -t142 + (t183 + (t40 + t165) * t93) * t112 + t178, t127 + t175, t119 + t197, t79 * t92 + t82 * t88, -t147 * t92 + t9 * t88 + t36 * t69 + t54 * t41 + ((-qJD(5) * t55 + t59) * t82 + t51 * t79 + t28 * t164) * t114 + ((-qJD(5) * t51 - t35) * t82 - t55 * t79 + t136) * t112, -t10 * t88 + t114 * t136 - t125 * t82 - t126 * t92 - t151 * t28 - t174 * t79 + t36 * t71 - t54 * t40, t112 * t140 + t13 * t149 - t15 * t79 - t2 * t92 + t22 * t41 - t5 * t82 + t6 * t69 - t7 * t88, -t14 * t41 - t15 * t40 - t3 * t69 + t5 * t71 - t138 * t87 + (-qJD(5) * t137 - t1 * t112 + t2 * t114) * t93, t1 * t92 - t114 * t140 + t13 * t151 + t14 * t79 + t22 * t40 + t3 * t82 - t6 * t71 + t8 * t88, t1 * t14 + t13 * t6 + t15 * t2 + t22 * t4 + t3 * t8 + t5 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, t160, 0, 0, 0, 0, 0, t117, t114 * t135 + t143 + t175, t117, t142 + (t146 * t93 - t183) * t112 + t178, t127 - t175, t118 * t93 + t13 * t88 - t137 * t87 + t4 * t92; 0, 0, 0, 0, 0, 0, -t159 * qJD(1) ^ 2, t133 * qJD(1), 0, 0, 0, 0, 0, 0.2e1 * t200 (t85 - t195) * qJD(4), 0, 0, 0, 0, 0, t130 - t182, -t114 * t191 - t184 - t74, -t172 * t82 - t182 + t76 (t40 + t186) * t114 + t112 * t148 + t176, t131 + t184, t139 * t112 + t114 * t198 - t13 * t194; 0, 0, 0, 0, 0, 0, 0, 0, -t194 * t85, t194 ^ 2 - t85 ^ 2 (-t85 - t195) * qJD(4), 0, 0, -t194 * t83 + t202 - t21, qJD(3) * t195 - t83 * t85, -t40 * t112 + t114 * t148 (-t40 + t186) * t114 - t71 * t172 + t176, t131 - t184, t130 + t182, -t82 * t194, -pkin(4) * t41 - t32 * t69 - t9 * t194 + (-t21 + (-t58 - t173) * t82) * t114 + (t196 * t82 + t128) * t112, pkin(4) * t40 + t177 * t82 + t10 * t194 - t32 * t71 + (t21 + t152) * t112 + t128 * t114, t112 * t193 + t141 * t114 + t12 * t82 - t179 * t69 + t194 * t7 + t95 * t41, t11 * t69 - t12 * t71 + ((qJD(5) * t71 - t41) * pkin(8) + t139) * t114 + (pkin(8) * t146 - t198) * t112, -t11 * t82 + t141 * t112 - t114 * t193 + t179 * t71 - t194 * t8 + t95 * t40, pkin(8) * t118 - t8 * t11 - t7 * t12 - t13 * t179 + t4 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t185, -t69 ^ 2 + t192, t17, t148 - t41, t79, -t28 * t71 - t147 + t187, t28 * t69 + t9 * t82 - t126, -t37 * t69 - t129 + t187 + 0.2e1 * t189, pkin(5) * t40 - t41 * qJ(6) + (-t10 + t8) * t71 + (t7 - t166) * t69, 0.2e1 * t167 - t13 * t69 + t37 * t71 + (0.2e1 * qJD(6) - t9) * t82 + t126, -t2 * pkin(5) + t1 * qJ(6) - t7 * t10 - t13 * t37 + t166 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t185 - t200, t17, -t191 - t192, t129 - t188 - t189;];
tauc_reg  = t16;
