% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauc_reg [5x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRR10_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR10_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:10:44
% EndTime: 2019-12-31 19:10:52
% DurationCPUTime: 2.22s
% Computational Cost: add. (2700->276), mult. (7173->386), div. (0->0), fcn. (5616->8), ass. (0->146)
t201 = qJD(4) + qJD(5);
t127 = cos(pkin(9));
t133 = cos(qJ(3));
t169 = t133 * t127;
t117 = qJD(1) * t169;
t126 = sin(pkin(9));
t130 = sin(qJ(3));
t174 = t126 * t130;
t157 = qJD(1) * t174;
t96 = t117 - t157;
t211 = t201 - t96;
t191 = pkin(6) + qJ(2);
t110 = t191 * t126;
t103 = qJD(1) * t110;
t111 = t191 * t127;
t104 = qJD(1) * t111;
t69 = -t130 * t103 + t133 * t104;
t210 = qJD(3) * t69;
t129 = sin(qJ(4));
t131 = cos(qJ(5));
t128 = sin(qJ(5));
t132 = cos(qJ(4));
t170 = t128 * t132;
t106 = t129 * t131 + t170;
t187 = t211 * t106;
t163 = t132 * qJD(3);
t102 = t126 * t133 + t127 * t130;
t97 = t102 * qJD(1);
t80 = t129 * t97 - t163;
t82 = qJD(3) * t129 + t132 * t97;
t146 = t128 * t80 - t131 * t82;
t27 = t128 * t82 + t131 * t80;
t209 = t146 * t27;
t208 = t146 ^ 2 - t27 ^ 2;
t164 = qJD(5) * t131;
t165 = qJD(5) * t128;
t167 = qJD(4) * t129;
t116 = qJD(3) * t117;
t88 = -qJD(3) * t157 + t116;
t38 = qJD(4) * t163 + t132 * t88 - t97 * t167;
t39 = t82 * qJD(4) + t129 * t88;
t8 = -t128 * t39 + t131 * t38 - t80 * t164 - t82 * t165;
t92 = qJD(4) - t96;
t87 = qJD(5) + t92;
t207 = t27 * t87 + t8;
t119 = -pkin(2) * t127 - pkin(1);
t109 = t119 * qJD(1) + qJD(2);
t43 = -pkin(3) * t96 - pkin(7) * t97 + t109;
t64 = qJD(3) * pkin(7) + t69;
t21 = t129 * t43 + t132 * t64;
t15 = -pkin(8) * t80 + t21;
t13 = t15 * t165;
t202 = -t103 * t133 - t130 * t104;
t63 = -qJD(3) * pkin(3) - t202;
t24 = pkin(4) * t80 + t63;
t206 = t24 * t27 + t13;
t135 = t146 * qJD(5) - t128 * t38 - t131 * t39;
t205 = -t146 * t87 + t135;
t101 = -t169 + t174;
t138 = t101 * qJD(2);
t31 = -qJD(1) * t138 + qJD(3) * t202;
t99 = t102 * qJD(3);
t89 = qJD(1) * t99;
t51 = pkin(3) * t89 - pkin(7) * t88;
t48 = t132 * t51;
t136 = -t21 * qJD(4) - t129 * t31 + t48;
t2 = pkin(4) * t89 - pkin(8) * t38 + t136;
t166 = qJD(4) * t132;
t140 = t129 * t51 + t132 * t31 + t43 * t166 - t64 * t167;
t5 = -pkin(8) * t39 + t140;
t159 = -t128 * t5 + t131 * t2;
t20 = -t129 * t64 + t132 * t43;
t14 = -pkin(8) * t82 + t20;
t10 = pkin(4) * t92 + t14;
t180 = t131 * t15;
t4 = t10 * t128 + t180;
t204 = -t4 * qJD(5) + t24 * t146 + t159;
t203 = t132 * t89 - t92 * t167;
t78 = t110 * t133 + t130 * t111;
t171 = t128 * t129;
t105 = -t131 * t132 + t171;
t188 = t211 * t105;
t200 = -t106 * t89 + t188 * t87;
t199 = t201 * t102;
t198 = pkin(7) + pkin(8);
t197 = t27 * t97;
t196 = t146 * t97;
t195 = t80 * t92;
t194 = t80 * t97;
t193 = t82 * t92;
t192 = t82 * t97;
t65 = pkin(3) * t97 - pkin(7) * t96;
t190 = t129 * t65 + t132 * t202;
t67 = pkin(3) * t101 - pkin(7) * t102 + t119;
t79 = -t110 * t130 + t111 * t133;
t72 = t132 * t79;
t189 = t129 * t67 + t72;
t185 = t129 * t38;
t183 = t129 * t89;
t182 = t129 * t96;
t98 = t101 * qJD(3);
t181 = t129 * t98;
t179 = t132 * t98;
t178 = t102 * t129;
t177 = t102 * t132;
t168 = t126 ^ 2 + t127 ^ 2;
t162 = qJD(1) * qJD(2);
t158 = qJD(4) * t198;
t156 = t102 * t166;
t155 = qJD(5) * t10 + t5;
t153 = t132 * t92;
t152 = t168 * qJD(1) ^ 2;
t32 = t102 * t162 + t210;
t151 = -t105 * t89 - t187 * t87;
t150 = pkin(7) * qJD(4) * t92 + t32;
t149 = -t69 + (t167 - t182) * pkin(4);
t112 = t198 * t129;
t148 = -pkin(8) * t182 + qJD(5) * t112 + t129 * t158 + t190;
t113 = t198 * t132;
t55 = t132 * t65;
t147 = pkin(4) * t97 + qJD(5) * t113 - t129 * t202 + t55 + (-pkin(8) * t96 + t158) * t132;
t144 = t92 * t182 + t203;
t143 = 0.2e1 * t168 * t162;
t142 = -pkin(7) * t89 + t92 * t63;
t141 = t156 - t181;
t44 = -t78 * qJD(3) - t138;
t66 = pkin(3) * t99 + pkin(7) * t98;
t139 = t129 * t66 + t132 * t44 + t67 * t166 - t79 * t167;
t45 = t102 * qJD(2) + t79 * qJD(3);
t122 = -pkin(4) * t132 - pkin(3);
t70 = t89 * t101;
t60 = t132 * t67;
t58 = t105 * t102;
t57 = t106 * t102;
t56 = t132 * t66;
t46 = pkin(4) * t178 + t78;
t23 = pkin(4) * t141 + t45;
t22 = -pkin(8) * t178 + t189;
t18 = pkin(4) * t101 - pkin(8) * t177 - t129 * t79 + t60;
t17 = pkin(4) * t39 + t32;
t12 = -t98 * t170 + (t201 * t177 - t181) * t131 - t171 * t199;
t11 = t105 * t98 - t106 * t199;
t7 = -pkin(8) * t141 + t139;
t6 = pkin(8) * t179 + pkin(4) * t99 - t129 * t44 + t56 + (-t72 + (pkin(8) * t102 - t67) * t129) * qJD(4);
t3 = t131 * t10 - t128 * t15;
t1 = [0, 0, 0, 0, 0, t143, qJ(2) * t143, t102 * t88 - t97 * t98, -t101 * t88 - t102 * t89 - t96 * t98 - t97 * t99, -t98 * qJD(3), -t99 * qJD(3), 0, -qJD(3) * t45 + t109 * t99 + t119 * t89, -qJD(3) * t44 - t109 * t98 + t119 * t88, -t82 * t179 + (t132 * t38 - t82 * t167) * t102, -(-t129 * t82 - t132 * t80) * t98 + (-t185 - t132 * t39 + (t129 * t80 - t132 * t82) * qJD(4)) * t102, t101 * t38 + t203 * t102 - t92 * t179 + t82 * t99, t92 * t181 - t101 * t39 - t80 * t99 + (-t92 * t166 - t183) * t102, t92 * t99 + t70, (-t79 * t166 + t56) * t92 + t60 * t89 + (-t64 * t166 + t48) * t101 + t20 * t99 + t45 * t80 + t78 * t39 + t63 * t156 + ((-qJD(4) * t67 - t44) * t92 - t79 * t89 + (-qJD(4) * t43 - t31) * t101 + t32 * t102 - t63 * t98) * t129, -t139 * t92 - t189 * t89 - t140 * t101 - t21 * t99 + t45 * t82 + t78 * t38 - t63 * t179 + (t32 * t132 - t63 * t167) * t102, -t11 * t146 - t58 * t8, -t11 * t27 + t12 * t146 - t135 * t58 - t57 * t8, t101 * t8 + t11 * t87 - t146 * t99 - t58 * t89, t101 * t135 - t12 * t87 - t27 * t99 - t57 * t89, t87 * t99 + t70, (-t128 * t7 + t131 * t6) * t87 + (-t128 * t22 + t131 * t18) * t89 + t159 * t101 + t3 * t99 + t23 * t27 - t46 * t135 + t17 * t57 + t24 * t12 + ((-t128 * t18 - t131 * t22) * t87 - t4 * t101) * qJD(5), t13 * t101 + t24 * t11 - t17 * t58 - t23 * t146 - t4 * t99 + t46 * t8 + (-(-qJD(5) * t22 + t6) * t87 - t18 * t89 - t2 * t101) * t128 + (-(qJD(5) * t18 + t7) * t87 - t22 * t89 - t155 * t101) * t131; 0, 0, 0, 0, 0, -t152, -qJ(2) * t152, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t97, t116 + (t96 - t157) * qJD(3), 0, 0, 0, 0, 0, t144 - t194, -t92 ^ 2 * t132 - t183 - t192, 0, 0, 0, 0, 0, t151 - t197, t196 + t200; 0, 0, 0, 0, 0, 0, 0, -t97 * t96, -t96 ^ 2 + t97 ^ 2, t116 + (-t96 - t157) * qJD(3), 0, 0, -t109 * t97 + t210 - t32, t101 * t162 - t109 * t96, t82 * t153 + t185, (t38 - t195) * t132 + (-t39 - t193) * t129, t92 * t153 + t183 - t192, t144 + t194, -t92 * t97, -pkin(3) * t39 - t20 * t97 - t55 * t92 - t69 * t80 - t150 * t132 + (t202 * t92 + t142) * t129, -pkin(3) * t38 + t150 * t129 + t142 * t132 + t190 * t92 + t21 * t97 - t69 * t82, t106 * t8 + t146 * t188, -t105 * t8 + t106 * t135 + t146 * t187 + t188 * t27, t196 - t200, t151 + t197, -t87 * t97, (-t112 * t131 - t113 * t128) * t89 - t122 * t135 + t17 * t105 - t3 * t97 + (t128 * t148 - t131 * t147) * t87 + t149 * t27 + t187 * t24, -(-t112 * t128 + t113 * t131) * t89 + t122 * t8 + t17 * t106 + t4 * t97 + (t128 * t147 + t131 * t148) * t87 - t149 * t146 - t188 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82 * t80, -t80 ^ 2 + t82 ^ 2, t38 + t195, t193 - t39, t89, t21 * t92 - t63 * t82 + t136, t20 * t92 + t63 * t80 - t140, -t209, t208, t207, t205, t89, -(-t128 * t14 - t180) * t87 + (t131 * t89 - t165 * t87 - t27 * t82) * pkin(4) + t204, (-t15 * t87 - t2) * t128 + (t14 * t87 - t155) * t131 + (-t128 * t89 + t146 * t82 - t164 * t87) * pkin(4) + t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t209, t208, t207, t205, t89, t4 * t87 + t204, -t128 * t2 - t131 * t155 + t3 * t87 + t206;];
tauc_reg = t1;
