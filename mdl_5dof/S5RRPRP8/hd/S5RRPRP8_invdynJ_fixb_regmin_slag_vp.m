% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% tau_reg [5x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRP8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:04:29
% EndTime: 2019-12-31 20:04:33
% DurationCPUTime: 1.64s
% Computational Cost: add. (1402->277), mult. (3048->343), div. (0->0), fcn. (1900->6), ass. (0->157)
t128 = sin(qJ(2));
t186 = qJD(1) * t128;
t101 = pkin(6) * t186;
t231 = -pkin(7) * t186 + qJD(3) + t101;
t127 = sin(qJ(4));
t131 = cos(qJ(2));
t130 = cos(qJ(4));
t181 = qJD(4) * t130;
t182 = qJD(4) * t127;
t183 = qJD(2) * t131;
t230 = t127 * t183 + t128 * t181 - t131 * t182;
t133 = -pkin(2) - pkin(3);
t172 = t133 * qJD(2);
t29 = t172 + t231;
t121 = qJD(2) * qJ(3);
t185 = qJD(1) * t131;
t102 = pkin(6) * t185;
t63 = -pkin(7) * t185 + t102;
t44 = t121 + t63;
t153 = t127 * t29 + t130 * t44;
t180 = qJD(1) * qJD(2);
t169 = t131 * t180;
t105 = t128 * qJDD(1);
t96 = pkin(6) * t105;
t174 = pkin(6) * t169 + qJDD(3) + t96;
t227 = t169 + t105;
t20 = -t227 * pkin(7) + t133 * qJDD(2) + t174;
t170 = t128 * t180;
t179 = t131 * qJDD(1);
t226 = t170 - t179;
t119 = qJDD(2) * qJ(3);
t120 = qJD(2) * qJD(3);
t97 = pkin(6) * t179;
t28 = -pkin(6) * t170 + t119 + t120 + t97;
t21 = t226 * pkin(7) + t28;
t229 = -t153 * qJD(4) - t127 * t21 + t130 * t20;
t118 = qJD(2) - qJD(4);
t196 = t127 * t128;
t54 = t130 * t131 + t196;
t40 = t54 * qJD(1);
t144 = t54 * qJD(4);
t8 = qJD(1) * t144 - t226 * t127 - t227 * t130;
t224 = t118 * t40 + t8;
t214 = g(3) * t54;
t45 = -qJD(1) * pkin(1) - pkin(2) * t185 - qJ(3) * t186;
t27 = pkin(3) * t185 - t45;
t129 = sin(qJ(1));
t191 = t131 * t127;
t151 = -t128 * t130 + t191;
t34 = t151 * t129;
t132 = cos(qJ(1));
t190 = t131 * t132;
t176 = t127 * t190;
t194 = t128 * t132;
t36 = -t130 * t194 + t176;
t42 = -t127 * t185 + t130 * t186;
t228 = g(1) * t36 + g(2) * t34 - t27 * t42 + t214 + t229;
t189 = t131 * pkin(2) + t128 * qJ(3);
t225 = -pkin(1) - t189;
t9 = t230 * qJD(1) + qJDD(1) * t54 - t130 * t170;
t223 = t118 * t42 + t9;
t65 = qJ(3) * t130 + t127 * t133;
t222 = qJD(4) * t65 + t231 * t127 + t130 * t63;
t160 = -qJ(3) * t127 + t130 * t133;
t221 = t160 * qJD(4) - t127 * t63 + t231 * t130;
t217 = t42 ^ 2;
t39 = t40 ^ 2;
t220 = -t39 + t217;
t219 = g(1) * t132 + g(2) * t129;
t199 = pkin(6) * qJDD(2);
t218 = (qJD(1) * t225 + t45) * qJD(2) - t199;
t216 = pkin(6) - pkin(7);
t166 = -t127 * t44 + t130 * t29;
t204 = qJ(5) * t42;
t6 = t166 - t204;
t5 = -pkin(4) * t118 + t6;
t215 = -t6 + t5;
t115 = g(1) * t129;
t212 = g(2) * t132;
t211 = t42 * t40;
t210 = -t204 + t221;
t205 = qJ(5) * t40;
t209 = t205 - t222;
t69 = t216 * t128;
t70 = t216 * t131;
t207 = t127 * t69 + t130 * t70;
t95 = pkin(4) * t130 + pkin(3);
t201 = t131 * t95;
t106 = t128 * qJD(3);
t200 = qJ(3) * t183 + t106;
t124 = qJDD(1) * pkin(1);
t197 = qJDD(2) * pkin(2);
t195 = t128 * t129;
t135 = qJD(1) ^ 2;
t193 = t128 * t135;
t192 = t129 * t131;
t122 = t128 ^ 2;
t123 = t131 ^ 2;
t187 = t122 - t123;
t184 = qJD(2) * t128;
t178 = pkin(4) * t196;
t177 = -g(1) * t194 - g(2) * t195 + g(3) * t131;
t175 = t131 * t193;
t173 = qJD(1) * t133;
t168 = t115 - t212;
t164 = -t127 * t70 + t130 * t69;
t162 = -t96 - t177;
t161 = -qJD(2) * pkin(2) + qJD(3);
t52 = t131 * pkin(3) - t225;
t159 = t132 * pkin(1) + pkin(2) * t190 + t129 * pkin(6) + qJ(3) * t194;
t158 = t118 ^ 2;
t157 = t128 * t172;
t156 = pkin(2) * t179 + t227 * qJ(3) + qJD(1) * t106 + t124;
t134 = qJD(2) ^ 2;
t155 = pkin(6) * t134 + t212;
t66 = t101 + t161;
t68 = t102 + t121;
t152 = t128 * t68 - t131 * t66;
t32 = t174 - t197;
t149 = -0.2e1 * pkin(1) * t180 - t199;
t148 = -t127 * t20 - t130 * t21 - t29 * t181 + t182 * t44;
t62 = t216 * t184;
t64 = qJD(2) * t70;
t147 = t127 * t64 - t130 * t62 + t69 * t181 - t182 * t70;
t25 = t157 + t200;
t143 = -t155 + 0.2e1 * t124;
t141 = -t207 * qJD(4) + t127 * t62 + t130 * t64;
t19 = pkin(2) * t170 - t156;
t38 = pkin(2) * t184 - t200;
t140 = -qJD(1) * t38 - qJDD(1) * t225 - t155 - t19;
t10 = pkin(3) * t179 + qJD(1) * t157 + t156;
t139 = -qJD(2) * t152 + t32 * t128 + t28 * t131;
t137 = pkin(4) * t9 + qJDD(5) + t10;
t35 = t54 * t129;
t37 = t54 * t132;
t136 = g(1) * t37 + g(2) * t35 - g(3) * t151 + t27 * t40 + t148;
t126 = -qJ(5) - pkin(7);
t117 = qJDD(2) - qJDD(4);
t111 = t132 * pkin(6);
t92 = qJ(3) * t185;
t87 = g(1) * t192;
t83 = qJ(3) * t190;
t81 = qJ(3) * t192;
t60 = -pkin(4) + t160;
t59 = pkin(2) * t186 - t92;
t33 = t128 * t173 + t92;
t23 = qJD(2) * t54 - t144;
t22 = -t130 * t184 + t230;
t15 = pkin(4) * t40 + qJD(5) + t27;
t14 = -qJ(5) * t54 + t207;
t13 = qJ(5) * t151 + t164;
t7 = t153 - t205;
t4 = -qJ(5) * t23 + qJD(5) * t151 + t141;
t3 = -qJ(5) * t22 - qJD(5) * t54 + t147;
t2 = -qJ(5) * t9 - qJD(5) * t40 - t148;
t1 = -pkin(4) * t117 + qJ(5) * t8 - qJD(5) * t42 + t229;
t11 = [qJDD(1), t168, t219, qJDD(1) * t122 + 0.2e1 * t128 * t169, 0.2e1 * t128 * t179 - 0.2e1 * t180 * t187, qJDD(2) * t128 + t131 * t134, qJDD(2) * t131 - t128 * t134, 0, t128 * t149 + t131 * t143 + t87, t149 * t131 + (-t143 - t115) * t128, t218 * t128 + t140 * t131 + t87, (t122 + t123) * qJDD(1) * pkin(6) + t139 - t219, -t218 * t131 + (t140 + t115) * t128, pkin(6) * t139 - g(1) * t111 - g(2) * t159 + t45 * t38 + (-t115 + t19) * t225, t151 * t8 + t23 * t42, t151 * t9 - t22 * t42 - t23 * t40 + t54 * t8, t117 * t151 - t118 * t23, t117 * t54 + t118 * t22, 0, g(1) * t35 - g(2) * t37 + t10 * t54 - t117 * t164 - t118 * t141 + t27 * t22 + t25 * t40 + t52 * t9, -g(1) * t34 + g(2) * t36 - t10 * t151 + t207 * t117 + t147 * t118 + t27 * t23 + t25 * t42 - t52 * t8, t1 * t151 + t13 * t8 - t14 * t9 - t2 * t54 - t22 * t7 - t23 * t5 - t3 * t40 - t4 * t42 + t219, t2 * t14 + t7 * t3 + t1 * t13 + t5 * t4 + t137 * (pkin(4) * t54 + t52) + t15 * (pkin(4) * t22 + t25) - g(1) * (t126 * t132 + t111) - g(2) * (t132 * t178 + t95 * t190 + t159) + (-g(1) * (t225 - t178 - t201) - g(2) * t126) * t129; 0, 0, 0, -t175, t187 * t135, t105, t179, qJDD(2), pkin(1) * t193 + t162, g(3) * t128 - t97 + (pkin(1) * t135 + t219) * t131, 0.2e1 * t197 - qJDD(3) + (-t128 * t45 + t131 * t59) * qJD(1) + t162, (-pkin(2) * t128 + qJ(3) * t131) * qJDD(1) + ((t68 - t121) * t128 + (t161 - t66) * t131) * qJD(1), 0.2e1 * t119 + 0.2e1 * t120 + t97 + (qJD(1) * t59 - g(3)) * t128 + (qJD(1) * t45 - t219) * t131, t28 * qJ(3) + t68 * qJD(3) - t32 * pkin(2) - t45 * t59 - g(1) * (-pkin(2) * t194 + t83) - g(2) * (-pkin(2) * t195 + t81) - g(3) * t189 + t152 * qJD(1) * pkin(6), -t211, -t220, t224, t223, t117, -t160 * t117 + t222 * t118 - t33 * t40 - t228, t65 * t117 + t221 * t118 - t33 * t42 - t136, t60 * t8 - t65 * t9 + (-t7 - t209) * t42 + (t5 - t210) * t40, t2 * t65 + t1 * t60 - t15 * (-pkin(4) * t42 + t92) - g(1) * (pkin(4) * t176 + t83) - g(2) * (pkin(4) * t129 * t191 + t81) - g(3) * (t189 + t201) + t210 * t7 + t209 * t5 + (-g(3) * pkin(4) * t127 - t15 * t173 + t219 * (pkin(2) + t95)) * t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(2) - t175, t105, -t122 * t135 - t134, -qJD(2) * t68 + t186 * t45 + t177 + t32, 0, 0, 0, 0, 0, -t117 * t130 - t127 * t158 - t186 * t40, t117 * t127 - t130 * t158 - t186 * t42, -t223 * t127 + t224 * t130, -t15 * t186 + (-t118 * t7 + t1) * t130 + (t118 * t5 + t2) * t127 + t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t211, t220, -t224, -t223, -t117, -t118 * t153 + t228, -t118 * t166 + t136, pkin(4) * t8 - t215 * t40, t215 * t7 + (-t15 * t42 + t219 * t151 + t1 + t214) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39 - t217, t40 * t7 + t42 * t5 + t137 + t168;];
tau_reg = t11;
