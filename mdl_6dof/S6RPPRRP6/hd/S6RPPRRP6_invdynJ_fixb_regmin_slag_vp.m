% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% 
% Output:
% tau_reg [6x27]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRRP6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP6_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:11:20
% EndTime: 2019-03-09 02:11:24
% DurationCPUTime: 2.34s
% Computational Cost: add. (2415->362), mult. (4164->434), div. (0->0), fcn. (2360->6), ass. (0->172)
t100 = cos(qJ(4));
t207 = t100 * pkin(8);
t97 = sin(qJ(4));
t137 = pkin(4) * t97 - t207;
t93 = qJDD(1) * pkin(1);
t162 = t93 - qJDD(2);
t87 = qJDD(1) * qJ(3);
t150 = -t87 - t162;
t136 = pkin(4) * t100 + pkin(8) * t97;
t49 = t136 * qJD(4) + qJD(3);
t16 = t49 * qJD(1) + t137 * qJDD(1) - t150;
t99 = cos(qJ(5));
t172 = qJD(5) * t99;
t167 = qJD(4) * t100;
t88 = qJD(1) * qJD(2);
t89 = qJ(2) * qJDD(1);
t149 = qJDD(3) + t88 + t89;
t55 = -pkin(7) * qJDD(1) + t149;
t74 = qJ(2) * qJD(1) + qJD(3);
t64 = -pkin(7) * qJD(1) + t74;
t24 = qJDD(4) * pkin(8) + t64 * t167 + t55 * t97;
t95 = pkin(1) + qJ(3);
t59 = t137 + t95;
t33 = t59 * qJD(1) - qJD(2);
t96 = sin(qJ(5));
t157 = -t96 * t16 - t33 * t172 - t99 * t24;
t174 = qJD(5) * t96;
t58 = t97 * t64;
t39 = qJD(4) * pkin(8) + t58;
t119 = -t39 * t174 - t157;
t161 = qJD(1) * qJD(4);
t147 = t100 * t161;
t163 = t97 * qJDD(1);
t48 = qJDD(5) + t147 + t163;
t187 = qJ(6) * t48;
t68 = qJD(1) * t97 + qJD(5);
t1 = qJD(6) * t68 + t119 + t187;
t11 = t33 * t99 - t39 * t96;
t170 = qJD(6) - t11;
t8 = -pkin(5) * t68 + t170;
t12 = t33 * t96 + t39 * t99;
t9 = qJ(6) * t68 + t12;
t138 = t8 * t96 + t9 * t99;
t101 = cos(qJ(1));
t98 = sin(qJ(1));
t190 = g(1) * t98 - g(2) * t101;
t146 = t99 * t16 - t39 * t172 - t33 * t174 - t96 * t24;
t214 = pkin(5) * t48;
t2 = qJDD(6) - t146 - t214;
t227 = -t138 * qJD(5) - t1 * t96 + t2 * t99 - t190;
t148 = t97 * t161;
t159 = t100 * qJDD(1);
t168 = qJD(1) * t100;
t177 = qJD(4) * t96;
t53 = t99 * t168 + t177;
t175 = qJD(5) * t53;
t18 = -t99 * qJDD(4) + (-t148 + t159) * t96 + t175;
t203 = t53 * t68;
t226 = t203 - t18;
t134 = g(1) * t101 + g(2) * t98;
t120 = t134 * t100;
t211 = g(3) * t97;
t225 = t211 - t120;
t173 = qJD(5) * t97;
t144 = qJD(1) + t173;
t152 = t99 * t167;
t171 = t99 * qJD(4);
t153 = t97 * t171;
t17 = (t100 * t174 + t153) * qJD(1) - qJD(5) * t171 - t96 * qJDD(4) - t99 * t159;
t178 = qJD(4) * t53;
t196 = t99 * t48;
t224 = (t144 * t96 - t152) * t68 + (t178 - t196) * t97 + t100 * t17;
t127 = t68 * t99;
t51 = t96 * t168 - t171;
t223 = t51 * t127;
t222 = qJD(1) * t95;
t184 = t100 * t64;
t40 = -qJD(4) * pkin(4) - t184;
t10 = pkin(5) * t51 - qJ(6) * t53 + t40;
t213 = pkin(8) * t48;
t221 = t68 * t10 - t213;
t201 = t96 * t48;
t122 = t68 * t172 + t201;
t199 = t97 * t99;
t158 = t68 * t199;
t185 = t100 * t53;
t218 = qJD(1) * (t158 + t185) + t122;
t65 = -qJD(2) + t222;
t94 = -pkin(7) + qJ(2);
t217 = (qJD(2) + t65 + t222) * qJD(4) + qJDD(4) * t94;
t216 = t53 ^ 2;
t215 = t68 ^ 2;
t82 = 0.2e1 * t88;
t210 = t68 * t9;
t208 = g(3) * t100;
t206 = t12 * t68;
t205 = t17 * t96;
t15 = t17 * t99;
t204 = t53 * t51;
t202 = t94 * t96;
t200 = t96 * t97;
t198 = t98 * t96;
t197 = t98 * t99;
t195 = t99 * t59;
t132 = pkin(5) * t96 - qJ(6) * t99;
t194 = -t96 * qJD(6) + t68 * t132 - t58;
t183 = t100 * t99;
t57 = t136 * qJD(1);
t193 = t64 * t183 + t96 * t57;
t192 = t94 * t199 + t96 * t59;
t191 = t101 * pkin(1) + t98 * qJ(2);
t92 = t100 ^ 2;
t189 = t97 ^ 2 - t92;
t188 = pkin(8) * qJD(5);
t182 = t101 * t97;
t181 = t101 * t99;
t179 = qJD(2) * t97;
t176 = qJD(4) * t97;
t102 = qJD(4) ^ 2;
t103 = qJD(1) ^ 2;
t169 = -t102 - t103;
t166 = qJDD(1) * t95;
t164 = qJDD(4) * t97;
t160 = qJD(3) * qJD(1);
t155 = t101 * qJ(3) + t191;
t154 = t68 * t188;
t151 = qJDD(2) - t190;
t145 = t94 * t152 + t59 * t172 + t99 * t179 + t96 * t49;
t143 = -0.2e1 * t147;
t142 = -t93 + t151;
t41 = t97 * t198 - t181;
t43 = t96 * t182 + t197;
t141 = -g(1) * t41 + g(2) * t43;
t42 = t101 * t96 + t97 * t197;
t44 = t97 * t181 - t198;
t140 = g(1) * t42 - g(2) * t44;
t139 = t8 * t99 - t9 * t96;
t133 = pkin(5) * t99 + qJ(6) * t96;
t130 = -t87 + t142;
t126 = pkin(4) + t133;
t125 = t132 - t94;
t121 = -t68 * t174 + t196;
t23 = -qJDD(4) * pkin(4) - t100 * t55 + t64 * t176;
t118 = -t134 + t82 + 0.2e1 * t89;
t116 = qJD(1) * t65 + t134;
t115 = t68 * t40 - t213;
t114 = t116 - t55;
t113 = g(1) * t43 + g(2) * t41 + t96 * t208 + t146;
t112 = -t134 * t97 - t208;
t111 = t51 * t168 + t215 * t96 - t196;
t3 = pkin(5) * t18 + qJ(6) * t17 - qJD(6) * t53 + t23;
t110 = -t120 - t3 - t154;
t109 = t139 * qJD(5) + t1 * t99 + t2 * t96;
t107 = t10 * t53 + qJDD(6) - t113;
t56 = -t150 + t160;
t106 = -t102 * t94 + t160 + t166 + t190 + t56;
t105 = -g(1) * t44 - g(2) * t42 - g(3) * t183 + t119;
t104 = -t48 * t200 - t100 * t18 + t51 * t176 + (-t144 * t99 - t96 * t167) * t68;
t81 = t101 * qJ(2);
t77 = qJDD(4) * t100;
t69 = g(3) * t199;
t31 = t125 * t100;
t28 = -t195 + (-pkin(5) + t202) * t97;
t27 = qJ(6) * t97 + t192;
t25 = pkin(5) * t53 + qJ(6) * t51;
t22 = -t99 * t57 + (-pkin(5) * qJD(1) + t64 * t96) * t100;
t19 = qJ(6) * t168 + t193;
t7 = -t125 * t176 + (t133 * qJD(5) - qJD(6) * t99 - qJD(2)) * t100;
t6 = t51 * t68 - t17;
t5 = -pkin(5) * t167 + (t94 * t173 - t49) * t99 + (qJD(5) * t59 + t94 * t167 + t179) * t96;
t4 = qJ(6) * t167 + (-t94 * t174 + qJD(6)) * t97 + t145;
t13 = [qJDD(1), t190, t134, -0.2e1 * t93 + t151, t118, t162 * pkin(1) - g(1) * (-pkin(1) * t98 + t81) - g(2) * t191 + (t82 + t89) * qJ(2), qJDD(3) + t118, -t130 + 0.2e1 * t160 + t166, t56 * t95 + t65 * qJD(3) + t149 * qJ(2) + t74 * qJD(2) - g(1) * (-t95 * t98 + t81) - g(2) * t155, qJDD(1) * t92 + t97 * t143, -0.2e1 * t97 * t159 + 0.2e1 * t189 * t161, -t102 * t97 + t77, -t100 * t102 - t164, 0, t217 * t100 + t106 * t97, t106 * t100 - t217 * t97, -t53 * t153 + (-t53 * t174 - t15) * t100 (t51 * t99 + t53 * t96) * t176 + (t205 - t18 * t99 + (t51 * t96 - t53 * t99) * qJD(5)) * t100 (-t68 * t171 - t17) * t97 + (t121 + t178) * t100 (t68 * t177 - t18) * t97 + (-qJD(4) * t51 - t122) * t100, t68 * t167 + t48 * t97 (-t59 * t174 + t99 * t49) * t68 + t48 * t195 + ((-qJD(2) * t96 - t94 * t172) * t68 - t94 * t201 + (-t40 * t96 + t51 * t94) * qJD(4) + t146) * t97 + (t40 * t172 - qJD(2) * t51 - t94 * t18 + t23 * t96 + (-t68 * t202 + t11) * qJD(4)) * t100 + t140, -t145 * t68 - t192 * t48 + ((t68 * t94 + t39) * t174 + (-t40 * t99 + t53 * t94) * qJD(4) + t157) * t97 + (-qJD(2) * t53 - qJD(4) * t12 + t17 * t94 - t40 * t174 + t23 * t99) * t100 + t141, t18 * t31 - t28 * t48 - t5 * t68 + t51 * t7 + (-t10 * t177 - t2) * t97 + (-qJD(4) * t8 + t10 * t172 + t3 * t96) * t100 + t140, t227 * t100 - t139 * t176 - t17 * t28 - t18 * t27 - t4 * t51 + t5 * t53, t17 * t31 + t27 * t48 + t4 * t68 - t53 * t7 + (t10 * t171 + t1) * t97 + (qJD(4) * t9 + t10 * t174 - t3 * t99) * t100 - t141, t1 * t27 + t9 * t4 + t3 * t31 + t10 * t7 + t2 * t28 + t8 * t5 - g(1) * (-pkin(5) * t42 - pkin(7) * t101 - qJ(6) * t41 + t81) - g(2) * (pkin(4) * t182 + pkin(5) * t44 + qJ(6) * t43 - t101 * t207 + t155) + (g(2) * pkin(7) + g(1) * t59) * t98; 0, 0, 0, qJDD(1), -t103, -qJ(2) * t103 + t142, -t103, -qJDD(1) (-qJD(3) - t74) * qJD(1) + t130, 0, 0, 0, 0, 0, t143 - t163, 0.2e1 * t148 - t159, 0, 0, 0, 0, 0, t111, t218, t111, -t226 * t96 - t15 + t223, -t218 (t10 * t100 - t138 * t97) * qJD(1) + t227; 0, 0, 0, 0, 0, 0, qJDD(1), -t103, -t116 + t149, 0, 0, 0, 0, 0, t169 * t97 + t77, t169 * t100 - t164, 0, 0, 0, 0, 0, t104, t224, t104 (t144 * t53 - t51 * t167 - t18 * t97) * t99 + (t144 * t51 + t53 * t167 - t17 * t97) * t96, -t224, t139 * qJD(1) + (qJD(4) * t138 - t3) * t100 + (qJD(4) * t10 + t109) * t97 - t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, t100 * t103 * t97, -t189 * t103, t159, -t163, qJDD(4), -t114 * t100 + t211, t114 * t97 + t208, t53 * t127 - t205, -t15 - t223 + (-t18 - t203) * t96 (t158 - t185) * qJD(1) + t122 (t100 * t51 - t68 * t200) * qJD(1) + t121, -t68 * t168, -t11 * t168 - t51 * t58 - pkin(4) * t18 + t69 + (t68 * t184 + t115) * t96 + (-t23 + (-t57 - t188) * t68 - t120) * t99, pkin(4) * t17 + t193 * t68 + t12 * t168 - t53 * t58 + t115 * t99 + (t154 + t23 - t225) * t96, t110 * t99 - t126 * t18 + t8 * t168 + t194 * t51 + t22 * t68 + t221 * t96 + t69, t19 * t51 - t22 * t53 + (t1 + t68 * t8 + (-t18 + t175) * pkin(8)) * t99 + (t2 - t210 + (qJD(5) * t51 - t17) * pkin(8)) * t96 + t112, -t9 * t168 - t17 * t126 - t19 * t68 - t194 * t53 - t221 * t99 + (t110 + t211) * t96, -t9 * t19 - t8 * t22 + t194 * t10 + (t109 + t112) * pkin(8) + (-t3 + t225) * t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t204, -t51 ^ 2 + t216, t6, t226, t48, -t40 * t53 + t113 + t206, t11 * t68 + t40 * t51 - t105, -t25 * t51 - t107 + t206 + 0.2e1 * t214, pkin(5) * t17 - qJ(6) * t18 + (-t12 + t9) * t53 + (t8 - t170) * t51, 0.2e1 * t187 - t10 * t51 + t25 * t53 + (0.2e1 * qJD(6) - t11) * t68 + t105, t1 * qJ(6) - t2 * pkin(5) - t10 * t25 - t8 * t12 - g(1) * (-pkin(5) * t43 + qJ(6) * t44) - g(2) * (-pkin(5) * t41 + qJ(6) * t42) + t170 * t9 + t132 * t208; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48 + t204, t6, -t215 - t216, t107 - t210 - t214;];
tau_reg  = t13;
