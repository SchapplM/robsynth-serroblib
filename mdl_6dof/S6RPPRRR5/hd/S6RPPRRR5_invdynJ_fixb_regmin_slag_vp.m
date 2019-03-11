% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPPRRR5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% 
% Output:
% tau_reg [6x30]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRRR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:29:04
% EndTime: 2019-03-09 02:29:10
% DurationCPUTime: 1.94s
% Computational Cost: add. (2332->306), mult. (4420->383), div. (0->0), fcn. (3032->10), ass. (0->177)
t131 = sin(qJ(5));
t132 = sin(qJ(4));
t135 = cos(qJ(5));
t136 = cos(qJ(4));
t66 = t131 * t136 + t135 * t132;
t62 = t66 * qJD(1);
t237 = qJD(6) + t62;
t118 = qJD(4) + qJD(5);
t130 = sin(qJ(6));
t134 = cos(qJ(6));
t198 = qJD(1) * t132;
t182 = t131 * t198;
t197 = qJD(1) * t136;
t61 = -t135 * t197 + t182;
t44 = -t134 * t118 - t130 * t61;
t238 = t237 * t44;
t190 = qJD(1) * qJD(4);
t180 = t136 * t190;
t186 = t132 * qJDD(1);
t236 = t180 + t186;
t181 = t132 * t190;
t185 = t136 * qJDD(1);
t235 = t181 - t185;
t133 = sin(qJ(1));
t137 = cos(qJ(1));
t161 = g(1) * t137 + g(2) * t133;
t38 = -t61 * pkin(5) + t62 * pkin(9);
t99 = t131 * pkin(4) + pkin(9);
t234 = (pkin(4) * t197 + qJD(6) * t99 + t38) * t237;
t208 = t62 * t118;
t233 = g(1) * t133 - g(2) * t137;
t129 = pkin(1) + qJ(3);
t232 = qJD(1) * t129;
t105 = qJ(2) * qJD(1) + qJD(3);
t79 = -pkin(7) * qJD(1) + t105;
t52 = -pkin(8) * t198 + t132 * t79;
t210 = t135 * t52;
t53 = -pkin(8) * t197 + t136 * t79;
t49 = qJD(4) * pkin(4) + t53;
t31 = t131 * t49 + t210;
t196 = qJD(4) * t132;
t121 = qJD(1) * qJD(2);
t122 = qJ(2) * qJDD(1);
t178 = qJDD(3) + t121 + t122;
t67 = -pkin(7) * qJDD(1) + t178;
t60 = t136 * t67;
t32 = qJDD(4) * pkin(4) + t235 * pkin(8) - t79 * t196 + t60;
t195 = qJD(4) * t136;
t35 = -t236 * pkin(8) + t132 * t67 + t79 * t195;
t231 = t31 * qJD(5) + t131 * t35 - t135 * t32;
t194 = qJD(5) * t131;
t230 = (qJD(5) * t49 + t35) * t135 + t131 * t32 - t52 * t194;
t152 = -qJD(5) * t182 - t235 * t131;
t167 = t118 * t136;
t26 = (qJD(1) * t167 + t186) * t135 + t152;
t22 = qJDD(6) + t26;
t211 = t134 * t22;
t229 = -t130 * t237 ^ 2 + t211;
t128 = -pkin(7) + qJ(2);
t80 = -qJD(2) + t232;
t228 = (qJD(2) + t80 + t232) * qJD(4) + qJDD(4) * t128;
t217 = pkin(8) - t128;
t72 = t217 * t132;
t73 = t217 * t136;
t42 = -t131 * t72 + t135 * t73;
t50 = t136 * qJD(2) + t217 * t196;
t51 = t132 * qJD(2) - qJD(4) * t73;
t14 = -t42 * qJD(5) + t131 * t50 + t135 * t51;
t117 = qJDD(4) + qJDD(5);
t64 = pkin(4) * t198 + t80;
t29 = t62 * pkin(5) + t61 * pkin(9) + t64;
t176 = t117 * pkin(9) + qJD(6) * t29 + t230;
t212 = t131 * t52;
t30 = t135 * t49 - t212;
t23 = -t118 * pkin(5) - t30;
t65 = t131 * t132 - t135 * t136;
t93 = t132 * pkin(4) + t129;
t37 = t66 * pkin(5) + t65 * pkin(9) + t93;
t193 = qJD(5) * t135;
t40 = -t131 * t195 - t132 * t193 - t135 * t196 - t136 * t194;
t43 = -t131 * t73 - t135 * t72;
t7 = -t117 * pkin(5) + t231;
t227 = -(qJD(6) * t37 + t14) * t237 - t176 * t66 - t43 * t22 + t23 * t40 - t7 * t65;
t111 = 0.2e1 * t121;
t126 = qJ(4) + qJ(5);
t107 = sin(t126);
t96 = g(3) * t107;
t108 = cos(t126);
t97 = g(3) * t108;
t224 = t23 * t62;
t223 = t23 * t65;
t222 = t37 * t22;
t221 = t237 * t61;
t220 = t61 * t44;
t46 = t130 * t118 - t134 * t61;
t219 = t61 * t46;
t218 = t61 * t62;
t216 = -t65 * t117 + t40 * t118;
t191 = qJD(6) * t134;
t192 = qJD(6) * t130;
t157 = -t131 * t186 + t135 * t185;
t25 = t157 - t208;
t12 = t130 * t117 + t118 * t191 + t134 * t25 + t61 * t192;
t215 = t12 * t130;
t214 = t130 * t22;
t209 = t61 * t118;
t207 = t133 * t130;
t206 = t133 * t134;
t205 = t137 * t130;
t204 = t137 * t134;
t202 = t137 * pkin(1) + t133 * qJ(2);
t125 = t136 ^ 2;
t200 = t132 ^ 2 - t125;
t138 = qJD(4) ^ 2;
t139 = qJD(1) ^ 2;
t199 = -t138 - t139;
t82 = pkin(4) * t195 + qJD(3);
t120 = qJD(3) * qJD(1);
t188 = qJDD(4) * t132;
t187 = t129 * qJDD(1);
t127 = qJDD(1) * pkin(1);
t184 = t127 - qJDD(2);
t183 = t65 * t192;
t179 = qJDD(2) - t233;
t24 = t118 * pkin(9) + t31;
t119 = qJDD(1) * qJ(3);
t163 = -t119 - t120 - t184;
t47 = t236 * pkin(4) - t163;
t9 = t26 * pkin(5) - t25 * pkin(9) + t47;
t177 = qJD(6) * t24 - t9;
t175 = -t134 * t117 + t130 * t25;
t169 = t134 * t237;
t168 = qJD(6) * t66 + qJD(1);
t166 = -0.2e1 * t180;
t10 = -t130 * t24 + t134 * t29;
t165 = t10 * t61 + t134 * t96 + t23 * t192;
t164 = -t127 + t179;
t33 = t131 * t53 + t210;
t162 = pkin(4) * t194 - t33;
t159 = -t65 * t12 + t40 * t46;
t158 = t22 * t65 - t237 * t40;
t41 = -t131 * t196 - t132 * t194 + t135 * t167;
t155 = -t66 * t117 - t41 * t118;
t154 = -t119 + t164;
t153 = -t176 + t97;
t11 = t130 * t29 + t134 * t24;
t151 = -t11 * t61 + t7 * t130 + t23 * t191 + (g(1) * t205 + g(2) * t207) * t108;
t150 = t161 * t108;
t149 = t111 + 0.2e1 * t122 - t161;
t147 = t80 * qJD(1) + t161;
t146 = t169 * t237 + t214;
t145 = -t150 - t7;
t34 = t135 * t53 - t212;
t144 = -t99 * t22 + t224 + (-pkin(4) * t193 + t34) * t237;
t143 = t161 * t107 + t64 * t62 - t230 + t97;
t142 = -t128 * t138 + t120 - t163 + t187 + t233;
t141 = (-t118 * t197 - t186) * t135 - t152;
t140 = t64 * t61 - t150 - t231 + t96;
t110 = t137 * qJ(2);
t106 = qJDD(4) * t136;
t100 = -t135 * pkin(4) - pkin(5);
t57 = t107 * t204 - t207;
t56 = -t107 * t205 - t206;
t55 = -t107 * t206 - t205;
t54 = t107 * t207 - t204;
t28 = t61 ^ 2 - t62 ^ 2;
t18 = t141 - t209;
t17 = t25 + t208;
t16 = t41 * pkin(5) - t40 * pkin(9) + t82;
t15 = t43 * qJD(5) + t131 * t51 - t135 * t50;
t13 = t46 * qJD(6) + t175;
t8 = t134 * t9;
t4 = t146 + t219;
t3 = -t220 + t229;
t2 = t46 * t169 + t215;
t1 = (t12 - t238) * t134 + (-t237 * t46 - t13) * t130;
t5 = [qJDD(1), t233, t161, -0.2e1 * t127 + t179, t149, t184 * pkin(1) - g(1) * (-t133 * pkin(1) + t110) - g(2) * t202 + (t111 + t122) * qJ(2), qJDD(3) + t149, 0.2e1 * t120 - t154 + t187, -t163 * t129 + t80 * qJD(3) + t178 * qJ(2) + t105 * qJD(2) - g(1) * (-t129 * t133 + t110) - g(2) * (t137 * qJ(3) + t202) t125 * qJDD(1) + t132 * t166, -0.2e1 * t132 * t185 + 0.2e1 * t200 * t190, -t138 * t132 + t106, -t138 * t136 - t188, 0, t142 * t132 + t228 * t136, -t228 * t132 + t142 * t136, -t25 * t65 - t61 * t40, -t25 * t66 + t65 * t26 - t40 * t62 + t61 * t41, t216, t155, 0, t107 * t233 - t42 * t117 - t15 * t118 + t93 * t26 + t64 * t41 + t47 * t66 + t82 * t62, t108 * t233 - t43 * t117 - t14 * t118 + t93 * t25 + t64 * t40 - t47 * t65 - t82 * t61, t134 * t159 + t46 * t183 (-t130 * t46 - t134 * t44) * t40 + (t215 + t13 * t134 + (-t130 * t44 + t134 * t46) * qJD(6)) * t65, t12 * t66 - t134 * t158 + t183 * t237 + t46 * t41, t191 * t237 * t65 - t13 * t66 + t130 * t158 - t44 * t41, t22 * t66 + t237 * t41, -g(1) * t55 - g(2) * t57 + t10 * t41 + t42 * t13 + t15 * t44 + t8 * t66 + (t16 * t237 + t222 + (-t237 * t43 - t24 * t66 - t223) * qJD(6)) * t134 + t227 * t130, -g(1) * t54 - g(2) * t56 - t11 * t41 + t42 * t12 + t15 * t46 + (-(-qJD(6) * t43 + t16) * t237 - t222 + t177 * t66 + qJD(6) * t223) * t130 + t227 * t134; 0, 0, 0, qJDD(1), -t139, -t139 * qJ(2) + t164, -t139, -qJDD(1), -t105 * qJD(1) - t120 + t154, 0, 0, 0, 0, 0, t166 - t186, 0.2e1 * t181 - t185, 0, 0, 0, 0, 0, t141 + t209, -t157 + 0.2e1 * t208, 0, 0, 0, 0, 0, -t220 - t229, t146 - t219; 0, 0, 0, 0, 0, 0, qJDD(1), -t139, -t147 + t178, 0, 0, 0, 0, 0, t199 * t132 + t106, t199 * t136 - t188, 0, 0, 0, 0, 0, -qJD(1) * t62 + t216, qJD(1) * t61 + t155, 0, 0, 0, 0, 0, -t66 * t214 + t65 * t13 - t40 * t44 + (-t130 * t41 - t168 * t134) * t237, -t66 * t211 + (t168 * t130 - t134 * t41) * t237 - t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, t136 * t139 * t132, -t200 * t139, t185, -t186, qJDD(4), g(3) * t132 - t136 * t147 + t60, g(3) * t136 + (t147 - t67) * t132, -t218, t28, t17, t18, t117, t33 * t118 + (t117 * t135 - t118 * t194 - t62 * t197) * pkin(4) + t140, t34 * t118 + (-t117 * t131 - t118 * t193 + t61 * t197) * pkin(4) + t143, t2, t1, t4, t3, t221, t100 * t13 + t162 * t44 + t144 * t130 + (t145 - t234) * t134 + t165, t100 * t12 + t162 * t46 + (-t96 + t234) * t130 + t144 * t134 + t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t218, t28, t17, t18, t117, t31 * t118 + t140, t30 * t118 + t143, t2, t1, t4, t3, t221, -pkin(5) * t13 - t31 * t44 + (-pkin(9) * t22 + t237 * t30 + t224) * t130 + ((-pkin(9) * qJD(6) - t38) * t237 + t145) * t134 + t165, -pkin(5) * t12 + (t130 * t38 + t134 * t30) * t237 - t31 * t46 + t134 * t224 - t130 * t96 + (t192 * t237 - t211) * pkin(9) + t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46 * t44, -t44 ^ 2 + t46 ^ 2, t12 + t238, -t175 + (-qJD(6) + t237) * t46, t22, -g(1) * t56 + g(2) * t54 + t11 * t237 + t130 * t153 - t24 * t191 - t23 * t46 + t8, g(1) * t57 - g(2) * t55 + t10 * t237 + t177 * t130 + t153 * t134 + t23 * t44;];
tau_reg  = t5;
