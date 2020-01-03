% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRR13
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% tau_reg [5x27]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRR13_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR13_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR13_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR13_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:15:31
% EndTime: 2019-12-31 19:15:38
% DurationCPUTime: 2.75s
% Computational Cost: add. (2059->333), mult. (4143->473), div. (0->0), fcn. (2824->10), ass. (0->178)
t119 = sin(qJ(1));
t123 = cos(qJ(1));
t148 = g(1) * t119 - g(2) * t123;
t118 = sin(qJ(3));
t193 = qJD(1) * t118;
t102 = qJD(4) + t193;
t122 = cos(qJ(3));
t223 = g(3) * t118;
t132 = t122 * t148 - t223;
t124 = -pkin(1) - pkin(6);
t100 = t124 * qJD(1) + qJD(2);
t190 = qJD(3) * t118;
t97 = t124 * qJDD(1) + qJDD(2);
t40 = -qJDD(3) * pkin(3) + t100 * t190 - t122 * t97;
t238 = qJD(4) * pkin(7) * t102 + t132 + t40;
t116 = sin(qJ(5));
t120 = cos(qJ(5));
t117 = sin(qJ(4));
t121 = cos(qJ(4));
t180 = t121 * qJD(3);
t192 = qJD(1) * t122;
t76 = t117 * t192 - t180;
t191 = qJD(3) * t117;
t78 = t121 * t192 + t191;
t147 = t116 * t76 - t120 * t78;
t31 = t116 * t78 + t120 * t76;
t237 = t147 * t31;
t236 = qJDD(2) - t148;
t235 = t147 ^ 2 - t31 ^ 2;
t182 = qJD(5) * t120;
t183 = qJD(5) * t116;
t184 = qJD(4) * t122;
t164 = t117 * t184;
t166 = t118 * t180;
t175 = t122 * qJDD(1);
t24 = t121 * t175 + qJD(4) * t180 + t117 * qJDD(3) + (-t164 - t166) * qJD(1);
t167 = t117 * t190;
t25 = -qJD(1) * t167 + qJD(4) * t78 - t121 * qJDD(3) + t117 * t175;
t6 = -t116 * t25 + t120 * t24 - t76 * t182 - t183 * t78;
t99 = qJD(5) + t102;
t234 = t31 * t99 + t6;
t115 = qJ(4) + qJ(5);
t111 = cos(t115);
t86 = pkin(3) * t118 - pkin(7) * t122 + qJ(2);
t59 = t86 * qJD(1);
t85 = t118 * t100;
t66 = qJD(3) * pkin(7) + t85;
t23 = t117 * t59 + t121 * t66;
t16 = -pkin(8) * t76 + t23;
t14 = t16 * t183;
t222 = g(3) * t122;
t204 = t100 * t122;
t67 = -qJD(3) * pkin(3) - t204;
t39 = pkin(4) * t76 + t67;
t110 = sin(t115);
t201 = t118 * t119;
t48 = t110 * t123 + t111 * t201;
t199 = t118 * t123;
t50 = -t110 * t119 + t111 * t199;
t233 = g(1) * t48 - g(2) * t50 + t111 * t222 + t31 * t39 + t14;
t150 = pkin(3) * t122 + pkin(7) * t118;
t74 = qJD(3) * t150 + qJD(2);
t29 = qJD(1) * t74 + qJDD(1) * t86;
t27 = t121 * t29;
t189 = qJD(3) * t122;
t41 = qJDD(3) * pkin(7) + t100 * t189 + t118 * t97;
t179 = qJD(1) * qJD(3);
t161 = t122 * t179;
t176 = t118 * qJDD(1);
t73 = qJDD(4) + t161 + t176;
t2 = t73 * pkin(4) - t24 * pkin(8) - qJD(4) * t23 - t117 * t41 + t27;
t185 = qJD(4) * t121;
t174 = -t117 * t29 - t121 * t41 - t59 * t185;
t187 = qJD(4) * t117;
t141 = -t187 * t66 - t174;
t3 = -pkin(8) * t25 + t141;
t170 = -t116 * t3 + t120 * t2;
t47 = -t110 * t201 + t111 * t123;
t49 = t110 * t199 + t111 * t119;
t22 = -t117 * t66 + t121 * t59;
t15 = -pkin(8) * t78 + t22;
t13 = pkin(4) * t102 + t15;
t214 = t120 * t16;
t5 = t116 * t13 + t214;
t232 = -g(1) * t47 - g(2) * t49 - qJD(5) * t5 + t110 * t222 + t39 * t147 + t170;
t129 = qJD(5) * t147 - t116 * t24 - t120 * t25;
t231 = -t147 * t99 + t129;
t81 = t116 * t121 + t117 * t120;
t52 = t81 * t122;
t208 = pkin(1) * qJDD(1);
t229 = t208 - t236;
t228 = -t116 * t187 - t117 * t183;
t227 = -qJD(5) * t121 - t185;
t226 = qJD(4) + qJD(5);
t225 = pkin(7) + pkin(8);
t69 = qJDD(5) + t73;
t80 = t116 * t117 - t120 * t121;
t221 = t80 * t69;
t220 = t81 * t69;
t139 = t80 * t118;
t219 = -qJD(1) * t139 - t226 * t80;
t137 = qJD(1) * t81;
t218 = t118 * t137 + t226 * t81;
t197 = t121 * t122;
t82 = t150 * qJD(1);
t217 = t100 * t197 + t117 * t82;
t200 = t118 * t121;
t98 = t124 * t200;
t216 = t117 * t86 + t98;
t215 = t117 * t73;
t213 = t121 * t73;
t212 = t122 * t24;
t211 = t24 * t117;
t210 = t76 * t102;
t209 = t78 * t102;
t126 = qJD(1) ^ 2;
t207 = qJ(2) * t126;
t206 = qJD(3) * t76;
t205 = qJD(3) * t78;
t203 = t117 * t122;
t202 = t117 * t124;
t198 = t119 * t121;
t196 = t121 * t123;
t114 = t122 ^ 2;
t195 = t118 ^ 2 - t114;
t125 = qJD(3) ^ 2;
t194 = -t125 - t126;
t188 = qJD(3) * t124;
t186 = qJD(4) * t118;
t178 = qJDD(1) * qJ(2);
t177 = qJDD(3) * t118;
t165 = t122 * t188;
t173 = t117 * t74 + t121 * t165 + t86 * t185;
t172 = pkin(8) * t200;
t169 = qJD(4) * t225;
t168 = t117 * t193;
t160 = pkin(4) - t202;
t159 = qJD(5) * t13 + t3;
t157 = -qJD(4) * t59 - t41;
t156 = t102 * t124 + t66;
t154 = qJD(1) + t186;
t153 = -t85 + (t168 + t187) * pkin(4);
t65 = t121 * t82;
t93 = t225 * t121;
t152 = qJD(5) * t93 - t100 * t203 + t65 + (pkin(4) * t122 + t172) * qJD(1) + t121 * t169;
t92 = t225 * t117;
t151 = pkin(8) * t168 + qJD(5) * t92 + t117 * t169 + t217;
t149 = g(1) * t123 + g(2) * t119;
t144 = t102 * t185 + t215;
t143 = -t102 * t187 + t213;
t140 = t99 * t80;
t138 = 0.2e1 * qJ(2) * t179 + qJDD(3) * t124;
t136 = -t121 * t184 + t167;
t135 = 0.2e1 * qJD(1) * qJD(2) - t149;
t134 = -pkin(7) * t73 + t102 * t67;
t133 = t148 - t97 + t207;
t130 = t135 + 0.2e1 * t178;
t128 = -t124 * t125 + t130;
t108 = qJDD(3) * t122;
t105 = -pkin(4) * t121 - pkin(3);
t75 = (pkin(4) * t117 - t124) * t122;
t72 = t121 * t86;
t63 = -t117 * t119 + t118 * t196;
t62 = t117 * t199 + t198;
t61 = t117 * t123 + t118 * t198;
t60 = -t117 * t201 + t196;
t56 = t121 * t74;
t53 = t80 * t122;
t42 = -pkin(4) * t136 + t118 * t188;
t38 = -pkin(8) * t203 + t216;
t28 = -pkin(8) * t197 + t118 * t160 + t72;
t12 = -t116 * t166 + t228 * t122 + (t197 * t226 - t167) * t120;
t11 = qJD(3) * t139 - t226 * t52;
t10 = pkin(4) * t25 + t40;
t9 = pkin(8) * t136 - t186 * t202 + t173;
t8 = t56 + (-t98 + (pkin(8) * t122 - t86) * t117) * qJD(4) + (t122 * t160 + t172) * qJD(3);
t4 = -t116 * t16 + t120 * t13;
t1 = [qJDD(1), t148, t149, -0.2e1 * t208 + t236, t130, t229 * pkin(1) + (t135 + t178) * qJ(2), qJDD(1) * t114 - 0.2e1 * t118 * t161, -0.2e1 * t118 * t175 + 0.2e1 * t179 * t195, -t118 * t125 + t108, -t122 * t125 - t177, 0, t118 * t128 + t122 * t138, -t118 * t138 + t122 * t128, -t78 * t164 + (-t190 * t78 + t212) * t121, (t117 * t78 + t121 * t76) * t190 + (-t211 - t121 * t25 + (t117 * t76 - t121 * t78) * qJD(4)) * t122, (-t102 * t180 + t24) * t118 + (t143 + t205) * t122, (t102 * t191 - t25) * t118 + (-t144 - t206) * t122, t102 * t189 + t118 * t73, -g(1) * t63 - g(2) * t61 + t56 * t102 + t72 * t73 + (-t156 * t185 + t188 * t76 + t27) * t118 + (qJD(3) * t22 - t124 * t25 + t185 * t67) * t122 + ((-qJD(4) * t86 - t165) * t102 + t40 * t122 + (-qJD(3) * t67 - t124 * t73 + t157) * t118) * t117, -t173 * t102 - t216 * t73 + g(1) * t62 - g(2) * t60 + (t156 * t187 + (-t121 * t67 + t124 * t78) * qJD(3) + t174) * t118 + (-qJD(3) * t23 + t40 * t121 - t124 * t24 - t187 * t67) * t122, -t11 * t147 - t53 * t6, -t11 * t31 + t12 * t147 - t129 * t53 - t52 * t6, t11 * t99 + t118 * t6 - t147 * t189 - t53 * t69, t118 * t129 - t12 * t99 - t189 * t31 - t52 * t69, t118 * t69 + t189 * t99, (-t116 * t9 + t120 * t8) * t99 + (-t116 * t38 + t120 * t28) * t69 + t170 * t118 + t4 * t189 + t42 * t31 - t75 * t129 + t10 * t52 + t39 * t12 - g(1) * t50 - g(2) * t48 + ((-t116 * t28 - t120 * t38) * t99 - t5 * t118) * qJD(5), -t5 * t189 + g(1) * t49 - g(2) * t47 - t10 * t53 + t39 * t11 + t14 * t118 - t42 * t147 + t75 * t6 + (-(-qJD(5) * t38 + t8) * t99 - t28 * t69 - t2 * t118) * t116 + (-(qJD(5) * t28 + t9) * t99 - t38 * t69 - t159 * t118) * t120; 0, 0, 0, qJDD(1), -t126, -t207 - t229, 0, 0, 0, 0, 0, t118 * t194 + t108, t122 * t194 - t177, 0, 0, 0, 0, 0, -t122 * t25 + (t206 - t215) * t118 + (-t117 * t189 - t121 * t154) * t102, -t212 + (t205 - t213) * t118 + (t117 * t154 - t122 * t180) * t102, 0, 0, 0, 0, 0, qJD(1) * t140 + (-qJD(3) * t81 * t99 + t129) * t122 + ((t120 * t227 - t228) * t99 - t220 + qJD(3) * t31) * t118, t99 * t137 + (qJD(3) * t140 - t6) * t122 + (-(t227 * t116 - t117 * t182 - t120 * t187) * t99 + t221 - qJD(3) * t147) * t118; 0, 0, 0, 0, 0, 0, t122 * t126 * t118, -t195 * t126, t175, -t176, qJDD(3), -t122 * t133 + t223, t118 * t133 + t222, t121 * t209 + t211, (t24 - t210) * t121 + (-t25 - t209) * t117, (t102 * t200 - t122 * t78) * qJD(1) + t144, (-t102 * t117 * t118 + t122 * t76) * qJD(1) + t143, -t102 * t192, -t22 * t192 - t76 * t85 - pkin(3) * t25 - t65 * t102 + (t102 * t204 + t134) * t117 - t238 * t121, -pkin(3) * t24 + t217 * t102 + t238 * t117 + t134 * t121 + t23 * t192 - t78 * t85, -t147 * t219 + t6 * t81, t129 * t81 + t147 * t218 - t219 * t31 - t6 * t80, t147 * t192 + t219 * t99 + t220, t31 * t192 - t218 * t99 - t221, -t99 * t192, (-t116 * t93 - t120 * t92) * t69 - t105 * t129 + t10 * t80 - t4 * t192 + (t116 * t151 - t120 * t152) * t99 + t218 * t39 + t153 * t31 - t132 * t111, -(-t116 * t92 + t120 * t93) * t69 + t105 * t6 + t10 * t81 + t5 * t192 + (t116 * t152 + t120 * t151) * t99 + t219 * t39 - t153 * t147 + t132 * t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78 * t76, -t76 ^ 2 + t78 ^ 2, t24 + t210, t209 - t25, t73, -t66 * t185 - g(1) * t60 - g(2) * t62 + t102 * t23 - t67 * t78 + t27 + (t157 + t222) * t117, g(1) * t61 - g(2) * t63 + g(3) * t197 + t102 * t22 + t67 * t76 - t141, -t237, t235, t234, t231, t69, -(-t116 * t15 - t214) * t99 + (t120 * t69 - t183 * t99 - t31 * t78) * pkin(4) + t232, (-t16 * t99 - t2) * t116 + (t15 * t99 - t159) * t120 + (-t116 * t69 + t147 * t78 - t182 * t99) * pkin(4) + t233; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t237, t235, t234, t231, t69, t5 * t99 + t232, -t116 * t2 - t120 * t159 + t4 * t99 + t233;];
tau_reg = t1;
