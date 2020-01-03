% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:09:31
% EndTime: 2020-01-03 12:09:36
% DurationCPUTime: 1.64s
% Computational Cost: add. (2753->267), mult. (4261->346), div. (0->0), fcn. (3045->14), ass. (0->181)
t171 = sin(qJ(2));
t246 = pkin(1) * t171;
t213 = qJD(1) * t246;
t175 = cos(qJ(2));
t245 = pkin(1) * t175;
t224 = -qJD(2) * t213 + qJDD(1) * t245;
t160 = qJDD(1) + qJDD(2);
t244 = pkin(2) * t160;
t101 = -t224 - t244;
t165 = qJ(1) + qJ(2);
t155 = sin(t165);
t156 = cos(t165);
t198 = g(2) * t156 + g(3) * t155;
t260 = t101 + t198;
t259 = qJ(4) + pkin(7);
t161 = qJD(3) + qJD(5);
t169 = sin(qJ(5));
t173 = cos(qJ(5));
t166 = sin(pkin(9));
t167 = cos(pkin(9));
t170 = sin(qJ(3));
t174 = cos(qJ(3));
t118 = -t166 * t170 + t167 * t174;
t162 = qJD(1) + qJD(2);
t97 = t118 * t162;
t85 = t173 * t97;
t119 = t166 * t174 + t167 * t170;
t98 = t119 * t162;
t48 = -t169 * t98 + t85;
t232 = t161 * t48;
t217 = qJD(5) * t169;
t111 = t119 * qJD(3);
t53 = -t111 * t162 + t118 * t160;
t218 = qJD(3) * t174;
t210 = t162 * t218;
t219 = qJD(3) * t170;
t211 = t162 * t219;
t54 = t119 * t160 - t166 * t211 + t167 * t210;
t8 = qJD(5) * t85 + t169 * t53 + t173 * t54 - t98 * t217;
t258 = t8 - t232;
t193 = t169 * t97 + t173 * t98;
t257 = t193 * t48;
t252 = g(2) * t155 - g(3) * t156;
t233 = t161 * t193;
t9 = qJD(5) * t193 + t169 * t54 - t173 * t53;
t256 = -t9 + t233;
t255 = t193 ^ 2 - t48 ^ 2;
t153 = qJ(3) + pkin(9) + qJ(5);
t142 = sin(t153);
t143 = cos(t153);
t248 = pkin(8) * t97;
t207 = t259 * t162 + t213;
t88 = t207 * t174;
t231 = t167 * t88;
t87 = t207 * t170;
t81 = qJD(3) * pkin(3) - t87;
t40 = t166 * t81 + t231;
t23 = t40 + t248;
t148 = pkin(3) * t174 + pkin(2);
t222 = qJD(1) * t175;
t215 = pkin(1) * t222;
t96 = -t148 * t162 + qJD(4) - t215;
t55 = -pkin(4) * t97 + t96;
t254 = g(1) * t142 + t252 * t143 + t23 * t217 - t55 * t48;
t154 = t174 * qJD(4);
t208 = qJD(3) * t259;
t108 = -t170 * t208 + t154;
t109 = -qJD(4) * t170 - t174 * t208;
t235 = -t108 * t166 + t167 * t109 + t119 * t215;
t234 = t167 * t108 + t166 * t109 - t118 * t215;
t177 = qJD(3) ^ 2;
t251 = pkin(7) * t177 - t244;
t250 = qJD(5) - t161;
t229 = t142 * t156;
t230 = t142 * t155;
t216 = qJDD(1) * t171;
t220 = qJD(2) * t175;
t102 = pkin(7) * t160 + (qJD(1) * t220 + t216) * pkin(1);
t186 = qJ(4) * t160 + qJD(4) * t162 + t102;
t190 = qJD(3) * t207;
t33 = qJDD(3) * pkin(3) - t170 * t186 - t174 * t190;
t36 = -t170 * t190 + t174 * t186;
t10 = -t166 * t36 + t167 * t33;
t4 = qJDD(3) * pkin(4) - pkin(8) * t54 + t10;
t11 = t166 * t33 + t167 * t36;
t5 = pkin(8) * t53 + t11;
t249 = -g(1) * t143 + g(2) * t230 - g(3) * t229 - t169 * t5 + t173 * t4 - t55 * t193;
t247 = pkin(8) * t98;
t243 = pkin(2) * t162;
t242 = pkin(3) * t166;
t112 = t118 * qJD(3);
t240 = pkin(8) * t112;
t239 = pkin(8) * t119;
t238 = g(1) * t174;
t147 = pkin(7) + t246;
t227 = -qJ(4) - t147;
t205 = qJD(3) * t227;
t214 = pkin(1) * t220;
t69 = t170 * t205 + t174 * t214 + t154;
t70 = (-qJD(4) - t214) * t170 + t174 * t205;
t35 = t166 * t70 + t167 * t69;
t75 = t166 * t88;
t42 = -t167 * t87 - t75;
t228 = t174 * t160;
t116 = t227 * t170;
t157 = t174 * qJ(4);
t117 = t147 * t174 + t157;
t62 = t166 * t116 + t167 * t117;
t134 = t259 * t170;
t135 = pkin(7) * t174 + t157;
t80 = -t166 * t134 + t167 * t135;
t226 = t155 * t148 - t156 * t259;
t163 = t170 ^ 2;
t223 = -t174 ^ 2 + t163;
t221 = qJD(2) * t171;
t151 = pkin(3) * t219;
t212 = t162 * t221;
t84 = pkin(4) * t111 + t151;
t34 = -t166 * t69 + t167 * t70;
t39 = t167 * t81 - t75;
t41 = t166 * t87 - t231;
t61 = t167 * t116 - t117 * t166;
t79 = -t167 * t134 - t135 * t166;
t206 = t156 * t148 + t155 * t259;
t60 = pkin(3) * t211 - t148 * t160 + qJDD(4) - t224;
t22 = -pkin(4) * t53 + t60;
t191 = t173 * t118 - t119 * t169;
t26 = qJD(5) * t191 - t111 * t169 + t112 * t173;
t64 = t118 * t169 + t119 * t173;
t204 = g(2) * t229 + g(3) * t230 + t22 * t64 + t55 * t26;
t203 = t162 * t213;
t128 = -t215 - t243;
t202 = t128 * t218 + t260 * t170;
t201 = t84 - t213;
t115 = t118 * pkin(8);
t59 = t115 + t80;
t200 = qJD(5) * t59 - t235 + t240;
t107 = t111 * pkin(8);
t58 = t79 - t239;
t199 = -qJD(5) * t58 + t107 - t234;
t21 = qJD(3) * pkin(4) - t247 + t39;
t196 = -t169 * t21 - t173 * t23;
t43 = t61 - t239;
t44 = t115 + t62;
t195 = -t169 * t44 + t173 * t43;
t194 = t169 * t43 + t173 * t44;
t192 = -t10 * t119 + t11 * t118 - t40 * t111 - t39 * t112 - t252;
t92 = -pkin(4) * t118 - t148;
t144 = pkin(3) * t167 + pkin(4);
t189 = t144 * t169 + t173 * t242;
t188 = t144 * t173 - t169 * t242;
t185 = -t128 * t162 - t102 + t252;
t27 = qJD(5) * t64 + t173 * t111 + t112 * t169;
t184 = -t143 * t198 - t191 * t22 + t55 * t27;
t183 = -t198 + t203;
t149 = -pkin(2) - t245;
t182 = pkin(1) * t212 + t147 * t177 + t149 * t160;
t179 = -pkin(7) * qJDD(3) + (t215 - t243) * qJD(3);
t178 = -qJDD(3) * t147 + (t149 * t162 - t214) * qJD(3);
t176 = cos(qJ(1));
t172 = sin(qJ(1));
t159 = qJDD(3) + qJDD(5);
t158 = t162 ^ 2;
t152 = pkin(1) * t221;
t133 = qJDD(3) * t174 - t170 * t177;
t132 = qJDD(3) * t170 + t174 * t177;
t113 = t128 * t219;
t103 = t160 * t163 + 0.2e1 * t170 * t210;
t83 = t92 - t245;
t74 = -0.2e1 * t223 * t162 * qJD(3) + 0.2e1 * t170 * t228;
t73 = t152 + t84;
t68 = pkin(3) * t162 * t170 + pkin(4) * t98;
t25 = t42 - t247;
t24 = t41 - t248;
t19 = -t107 + t35;
t18 = t34 - t240;
t13 = t159 * t191 - t161 * t27;
t12 = t159 * t64 + t161 * t26;
t2 = t193 * t26 + t64 * t8;
t1 = t191 * t8 - t193 * t27 + t26 * t48 - t64 * t9;
t3 = [qJDD(1), -g(2) * t176 - g(3) * t172, g(2) * t172 - g(3) * t176, t160, (t160 * t175 - t212) * pkin(1) - t198 + t224, ((-qJDD(1) - t160) * t171 + (-qJD(1) - t162) * t220) * pkin(1) + t252, t103, t74, t132, t133, 0, t113 + t178 * t170 + (-t182 - t260) * t174, t170 * t182 + t174 * t178 + t202, -t34 * t98 + t35 * t97 + t53 * t62 - t54 * t61 + t192, t11 * t62 + t40 * t35 + t10 * t61 + t39 * t34 + t60 * (-t148 - t245) + t96 * (t152 + t151) - g(2) * (pkin(1) * t176 + t206) - g(3) * (pkin(1) * t172 + t226), t2, t1, t12, t13, 0, -t73 * t48 + t83 * t9 + (-qJD(5) * t194 - t169 * t19 + t173 * t18) * t161 + t195 * t159 + t184, t73 * t193 + t83 * t8 - (qJD(5) * t195 + t169 * t18 + t173 * t19) * t161 - t194 * t159 + t204; 0, 0, 0, t160, t183 + t224, (-t216 + (-qJD(2) + t162) * t222) * pkin(1) + t252, t103, t74, t132, t133, 0, t113 + t179 * t170 + (-t101 + t183 - t251) * t174, t179 * t174 + (-t203 + t251) * t170 + t202, t234 * t97 - t235 * t98 + t53 * t80 - t54 * t79 + t192, t11 * t80 + t10 * t79 - t60 * t148 - g(2) * t206 - g(3) * t226 + (-t213 + t151) * t96 + t234 * t40 + t235 * t39, t2, t1, t12, t13, 0, t92 * t9 + (-t169 * t59 + t173 * t58) * t159 - t201 * t48 + (t169 * t199 - t173 * t200) * t161 + t184, t92 * t8 - (t169 * t58 + t173 * t59) * t159 + t201 * t193 + (t169 * t200 + t173 * t199) * t161 + t204; 0, 0, 0, 0, 0, 0, -t170 * t158 * t174, t223 * t158, t170 * t160, t228, qJDD(3), t170 * t185 - t238, g(1) * t170 + t174 * t185, (t40 + t41) * t98 + (t39 - t42) * t97 + (t166 * t53 - t167 * t54) * pkin(3), -t39 * t41 - t40 * t42 + (-t238 + t10 * t167 + t11 * t166 + (-t162 * t96 + t252) * t170) * pkin(3), -t257, t255, t258, t256, t159, t188 * t159 + t68 * t48 - (-t169 * t25 + t173 * t24) * t161 + (-t161 * t189 + t196) * qJD(5) + t249, -t189 * t159 - t173 * t5 - t169 * t4 - t68 * t193 + (t169 * t24 + t173 * t25) * t161 + (-t161 * t188 - t173 * t21) * qJD(5) + t254; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97 ^ 2 - t98 ^ 2, t39 * t98 - t40 * t97 + t198 + t60, 0, 0, 0, 0, 0, t9 + t233, t8 + t232; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t257, t255, t258, t256, t159, t196 * t250 + t249, (-t161 * t23 - t4) * t169 + (-t21 * t250 - t5) * t173 + t254;];
tau_reg = t3;
