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
% Datum: 2019-12-05 18:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:43:00
% EndTime: 2019-12-05 18:43:05
% DurationCPUTime: 1.64s
% Computational Cost: add. (2753->267), mult. (4261->347), div. (0->0), fcn. (3045->14), ass. (0->180)
t258 = qJ(4) + pkin(7);
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
t234 = t161 * t48;
t219 = qJD(5) * t169;
t111 = t119 * qJD(3);
t160 = qJDD(1) + qJDD(2);
t53 = -t162 * t111 + t118 * t160;
t220 = qJD(3) * t174;
t211 = t162 * t220;
t221 = qJD(3) * t170;
t212 = t162 * t221;
t54 = t119 * t160 - t166 * t212 + t167 * t211;
t8 = qJD(5) * t85 + t169 * t53 + t173 * t54 - t98 * t219;
t257 = t8 - t234;
t194 = t169 * t97 + t173 * t98;
t256 = t194 * t48;
t165 = qJ(1) + qJ(2);
t155 = sin(t165);
t156 = cos(t165);
t255 = g(2) * t156 + g(3) * t155;
t250 = -g(2) * t155 + g(3) * t156;
t235 = t161 * t194;
t9 = t194 * qJD(5) + t169 * t54 - t173 * t53;
t254 = -t9 + t235;
t253 = t194 ^ 2 - t48 ^ 2;
t153 = qJ(3) + pkin(9) + qJ(5);
t140 = sin(t153);
t247 = pkin(8) * t97;
t171 = sin(qJ(2));
t245 = pkin(1) * t171;
t215 = qJD(1) * t245;
t208 = t258 * t162 + t215;
t88 = t208 * t174;
t233 = t167 * t88;
t87 = t208 * t170;
t81 = qJD(3) * pkin(3) - t87;
t40 = t166 * t81 + t233;
t23 = t40 + t247;
t141 = cos(t153);
t230 = t141 * t156;
t231 = t141 * t155;
t148 = pkin(3) * t174 + pkin(2);
t175 = cos(qJ(2));
t224 = qJD(1) * t175;
t217 = pkin(1) * t224;
t96 = -t148 * t162 + qJD(4) - t217;
t55 = -pkin(4) * t97 + t96;
t252 = g(1) * t140 - g(2) * t231 + g(3) * t230 + t23 * t219 - t55 * t48;
t154 = t174 * qJD(4);
t209 = qJD(3) * t258;
t108 = -t170 * t209 + t154;
t109 = -qJD(4) * t170 - t174 * t209;
t237 = -t108 * t166 + t167 * t109 + t119 * t217;
t236 = t167 * t108 + t166 * t109 - t118 * t217;
t249 = qJD(5) - t161;
t218 = qJDD(1) * t171;
t222 = qJD(2) * t175;
t102 = pkin(7) * t160 + (qJD(1) * t222 + t218) * pkin(1);
t186 = qJ(4) * t160 + qJD(4) * t162 + t102;
t190 = qJD(3) * t208;
t33 = qJDD(3) * pkin(3) - t186 * t170 - t174 * t190;
t36 = -t170 * t190 + t186 * t174;
t10 = -t166 * t36 + t167 * t33;
t4 = qJDD(3) * pkin(4) - pkin(8) * t54 + t10;
t11 = t166 * t33 + t167 * t36;
t5 = pkin(8) * t53 + t11;
t248 = -g(1) * t141 + t250 * t140 - t169 * t5 + t173 * t4 - t55 * t194;
t246 = pkin(8) * t98;
t244 = pkin(1) * t175;
t243 = pkin(2) * t160;
t242 = pkin(2) * t162;
t241 = pkin(3) * t166;
t112 = t118 * qJD(3);
t240 = pkin(8) * t112;
t239 = pkin(8) * t119;
t238 = g(1) * t174;
t147 = pkin(7) + t245;
t228 = -qJ(4) - t147;
t206 = qJD(3) * t228;
t216 = pkin(1) * t222;
t69 = t170 * t206 + t174 * t216 + t154;
t70 = (-qJD(4) - t216) * t170 + t174 * t206;
t35 = t166 * t70 + t167 * t69;
t75 = t166 * t88;
t42 = -t167 * t87 - t75;
t226 = -qJD(2) * t215 + qJDD(1) * t244;
t101 = -t226 - t243;
t128 = -t217 - t242;
t232 = t101 * t170 + t128 * t220;
t229 = t174 * t160;
t116 = t228 * t170;
t157 = t174 * qJ(4);
t117 = t147 * t174 + t157;
t62 = t166 * t116 + t167 * t117;
t132 = t258 * t170;
t133 = pkin(7) * t174 + t157;
t80 = -t166 * t132 + t167 * t133;
t163 = t170 ^ 2;
t225 = -t174 ^ 2 + t163;
t223 = qJD(2) * t171;
t151 = pkin(3) * t221;
t214 = t128 * t221 + t255 * t174;
t213 = t162 * t223;
t84 = pkin(4) * t111 + t151;
t34 = -t166 * t69 + t167 * t70;
t39 = t167 * t81 - t75;
t41 = t166 * t87 - t233;
t61 = t167 * t116 - t117 * t166;
t79 = -t167 * t132 - t133 * t166;
t207 = -t148 * t156 - t155 * t258;
t192 = t173 * t118 - t119 * t169;
t60 = pkin(3) * t212 - t148 * t160 + qJDD(4) - t226;
t22 = -pkin(4) * t53 + t60;
t64 = t118 * t169 + t119 * t173;
t27 = t64 * qJD(5) + t173 * t111 + t112 * t169;
t205 = g(2) * t230 + g(3) * t231 - t192 * t22 + t55 * t27;
t204 = t162 * t215;
t203 = -t226 - t255;
t202 = t84 - t215;
t115 = t118 * pkin(8);
t59 = t115 + t80;
t201 = qJD(5) * t59 - t237 + t240;
t107 = t111 * pkin(8);
t58 = t79 - t239;
t200 = -qJD(5) * t58 + t107 - t236;
t21 = qJD(3) * pkin(4) - t246 + t39;
t197 = -t169 * t21 - t173 * t23;
t43 = t61 - t239;
t44 = t115 + t62;
t196 = -t169 * t44 + t173 * t43;
t195 = t169 * t43 + t173 * t44;
t193 = -t10 * t119 + t11 * t118 - t40 * t111 - t39 * t112 - t250;
t191 = -t148 * t155 + t156 * t258;
t92 = -pkin(4) * t118 - t148;
t142 = pkin(3) * t167 + pkin(4);
t189 = t142 * t169 + t173 * t241;
t188 = t142 * t173 - t169 * t241;
t185 = -t128 * t162 - t102 + t250;
t26 = t192 * qJD(5) - t111 * t169 + t112 * t173;
t184 = -t140 * t255 + t22 * t64 + t55 * t26;
t177 = qJD(3) ^ 2;
t183 = -pkin(7) * t177 + t204 + t243;
t149 = -pkin(2) - t244;
t182 = -pkin(1) * t213 - t147 * t177 - t149 * t160;
t179 = -pkin(7) * qJDD(3) + (t217 - t242) * qJD(3);
t178 = -qJDD(3) * t147 + (t149 * t162 - t216) * qJD(3);
t176 = cos(qJ(1));
t172 = sin(qJ(1));
t159 = qJDD(3) + qJDD(5);
t158 = t162 ^ 2;
t152 = pkin(1) * t223;
t131 = qJDD(3) * t174 - t170 * t177;
t130 = qJDD(3) * t170 + t174 * t177;
t103 = t160 * t163 + 0.2e1 * t170 * t211;
t83 = t92 - t244;
t74 = -0.2e1 * t225 * t162 * qJD(3) + 0.2e1 * t170 * t229;
t73 = t152 + t84;
t68 = pkin(3) * t162 * t170 + pkin(4) * t98;
t25 = t42 - t246;
t24 = t41 - t247;
t19 = -t107 + t35;
t18 = t34 - t240;
t13 = t159 * t192 - t161 * t27;
t12 = t159 * t64 + t161 * t26;
t2 = t194 * t26 + t64 * t8;
t1 = t192 * t8 - t194 * t27 + t26 * t48 - t64 * t9;
t3 = [qJDD(1), g(2) * t176 + g(3) * t172, -g(2) * t172 + g(3) * t176, t160, (t160 * t175 - t213) * pkin(1) - t203, ((-qJDD(1) - t160) * t171 + (-qJD(1) - t162) * t222) * pkin(1) + t250, t103, t74, t130, t131, 0, t178 * t170 + (-t101 + t182) * t174 + t214, t178 * t174 + (-t182 - t255) * t170 + t232, -t34 * t98 + t35 * t97 + t53 * t62 - t54 * t61 + t193, t11 * t62 + t40 * t35 + t10 * t61 + t39 * t34 + t60 * (-t148 - t244) + t96 * (t152 + t151) - g(2) * (-pkin(1) * t176 + t207) - g(3) * (-pkin(1) * t172 + t191), t2, t1, t12, t13, 0, -t73 * t48 + t83 * t9 + (-qJD(5) * t195 - t169 * t19 + t173 * t18) * t161 + t196 * t159 + t205, t73 * t194 + t83 * t8 - (qJD(5) * t196 + t169 * t18 + t173 * t19) * t161 - t195 * t159 + t184; 0, 0, 0, t160, -t203 + t204, (-t218 + (-qJD(2) + t162) * t224) * pkin(1) + t250, t103, t74, t130, t131, 0, t179 * t170 + (-t101 + t183) * t174 + t214, t179 * t174 + (-t183 - t255) * t170 + t232, t236 * t97 - t237 * t98 + t53 * t80 - t54 * t79 + t193, t11 * t80 + t10 * t79 - t60 * t148 - g(2) * t207 - g(3) * t191 + (-t215 + t151) * t96 + t236 * t40 + t237 * t39, t2, t1, t12, t13, 0, t92 * t9 + (-t169 * t59 + t173 * t58) * t159 - t202 * t48 + (t169 * t200 - t173 * t201) * t161 + t205, t92 * t8 - (t169 * t58 + t173 * t59) * t159 + t202 * t194 + (t169 * t201 + t173 * t200) * t161 + t184; 0, 0, 0, 0, 0, 0, -t170 * t158 * t174, t225 * t158, t170 * t160, t229, qJDD(3), t170 * t185 - t238, g(1) * t170 + t174 * t185, (t40 + t41) * t98 + (t39 - t42) * t97 + (t166 * t53 - t167 * t54) * pkin(3), -t39 * t41 - t40 * t42 + (-t238 + t10 * t167 + t11 * t166 + (-t162 * t96 + t250) * t170) * pkin(3), -t256, t253, t257, t254, t159, t188 * t159 + t68 * t48 - (-t169 * t25 + t173 * t24) * t161 + (-t161 * t189 + t197) * qJD(5) + t248, -t189 * t159 - t173 * t5 - t169 * t4 - t68 * t194 + (t169 * t24 + t173 * t25) * t161 + (-t161 * t188 - t173 * t21) * qJD(5) + t252; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97 ^ 2 - t98 ^ 2, t39 * t98 - t40 * t97 - t255 + t60, 0, 0, 0, 0, 0, t9 + t235, t8 + t234; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t256, t253, t257, t254, t159, t249 * t197 + t248, (-t161 * t23 - t4) * t169 + (-t249 * t21 - t5) * t173 + t252;];
tau_reg = t3;
