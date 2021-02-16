% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRPP4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% tau_reg [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:26
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:24:49
% EndTime: 2021-01-15 22:24:59
% DurationCPUTime: 2.29s
% Computational Cost: add. (4271->328), mult. (10101->406), div. (0->0), fcn. (7093->12), ass. (0->186)
t158 = qJ(2) + qJ(3);
t149 = sin(t158);
t150 = cos(t158);
t163 = sin(qJ(1));
t166 = cos(qJ(1));
t204 = g(1) * t166 + g(2) * t163;
t265 = -g(3) * t150 + t149 * t204;
t155 = qJD(2) + qJD(3);
t159 = sin(pkin(8));
t160 = cos(pkin(8));
t164 = cos(qJ(3));
t165 = cos(qJ(2));
t226 = qJD(1) * t165;
t214 = t164 * t226;
t161 = sin(qJ(3));
t162 = sin(qJ(2));
t227 = qJD(1) * t162;
t215 = t161 * t227;
t85 = -t214 + t215;
t87 = -t161 * t226 - t164 * t227;
t56 = t159 * t87 - t160 * t85;
t264 = t56 * t155;
t223 = qJD(1) * qJD(2);
t211 = t165 * t223;
t222 = t162 * qJDD(1);
t262 = t211 + t222;
t195 = -t159 * t85 - t160 * t87;
t261 = t195 ^ 2;
t256 = pkin(6) + pkin(7);
t151 = t165 * pkin(2);
t246 = pkin(1) + t151;
t114 = t256 * t162;
t100 = qJD(1) * t114;
t115 = t256 * t165;
t102 = qJD(1) * t115;
t231 = t164 * t102;
t193 = t161 * t100 - t231;
t235 = t85 * qJ(4);
t183 = t193 + t235;
t220 = t159 * t161 * pkin(2);
t224 = qJD(3) * t164;
t88 = t161 * t102;
t241 = -t164 * t100 - t88;
t81 = t87 * qJ(4);
t50 = t81 + t241;
t242 = -qJD(3) * t220 - t159 * t183 + (pkin(2) * t224 - t50) * t160;
t239 = qJD(2) * pkin(2);
t93 = -t100 + t239;
t206 = t164 * t93 - t88;
t46 = t206 + t81;
t148 = pkin(8) + t158;
t135 = sin(t148);
t136 = cos(t148);
t197 = t136 * pkin(4) + t135 * qJ(5);
t260 = g(1) * t163 - g(2) * t166;
t153 = qJDD(2) + qJDD(3);
t259 = -t153 * pkin(4) + qJDD(5);
t194 = -t161 * t93 - t231;
t66 = qJDD(2) * pkin(2) - t256 * t262;
t212 = t162 * t223;
t221 = t165 * qJDD(1);
t68 = t256 * (-t212 + t221);
t175 = qJD(3) * t194 - t161 * t68 + t164 * t66;
t48 = qJD(3) * t214 - t155 * t215 + t161 * t221 + t262 * t164;
t13 = t153 * pkin(3) - t48 * qJ(4) + t87 * qJD(4) + t175;
t225 = qJD(3) * t161;
t257 = (qJD(3) * t93 + t68) * t164 - t102 * t225 + t161 * t66;
t196 = t161 * t222 - t164 * t221;
t96 = t161 * t165 + t164 * t162;
t65 = t155 * t96;
t49 = qJD(1) * t65 + t196;
t16 = -t49 * qJ(4) - t85 * qJD(4) + t257;
t245 = -t160 * t13 + t159 * t16;
t190 = -g(3) * t136 + t204 * t135 - t245;
t113 = t246 * qJD(1);
t67 = t85 * pkin(3) + qJD(4) - t113;
t26 = -pkin(4) * t56 - qJ(5) * t195 + t67;
t173 = -t195 * t26 + t190 - t259;
t192 = -t164 * t114 - t161 * t115;
t182 = -t96 * qJ(4) + t192;
t240 = -t161 * t114 + t164 * t115;
t95 = t161 * t162 - t164 * t165;
t53 = -t95 * qJ(4) + t240;
t34 = t159 * t182 + t160 * t53;
t217 = qJD(2) * t256;
t101 = t162 * t217;
t103 = t165 * t217;
t172 = -t240 * qJD(3) + t161 * t101 - t164 * t103;
t64 = t155 * t95;
t169 = t64 * qJ(4) - t96 * qJD(4) + t172;
t188 = -t164 * t101 - t161 * t103 - t114 * t224 - t115 * t225;
t27 = -t65 * qJ(4) - t95 * qJD(4) + t188;
t8 = t159 * t169 + t160 * t27;
t258 = t135 * t260 + t34 * t153 + t8 * t155;
t255 = t87 * pkin(3);
t254 = pkin(3) * t149;
t253 = pkin(4) * t135;
t47 = -t194 - t235;
t42 = t160 * t47;
t21 = t159 * t46 + t42;
t248 = t21 * t195;
t247 = t87 * t85;
t4 = t159 * t13 + t160 * t16;
t41 = t155 * pkin(3) + t46;
t20 = t159 * t41 + t42;
t244 = qJD(5) + t242;
t232 = t160 * t161;
t243 = -t159 * t50 + t160 * t183 + (t159 * t164 + t232) * qJD(3) * pkin(2);
t238 = t159 * t47;
t237 = t21 * t155;
t22 = t160 * t46 - t238;
t236 = t22 * t155;
t234 = t136 * t163;
t233 = t136 * t166;
t230 = qJD(5) - t22;
t141 = t164 * pkin(2) + pkin(3);
t80 = pkin(2) * t232 + t159 * t141;
t139 = pkin(3) * t150;
t229 = t139 + t151;
t156 = t162 ^ 2;
t228 = -t165 ^ 2 + t156;
t145 = t162 * t239;
t219 = t139 + t197;
t216 = t243 * t195;
t59 = t65 * pkin(3) + t145;
t213 = t243 * t155;
t104 = -t162 * pkin(2) - t254;
t210 = t104 - t253;
t23 = t159 * t48 + t160 * t49;
t70 = t95 * pkin(3) - t246;
t205 = -t253 - t254;
t19 = t160 * t41 - t238;
t17 = -t155 * pkin(4) + qJD(5) - t19;
t18 = t155 * qJ(5) + t20;
t202 = -t17 * t56 + t18 * t195;
t201 = t19 * t56 + t195 * t20;
t200 = -t56 ^ 2 - t261;
t199 = g(1) * t233 + g(2) * t234 + g(3) * t135 - t4;
t24 = -t159 * t49 + t160 * t48;
t191 = t155 * t195 + t23;
t79 = t160 * t141 - t220;
t189 = -0.2e1 * pkin(1) * t223 - pkin(6) * qJDD(2);
t82 = pkin(2) * t212 - qJDD(1) * t246;
t31 = pkin(4) * t195 - t56 * qJ(5) - t255;
t185 = -t67 * t56 + t199;
t184 = t24 + t264;
t143 = t153 * qJ(5);
t181 = t26 * t56 + t143 - t199;
t180 = -t195 * t67 + t190;
t33 = t159 * t53 - t160 * t182;
t7 = t159 * t27 - t160 * t169;
t178 = g(1) * t234 - g(2) * t233 - t33 * t153 - t7 * t155;
t167 = qJD(2) ^ 2;
t177 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t167 + t260;
t168 = qJD(1) ^ 2;
t176 = pkin(1) * t168 - pkin(6) * qJDD(1) + t204;
t35 = t49 * pkin(3) + qJDD(4) + t82;
t174 = t195 * t7 - t34 * t23 + t33 * t24 + t56 * t8 - t204;
t171 = g(3) * t149 - t113 * t85 + t204 * t150 - t257;
t5 = t23 * pkin(4) - t24 * qJ(5) - qJD(5) * t195 + t35;
t170 = -t113 * t87 + t175 + t265;
t154 = -qJ(4) - t256;
t146 = t155 * qJD(5);
t144 = pkin(2) * t227;
t137 = -t160 * pkin(3) - pkin(4);
t134 = t159 * pkin(3) + qJ(5);
t106 = qJ(5) * t233;
t105 = qJ(5) * t234;
t99 = pkin(1) + t229;
t92 = t166 * t99;
t75 = -pkin(4) - t79;
t74 = qJ(5) + t80;
t69 = t144 - t255;
t62 = -t159 * t95 + t160 * t96;
t61 = t159 * t96 + t160 * t95;
t51 = -t85 ^ 2 + t87 ^ 2;
t39 = -t159 * t65 - t160 * t64;
t38 = -t159 * t64 + t160 * t65;
t37 = -t196 + (-qJD(1) * t96 - t87) * t155;
t36 = t85 * t155 + t48;
t32 = t61 * pkin(4) - t62 * qJ(5) + t70;
t30 = t144 + t31;
t9 = t38 * pkin(4) - t39 * qJ(5) - t62 * qJD(5) + t59;
t2 = t245 + t259;
t1 = t143 + t146 + t4;
t3 = [qJDD(1), t260, t204, t156 * qJDD(1) + 0.2e1 * t162 * t211, 0.2e1 * t162 * t221 - 0.2e1 * t228 * t223, qJDD(2) * t162 + t167 * t165, qJDD(2) * t165 - t167 * t162, 0, t162 * t189 + t165 * t177, -t162 * t177 + t165 * t189, t48 * t96 + t87 * t64, -t48 * t95 - t96 * t49 + t64 * t85 + t87 * t65, t96 * t153 - t64 * t155, -t95 * t153 - t65 * t155, 0, -t113 * t65 + t145 * t85 + t150 * t260 + t153 * t192 + t155 * t172 - t246 * t49 + t82 * t95, t113 * t64 - t145 * t87 - t149 * t260 - t153 * t240 - t155 * t188 - t246 * t48 + t82 * t96, t70 * t23 + t35 * t61 + t67 * t38 - t56 * t59 + t178, t195 * t59 + t70 * t24 + t35 * t62 + t67 * t39 - t258, -t19 * t39 - t20 * t38 + t245 * t62 - t4 * t61 + t174, t4 * t34 + t20 * t8 + t245 * t33 - t19 * t7 + t35 * t70 + t67 * t59 - g(1) * (-t166 * t154 - t163 * t99) - g(2) * (-t163 * t154 + t92), t32 * t23 + t26 * t38 + t5 * t61 - t56 * t9 + t178, -t1 * t61 + t17 * t39 - t18 * t38 + t2 * t62 + t174, -t195 * t9 - t32 * t24 - t26 * t39 - t5 * t62 + t258, -g(2) * t92 + t1 * t34 + t17 * t7 + t18 * t8 + t2 * t33 + t26 * t9 + t5 * t32 + (g(1) * t154 - g(2) * t197) * t166 + (-g(1) * (-t197 - t99) + g(2) * t154) * t163; 0, 0, 0, -t162 * t168 * t165, t228 * t168, t222, t221, qJDD(2), -g(3) * t165 + t162 * t176, g(3) * t162 + t165 * t176, -t247, t51, t36, t37, t153, -t193 * t155 + (t164 * t153 - t155 * t225 - t227 * t85) * pkin(2) + t170, t241 * t155 + (-t161 * t153 - t155 * t224 + t227 * t87) * pkin(2) + t171, t79 * t153 + t56 * t69 + t180 - t213, -t80 * t153 - t155 * t242 - t195 * t69 + t185, -t80 * t23 - t79 * t24 + t242 * t56 + t201 + t216, -g(3) * t229 - t104 * t204 - t19 * t243 + t20 * t242 - t245 * t79 + t4 * t80 - t67 * t69, -t75 * t153 + t30 * t56 + t173 - t213, -t74 * t23 + t75 * t24 + t244 * t56 + t202 + t216, t74 * t153 + t155 * t244 + t195 * t30 + t146 + t181, t1 * t74 + t2 * t75 - t26 * t30 - g(1) * (t166 * t210 + t106) - g(2) * (t163 * t210 + t105) - g(3) * (t151 + t219) + t244 * t18 + t243 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t247, t51, t36, t37, t153, -t155 * t194 + t170, t155 * t206 + t171, t237 + (t153 * t160 - t56 * t87) * pkin(3) + t180, t236 + (-t153 * t159 + t195 * t87) * pkin(3) + t185, -t248 - t22 * t56 + (-t159 * t23 - t160 * t24) * pkin(3) + t201, t19 * t21 - t20 * t22 + (t159 * t4 - t160 * t245 + t67 * t87 + t265) * pkin(3), -t137 * t153 + t31 * t56 + t173 + t237, -t134 * t23 + t137 * t24 + t230 * t56 + t202 - t248, t134 * t153 + t195 * t31 + 0.2e1 * t146 + t181 - t236, t1 * t134 + t2 * t137 - t26 * t31 - t17 * t21 - g(1) * (t166 * t205 + t106) - g(2) * (t163 * t205 + t105) - g(3) * t219 + t230 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t191, t184, t200, t19 * t195 - t20 * t56 - t260 + t35, t191, t200, -t184, -t17 * t195 - t18 * t56 - t260 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t195 * t56 - t153, t24 - t264, -t155 ^ 2 - t261, -t18 * t155 - t173;];
tau_reg = t3;
