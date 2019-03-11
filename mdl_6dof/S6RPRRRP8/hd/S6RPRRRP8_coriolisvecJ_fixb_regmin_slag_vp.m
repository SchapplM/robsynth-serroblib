% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% tauc_reg [6x31]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRRP8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:24:48
% EndTime: 2019-03-09 06:24:53
% DurationCPUTime: 2.67s
% Computational Cost: add. (5605->353), mult. (11589->461), div. (0->0), fcn. (7654->6), ass. (0->191)
t147 = sin(qJ(3));
t149 = cos(qJ(3));
t252 = sin(qJ(4));
t253 = cos(qJ(4));
t112 = t253 * t147 + t252 * t149;
t105 = t112 * qJD(1);
t268 = qJD(5) + t105;
t273 = t268 ^ 2;
t198 = t252 * t147;
t184 = qJD(1) * t198;
t199 = t253 * t149;
t104 = -qJD(1) * t199 + t184;
t146 = sin(qJ(5));
t148 = cos(qJ(5));
t205 = qJD(3) + qJD(4);
t189 = t148 * t205;
t82 = -t146 * t104 - t189;
t272 = t268 * t82;
t166 = t148 * t104 - t146 * t205;
t271 = t166 * t268;
t208 = qJD(5) * t148;
t221 = t105 * t148;
t270 = t208 + t221;
t209 = qJD(5) * t146;
t222 = t105 * t146;
t269 = t209 + t222;
t111 = t198 - t199;
t188 = qJD(5) * t112 + qJD(1);
t194 = t253 * qJD(4);
t164 = t253 * qJD(3) + t194;
t162 = t164 * t149;
t216 = t205 * t184;
t73 = qJD(1) * t162 - t216;
t234 = t148 * t73;
t260 = t105 * t205;
t43 = -qJD(5) * t189 - t104 * t209 + t148 * t260;
t196 = qJD(4) * t252;
t258 = -qJD(3) * t252 - t196;
t78 = -t164 * t147 + t258 * t149;
t79 = t258 * t147 + t162;
t267 = (t146 * t188 - t148 * t79) * t268 - t111 * t43 - t112 * t234 + t166 * t78;
t150 = -pkin(1) - pkin(7);
t122 = t150 * qJD(1) + qJD(2);
t212 = qJD(1) * t149;
t100 = -pkin(8) * t212 + t149 * t122;
t95 = qJD(3) * pkin(3) + t100;
t191 = pkin(8) * qJD(1) - t122;
t211 = qJD(3) * t147;
t96 = t191 * t211;
t210 = qJD(3) * t149;
t97 = t191 * t210;
t213 = qJD(1) * t147;
t99 = -pkin(8) * t213 + t122 * t147;
t30 = t95 * t194 - t99 * t196 + t252 * t96 - t253 * t97;
t142 = qJD(1) * qJD(2);
t206 = qJD(1) * qJD(3);
t193 = t149 * t206;
t114 = pkin(3) * t193 + t142;
t34 = t73 * pkin(4) + pkin(9) * t260 + t114;
t94 = t253 * t99;
t63 = t252 * t95 + t94;
t56 = t205 * pkin(9) + t63;
t118 = pkin(3) * t213 + qJD(1) * qJ(2);
t64 = pkin(4) * t105 + pkin(9) * t104 + t118;
t168 = t146 * t34 + t148 * t30 + t64 * t208 - t56 * t209;
t226 = t73 * qJ(6);
t2 = qJD(6) * t268 + t168 + t226;
t192 = t146 * t30 - t148 * t34 + t56 * t208 + t64 * t209;
t254 = pkin(5) * t73;
t4 = t192 - t254;
t266 = t4 * t146 + t2 * t148;
t265 = (t209 * t268 - t234) * pkin(9);
t264 = t79 * t205;
t26 = t146 * t64 + t148 * t56;
t17 = qJ(6) * t268 + t26;
t93 = t252 * t99;
t62 = t253 * t95 - t93;
t55 = -t205 * pkin(4) - t62;
t22 = t82 * pkin(5) + qJ(6) * t166 + t55;
t263 = t17 * t104 - t22 * t221;
t66 = t252 * t100 + t94;
t175 = pkin(3) * t196 - t66;
t31 = t99 * t194 + t95 * t196 - t252 * t97 - t253 * t96;
t44 = -t166 * qJD(5) - t146 * t260;
t6 = pkin(5) * t44 + qJ(6) * t43 + qJD(6) * t166 + t31;
t262 = t6 * t148 - t22 * t209;
t261 = -t31 * t148 + t55 * t209;
t259 = pkin(5) * t269 - qJ(6) * t270 - t146 * qJD(6);
t246 = pkin(8) - t150;
t115 = t246 * t147;
t116 = t246 * t149;
t80 = -t252 * t115 + t253 * t116;
t256 = t166 ^ 2;
t255 = 0.2e1 * t142;
t251 = pkin(5) * t104;
t250 = t22 * t166;
t249 = t6 * t146;
t247 = t166 * t82;
t245 = -t259 - t175;
t74 = -pkin(4) * t104 + pkin(9) * t105;
t244 = t146 * t74 + t148 * t62;
t67 = t253 * t100 - t93;
t68 = pkin(3) * t212 + t74;
t243 = t146 * t68 + t148 * t67;
t133 = t147 * pkin(3) + qJ(2);
t75 = pkin(4) * t112 + pkin(9) * t111 + t133;
t81 = -t253 * t115 - t252 * t116;
t242 = t146 * t75 + t148 * t81;
t241 = pkin(3) * qJD(4);
t135 = t252 * pkin(3) + pkin(9);
t240 = t135 * t73;
t239 = t146 * t73;
t238 = t146 * t78;
t237 = t146 * t82;
t236 = t146 * t166;
t235 = t148 * t44;
t233 = t148 * t78;
t231 = t148 * t82;
t230 = t148 * t166;
t227 = t43 * t146;
t225 = -t63 + t259;
t224 = t268 * t104;
t223 = t104 * t105;
t151 = qJD(3) ^ 2;
t220 = t151 * t147;
t219 = t151 * t149;
t152 = qJD(1) ^ 2;
t218 = t152 * qJ(2);
t25 = -t146 * t56 + t148 * t64;
t217 = qJD(6) - t25;
t215 = t147 ^ 2 - t149 ^ 2;
t214 = -t151 - t152;
t123 = pkin(3) * t210 + qJD(2);
t204 = 0.2e1 * qJD(1);
t203 = t253 * pkin(3);
t202 = t146 * t253;
t201 = t148 * t253;
t190 = t148 * t268;
t187 = pkin(3) * t194;
t185 = -t26 * t104 + t31 * t146 + t55 * t208;
t182 = -t148 * pkin(5) - t146 * qJ(6);
t181 = pkin(5) * t146 - qJ(6) * t148;
t180 = t105 * t55 - t240;
t16 = -pkin(5) * t268 + t217;
t179 = -t146 * t17 + t148 * t16;
t178 = t146 * t16 + t148 * t17;
t177 = -t146 * t62 + t148 * t74;
t176 = -t230 + t237;
t117 = -pkin(4) + t182;
t174 = -t16 * t104 - t262;
t173 = t25 * t104 + t261;
t172 = -t22 * t208 - t249;
t171 = t26 * t268 - t192;
t170 = t118 * t104 - t31;
t41 = pkin(4) * t79 - pkin(9) * t78 + t123;
t109 = t246 * t211;
t110 = qJD(3) * t116;
t45 = -t80 * qJD(4) + t252 * t109 - t253 * t110;
t167 = t146 * t41 + t148 * t45 + t75 * t208 - t81 * t209;
t165 = (-t208 * t268 - t239) * pkin(9);
t161 = -t135 * t209 + t148 * t187;
t160 = t16 * t270 - t17 * t269 + t266;
t159 = t179 * qJD(5) + t266;
t158 = t176 * qJD(5) - t227 - t235;
t157 = -t164 * t212 + t216;
t156 = t118 * t105 - t30;
t155 = -t112 * t239 + t111 * t44 - t78 * t82 + (-t146 * t79 - t188 * t148) * t268;
t46 = t81 * qJD(4) - t253 * t109 - t252 * t110;
t136 = -t203 - pkin(4);
t107 = -t203 + t117;
t101 = t104 * qJ(6);
t76 = t78 * t205;
t58 = t104 ^ 2 - t105 ^ 2;
t54 = -t104 * t205 + t157;
t49 = -pkin(5) * t166 + qJ(6) * t82;
t47 = -t181 * t111 + t80;
t37 = -pkin(5) * t112 + t146 * t81 - t148 * t75;
t36 = qJ(6) * t112 + t242;
t24 = -t177 + t251;
t23 = -t101 + t244;
t21 = t146 * t67 - t148 * t68 + t251;
t20 = -t101 + t243;
t18 = -t43 + t272;
t12 = -t104 * t166 + t190 * t268 + t239;
t11 = -t82 * t104 - t146 * t273 + t234;
t10 = -t166 * t190 - t227;
t9 = t181 * t78 + (t182 * qJD(5) + qJD(6) * t148) * t111 + t46;
t8 = -t79 * pkin(5) + t242 * qJD(5) + t146 * t45 - t148 * t41;
t7 = qJ(6) * t79 + qJD(6) * t112 + t167;
t5 = (-t43 - t272) * t148 + (-t44 + t271) * t146;
t1 = [0, 0, 0, 0, t255, qJ(2) * t255, -0.2e1 * t147 * t193, 0.2e1 * t215 * t206, -t220, -t219, 0, -t150 * t220 + (qJ(2) * t210 + qJD(2) * t147) * t204, -t150 * t219 + (-qJ(2) * t211 + qJD(2) * t149) * t204, -t104 * t78 + t111 * t260, t104 * t79 - t78 * t105 + t111 * t73 + t112 * t260, t76, -t264, 0, t123 * t105 + t114 * t112 + t118 * t79 + t133 * t73 - t46 * t205, -t123 * t104 - t114 * t111 + t118 * t78 - t133 * t260 - t45 * t205, -t78 * t230 + (t43 * t148 - t166 * t209) * t111 (-t231 + t236) * t78 + (-t227 + t235 + (-t230 - t237) * qJD(5)) * t111, -t111 * t234 - t43 * t112 - t166 * t79 + (t111 * t209 + t233) * t268, t111 * t239 - t44 * t112 - t82 * t79 + (t111 * t208 - t238) * t268, t112 * t73 + t268 * t79, -t192 * t112 + t25 * t79 + t46 * t82 + t80 * t44 + ((-qJD(5) * t81 + t41) * t268 + t75 * t73 - t55 * qJD(5) * t111) * t148 + ((-qJD(5) * t75 - t45) * t268 - t81 * t73 - t31 * t111 + t55 * t78) * t146, t261 * t111 - t168 * t112 - t166 * t46 - t167 * t268 + t55 * t233 - t242 * t73 - t26 * t79 - t80 * t43, t111 * t172 - t4 * t112 - t16 * t79 + t22 * t238 - t268 * t8 - t37 * t73 + t47 * t44 + t9 * t82, -t36 * t44 - t37 * t43 - t7 * t82 - t8 * t166 + t179 * t78 + (qJD(5) * t178 + t146 * t2 - t148 * t4) * t111, t262 * t111 + t2 * t112 + t166 * t9 + t17 * t79 - t22 * t233 + t268 * t7 + t36 * t73 + t47 * t43, t16 * t8 + t17 * t7 + t2 * t36 + t22 * t9 + t37 * t4 + t47 * t6; 0, 0, 0, 0, -t152, -t218, 0, 0, 0, 0, 0, t214 * t147, t214 * t149, 0, 0, 0, 0, 0, -qJD(1) * t105 + t76, qJD(1) * t104 - t264, 0, 0, 0, 0, 0, t155, t267, t155 (-t231 - t236) * t79 + t176 * qJD(1) + t158 * t112, -t267, qJD(1) * t179 + t6 * t111 + t112 * t159 + t178 * t79 - t22 * t78; 0, 0, 0, 0, 0, 0, t149 * t152 * t147, -t215 * t152, 0, 0, 0, -t149 * t218, t147 * t218, -t223, t58, 0, t54, 0, t66 * t205 + (-t105 * t212 - t205 * t196) * pkin(3) + t170, t67 * t205 + (t104 * t212 - t205 * t194) * pkin(3) + t156, t10, t5, t12, t11, t224, t136 * t44 + t175 * t82 + t180 * t146 + ((-qJD(5) * t135 - t68) * t148 + (-t187 + t67) * t146) * t268 + t173, -t136 * t43 - t175 * t166 + t180 * t148 + (-t161 + t243) * t268 + t185, t107 * t44 - t245 * t82 + (t105 * t22 - t240) * t146 + (-t135 * t208 - t146 * t187 + t21) * t268 + t174, t20 * t82 + t21 * t166 + (-t166 * t202 - t201 * t82) * t241 + t158 * t135 + t160, t107 * t43 - t249 - t245 * t166 + (-qJD(5) * t22 + t240) * t148 + (-t20 + t161) * t268 + t263, t6 * t107 - t16 * t21 - t17 * t20 - t245 * t22 + (t16 * t202 + t17 * t201) * t241 + t159 * t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t223, t58, 0, t54, 0, t63 * t205 + t170, t62 * t205 + t156, t10, t5, t12, t11, t224, -pkin(4) * t44 - t177 * t268 + t222 * t55 - t63 * t82 + t165 + t173, pkin(4) * t43 + t166 * t63 + t55 * t221 + t244 * t268 + t185 + t265, t117 * t44 + t22 * t222 + t225 * t82 + t24 * t268 + t165 + t174, pkin(9) * t158 + t166 * t24 + t23 * t82 + t160, t117 * t43 + t166 * t225 - t23 * t268 + t172 + t263 - t265, pkin(9) * t159 + t6 * t117 - t16 * t24 - t17 * t23 + t22 * t225; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t247, -t82 ^ 2 + t256, t18, -t44 - t271, t73, t166 * t55 + t171, t25 * t268 + t55 * t82 - t168, -t49 * t82 + t171 + t250 + 0.2e1 * t254, pkin(5) * t43 - t44 * qJ(6) - (t17 - t26) * t166 + (t16 - t217) * t82, 0.2e1 * t226 - t22 * t82 - t49 * t166 + (0.2e1 * qJD(6) - t25) * t268 + t168, -t4 * pkin(5) + t2 * qJ(6) - t16 * t26 + t17 * t217 - t22 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t157 - t247, t18, -t256 - t273, -t17 * t268 - t250 + t4;];
tauc_reg  = t1;
