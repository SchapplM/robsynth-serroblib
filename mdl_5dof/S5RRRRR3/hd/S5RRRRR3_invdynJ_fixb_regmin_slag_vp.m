% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRRR3
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
% 
% Output:
% tau_reg [5x31]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:19
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:18:46
% EndTime: 2019-07-18 17:19:02
% DurationCPUTime: 4.91s
% Computational Cost: add. (3882->346), mult. (8692->491), div. (0->0), fcn. (7088->14), ass. (0->201)
t148 = qJ(2) + qJ(3);
t140 = sin(t148);
t153 = sin(qJ(1));
t157 = cos(qJ(1));
t182 = g(1) * t157 + g(2) * t153;
t174 = t182 * t140;
t142 = cos(t148);
t259 = g(3) * t142;
t288 = t259 - t174;
t150 = sin(qJ(4));
t154 = cos(qJ(5));
t231 = t150 * t154;
t149 = sin(qJ(5));
t155 = cos(qJ(4));
t233 = t149 * t155;
t102 = t231 + t233;
t147 = qJ(4) + qJ(5);
t139 = sin(t147);
t151 = sin(qJ(3));
t262 = pkin(1) * t151;
t138 = qJD(2) * t262;
t144 = qJD(2) + qJD(3);
t105 = pkin(5) * t144 + t138;
t156 = cos(qJ(2));
t224 = qJD(1) * t156;
t207 = pkin(1) * t224;
t263 = cos(qJ(3));
t197 = qJD(1) * t263;
t152 = sin(qJ(2));
t225 = qJD(1) * t152;
t285 = -t151 * t225 + t156 * t197;
t97 = -t151 * t224 - t152 * t197;
t57 = -pkin(2) * t285 + pkin(5) * t97 - t207;
t44 = t105 * t155 + t150 * t57;
t244 = t154 * t44;
t43 = -t105 * t150 + t155 * t57;
t92 = qJD(4) - t285;
t29 = pkin(3) * t92 + t43;
t17 = t149 * t29 + t244;
t143 = qJDD(2) + qJDD(3);
t192 = qJDD(1) * t263;
t212 = t156 * qJDD(1);
t286 = t144 * t285;
t49 = t151 * t212 + t152 * t192 + t286;
t73 = t144 * t150 - t155 * t97;
t28 = qJD(4) * t73 - t155 * t143 + t150 * t49;
t208 = t263 * pkin(1);
t227 = -qJD(3) * t138 + qJDD(2) * t208;
t84 = -pkin(2) * t143 - t227;
t23 = pkin(3) * t28 + t84;
t234 = t149 * t150;
t270 = t154 * t155;
t100 = t234 - t270;
t265 = qJD(4) + qJD(5);
t65 = t265 * t100;
t251 = -t100 * t285 + t65;
t196 = t263 * qJD(2);
t185 = pkin(1) * t196;
t106 = -t144 * pkin(2) - t185;
t71 = -t155 * t144 - t150 * t97;
t58 = t71 * pkin(3) + t106;
t287 = t23 * t102 + t288 * t139 - t17 * t97 - t251 * t58;
t221 = qJD(4) * t150;
t284 = (-t150 * t285 + t221) * pkin(3);
t283 = -qJD(5) * t150 - t221;
t278 = t265 * t102;
t281 = -t102 * t285 + t278;
t178 = t149 * t71 - t154 * t73;
t39 = t149 * t73 + t154 * t71;
t280 = t178 * t39;
t279 = t84 + t259;
t222 = qJD(3) * t151;
t206 = pkin(1) * t222;
t277 = t206 + t284;
t276 = t178 ^ 2 - t39 ^ 2;
t217 = qJD(5) * t154;
t219 = qJD(5) * t149;
t220 = qJD(4) * t155;
t27 = t150 * t143 + t144 * t220 + t155 * t49 + t221 * t97;
t8 = -t149 * t28 + t154 * t27 - t71 * t217 - t219 * t73;
t86 = qJD(5) + t92;
t275 = t39 * t86 + t8;
t131 = g(3) * t140;
t141 = cos(t147);
t36 = t44 * t219;
t236 = t142 * t153;
t77 = t139 * t157 - t141 * t236;
t235 = t142 * t157;
t79 = t139 * t153 + t141 * t235;
t274 = g(1) * t79 - g(2) * t77 + t141 * t131 + t58 * t39 + t36;
t215 = qJD(1) * qJD(2);
t194 = t152 * t215;
t213 = t152 * qJDD(1);
t179 = t151 * t213 - t156 * t192;
t103 = t151 * t156 + t263 * t152;
t68 = t144 * t103;
t240 = qJD(1) * t68;
t50 = t179 + t240;
t22 = t50 * pkin(2) - t49 * pkin(5) + (t194 - t212) * pkin(1);
t195 = t263 * qJD(3);
t214 = qJDD(2) * t151;
t85 = t143 * pkin(5) + (qJD(2) * t195 + t214) * pkin(1);
t10 = t43 * qJD(4) + t150 * t22 + t155 * t85;
t21 = t155 * t22;
t47 = qJDD(4) + t50;
t5 = t47 * pkin(3) - qJD(4) * t44 - t150 * t85 + t21;
t4 = t154 * t5;
t76 = t139 * t236 + t141 * t157;
t78 = -t139 * t235 + t141 * t153;
t273 = -g(1) * t78 + g(2) * t76 - qJD(5) * t17 - t149 * t10 + t139 * t131 + t58 * t178 + t4;
t163 = qJD(5) * t178 - t149 * t27 - t154 * t28;
t272 = -t178 * t86 + t163;
t134 = pkin(5) + t262;
t62 = -pkin(2) * t97 - pkin(5) * t285;
t56 = pkin(1) * t225 + t62;
t271 = (qJD(4) * t134 + t56) * t92;
t268 = t283 * t149;
t267 = -qJD(5) * t155 - t220;
t266 = qJD(1) * t103;
t172 = -t151 * t152 + t263 * t156;
t190 = qJD(5) * t29 + t10;
t223 = qJD(2) * t152;
t67 = t144 * t172;
t34 = pkin(1) * t223 + pkin(2) * t68 - pkin(5) * t67;
t46 = qJDD(5) + t47;
t63 = -pkin(1) * t156 - pkin(2) * t172 - pkin(5) * t103;
t48 = -pkin(3) * t172 + t155 * t63;
t264 = t172 * t190 - (qJD(5) * t48 + t150 * t34 + t220 * t63) * t86 - t150 * t63 * t46;
t258 = t155 * pkin(3);
t257 = t71 * t92;
t256 = t73 * t92;
t255 = t86 * t97;
t254 = t92 * t97;
t253 = t97 * t285;
t250 = pkin(1) * qJD(3);
t249 = t100 * t46;
t248 = t102 * t46;
t247 = t106 * t285;
t246 = t150 * t47;
t245 = t150 * t67;
t243 = t155 * t47;
t242 = t155 * t73;
t188 = t155 * t92;
t241 = t27 * t150;
t239 = qJD(4) * t92;
t238 = t103 * t106;
t232 = t150 * t153;
t230 = t150 * t157;
t229 = t153 * t155;
t228 = t155 * t157;
t145 = t152 ^ 2;
t226 = -t156 ^ 2 + t145;
t210 = pkin(5) * t239;
t204 = t150 * t263;
t203 = t155 * t263;
t202 = t263 * t144;
t191 = -qJD(4) * t57 - t85;
t186 = g(1) * t235 + g(2) * t236 + t207 * t285 + t131;
t135 = -t208 - pkin(2);
t184 = -t138 + t284;
t183 = -pkin(5) * t47 - t247;
t181 = g(1) * t153 - g(2) * t157;
t180 = t34 * t92 + t63 * t47;
t177 = t106 * t220 + t279 * t150 - t44 * t97;
t176 = t106 * t221 + t155 * t174 + t43 * t97;
t175 = t191 + t131;
t169 = -t97 * t207 + t227 - t288;
t167 = t84 * t103 + t106 * t67 - t63 * t239;
t165 = -t48 * t46 - (pkin(3) * t68 + t155 * t34 + t283 * t63) * t86;
t164 = -pkin(1) * t195 * t92 - t134 * t47 - t247;
t16 = -t149 * t44 + t154 * t29;
t162 = t23 * t100 - t141 * t288 + t16 * t97 + t281 * t58;
t159 = qJD(1) ^ 2;
t158 = qJD(2) ^ 2;
t136 = -pkin(2) - t258;
t112 = t135 - t258;
t91 = t97 * pkin(3);
t90 = t142 * t228 + t232;
t89 = -t142 * t230 + t229;
t88 = -t142 * t229 + t230;
t87 = t142 * t232 + t228;
t61 = t100 * t103;
t60 = t102 * t103;
t59 = t155 * t62;
t52 = t150 * t62 + t155 * t185;
t51 = -t285 ^ 2 + t97 ^ 2;
t37 = t155 * t56 - t91;
t35 = -t150 * t185 + t59 - t91;
t31 = -t179 + (-t266 - t97) * t144;
t30 = t49 - t286;
t15 = t245 * t154 + t67 * t233 + (t265 * t270 + t268) * t103;
t14 = -t100 * t67 - t103 * t278;
t13 = t188 * t92 + t73 * t97 + t246;
t12 = -t92 ^ 2 * t150 - t71 * t97 + t243;
t11 = t188 * t73 + t241;
t7 = -t281 * t86 - t39 * t97 - t249;
t6 = -t178 * t97 - t251 * t86 + t248;
t3 = (t27 - t257) * t155 + (-t28 - t256) * t150;
t2 = t8 * t102 + t178 * t251;
t1 = -t8 * t100 + t102 * t163 + t178 * t281 + t251 * t39;
t9 = [qJDD(1), t181, t182, qJDD(1) * t145 + 0.2e1 * t156 * t194, 0.2e1 * t152 * t212 - 0.2e1 * t215 * t226, qJDD(2) * t152 + t156 * t158, qJDD(2) * t156 - t152 * t158, 0, t181 * t156, -t181 * t152, t103 * t49 - t67 * t97, -t103 * t50 + t172 * t49 + t285 * t67 + t68 * t97, t103 * t143 + t144 * t67, t143 * t172 - t144 * t68, 0, t181 * t142 + ((-qJD(1) * t172 - t285) * t223 + (qJDD(1) * t172 - t240 - t50) * t156) * pkin(1), -t181 * t140 + ((-t97 + t266) * t223 + (-qJD(1) * t67 - qJDD(1) * t103 - t49) * t156) * pkin(1), t67 * t242 + (t27 * t155 - t221 * t73) * t103, (-t150 * t73 - t155 * t71) * t67 + (-t241 - t155 * t28 + (t150 * t71 - t242) * qJD(4)) * t103, t67 * t188 - t27 * t172 + t73 * t68 + (-t221 * t92 + t243) * t103, -t92 * t245 + t28 * t172 - t71 * t68 + (-t220 * t92 - t246) * t103, -t172 * t47 + t68 * t92, -g(1) * t88 - g(2) * t90 - t21 * t172 + t43 * t68 + ((t105 * t172 + t238) * qJD(4) + t180) * t155 + (-t172 * t191 + t167) * t150, -g(1) * t87 - g(2) * t89 + t10 * t172 - t44 * t68 + t167 * t155 + (-qJD(4) * t238 - t180) * t150, -t14 * t178 - t61 * t8, -t14 * t39 + t15 * t178 - t163 * t61 - t60 * t8, t14 * t86 - t172 * t8 - t178 * t68 - t46 * t61, -t15 * t86 - t163 * t172 - t39 * t68 - t46 * t60, -t172 * t46 + t68 * t86, -g(1) * t77 - g(2) * t79 - t4 * t172 + t58 * t15 + t16 * t68 + t23 * t60 + (qJD(5) * t172 * t44 - t165) * t154 + t264 * t149 + (t39 * t245 + (-t150 * t163 + t220 * t39) * t103) * pkin(3), -g(1) * t76 - g(2) * t78 - t36 * t172 + t58 * t14 - t17 * t68 - t23 * t61 + (t172 * t5 + t165) * t149 + t264 * t154 + (-t178 * t245 + (t150 * t8 - t178 * t220) * t103) * pkin(3); 0, 0, 0, -t152 * t159 * t156, t226 * t159, t213, t212, qJDD(2), -g(3) * t156 + t152 * t182, g(3) * t152 + t156 * t182, t253, t51, t30, t31, t143, (t263 * t143 - t144 * t222 + t225 * t285) * pkin(1) + t169, (t97 * t225 + (-qJDD(2) - t143) * t151 + (-t196 - t202) * qJD(3)) * pkin(1) + t186, t11, t3, t13, t12, t254, t71 * t206 + t135 * t28 + t164 * t150 + (-t279 - t271) * t155 + t176, t73 * t206 + t135 * t27 + t164 * t155 + (-t174 + t271) * t150 + t177, t2, t1, t6, t7, t255, -t112 * t163 - t134 * t248 + t162 + ((-t149 * t203 - t154 * t204) * t250 + t134 * t65 - t154 * t37 + t234 * t56) * t86 + t277 * t39, t112 * t8 + t134 * t249 + (-(-t149 * t204 + t154 * t203) * t250 + t278 * t134 + t149 * t37 + t231 * t56) * t86 - t277 * t178 + t287; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t253, t51, t30, t31, t143, t138 * t144 + t169, (-t214 + (-t195 + t202) * qJD(2)) * pkin(1) + t186, t11, t3, t13, t12, t254, -t71 * t138 - pkin(2) * t28 - t59 * t92 + (t185 * t92 + t183) * t150 + (-t279 - t210) * t155 + t176, -t73 * t138 - pkin(2) * t27 + t52 * t92 + t183 * t155 + (-t174 + t210) * t150 + t177, t2, t1, t6, t7, t255, -t136 * t163 - (-t149 * t52 + t154 * t35) * t86 + t184 * t39 + ((t267 * t154 - t268) * t86 - t248) * pkin(5) + t162, t136 * t8 + (t149 * t35 + t154 * t52) * t86 - t184 * t178 + (-(t267 * t149 - t150 * t217 - t154 * t221) * t86 + t249) * pkin(5) + t287; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73 * t71, -t71 ^ 2 + t73 ^ 2, t27 + t257, t256 - t28, t47, -g(1) * t89 + g(2) * t87 - t105 * t220 - t106 * t73 + t150 * t175 + t44 * t92 + t21, g(1) * t90 - g(2) * t88 + t106 * t71 + t43 * t92 + (qJD(4) * t105 - t22) * t150 + t175 * t155, -t280, t276, t275, t272, t46, -(-t149 * t43 - t244) * t86 + (t154 * t46 - t219 * t86 - t39 * t73) * pkin(3) + t273, (-t44 * t86 - t5) * t149 + (t43 * t86 - t190) * t154 + (-t149 * t46 + t178 * t73 - t217 * t86) * pkin(3) + t274; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t280, t276, t275, t272, t46, t17 * t86 + t273, -t149 * t5 - t154 * t190 + t16 * t86 + t274;];
tau_reg  = t9;
