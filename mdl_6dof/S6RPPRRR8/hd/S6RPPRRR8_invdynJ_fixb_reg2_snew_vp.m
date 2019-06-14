% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPPRRR8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 16:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPPRRR8_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR8_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR8_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR8_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:14:55
% EndTime: 2019-05-05 16:15:11
% DurationCPUTime: 6.66s
% Computational Cost: add. (28130->432), mult. (64445->611), div. (0->0), fcn. (47276->10), ass. (0->273)
t233 = sin(qJ(4));
t228 = sin(pkin(10));
t237 = cos(qJ(4));
t229 = cos(pkin(10));
t287 = t229 * t233;
t255 = t228 * t237 + t287;
t210 = t255 * qJD(1);
t288 = t228 * t233;
t212 = (t229 * t237 - t288) * qJD(1);
t289 = t212 * t210;
t320 = qJDD(4) - t289;
t322 = t233 * t320;
t321 = t237 * t320;
t231 = sin(qJ(6));
t232 = sin(qJ(5));
t236 = cos(qJ(5));
t193 = -t236 * qJD(4) + t232 * t212;
t195 = t232 * qJD(4) + t236 * t212;
t235 = cos(qJ(6));
t162 = t235 * t193 + t231 * t195;
t164 = -t231 * t193 + t235 * t195;
t128 = t164 * t162;
t273 = t212 * qJD(4);
t308 = t255 * qJDD(1);
t186 = -t308 - t273;
t177 = qJDD(5) - t186;
t176 = qJDD(6) + t177;
t311 = -t128 + t176;
t319 = t231 * t311;
t169 = t195 * t193;
t309 = -t169 + t177;
t318 = t232 * t309;
t317 = t235 * t311;
t316 = t236 * t309;
t239 = qJD(1) ^ 2;
t234 = sin(qJ(1));
t238 = cos(qJ(1));
t256 = t234 * g(1) - t238 * g(2);
t251 = qJDD(2) - t256;
t246 = -t239 * qJ(2) + t251;
t266 = -0.2e1 * qJD(3) * qJD(1);
t300 = pkin(1) + qJ(3);
t315 = -t300 * qJDD(1) + t246 + t266;
t268 = t228 * qJDD(1);
t224 = t228 ^ 2;
t225 = t229 ^ 2;
t275 = t224 + t225;
t314 = pkin(3) * t268 - (t275 * pkin(7) + t300) * t239;
t313 = -pkin(7) - t300;
t312 = t275 * t239;
t267 = t229 * qJDD(1);
t209 = -t233 * t268 + t237 * t267;
t274 = t210 * qJD(4);
t188 = t209 - t274;
t254 = -t232 * qJDD(4) - t236 * t188;
t155 = -t193 * qJD(5) - t254;
t259 = -t236 * qJDD(4) + t232 * t188;
t250 = t195 * qJD(5) + t259;
t102 = -t162 * qJD(6) + t235 * t155 - t231 * t250;
t205 = qJD(5) + t210;
t202 = qJD(6) + t205;
t151 = t202 * t162;
t310 = -t151 + t102;
t174 = t205 * t193;
t133 = t155 + t174;
t260 = t231 * t155 + t235 * t250;
t81 = (qJD(6) - t202) * t164 + t260;
t129 = (qJD(5) - t205) * t195 + t259;
t160 = t162 ^ 2;
t161 = t164 ^ 2;
t305 = t193 ^ 2;
t192 = t195 ^ 2;
t201 = t202 ^ 2;
t204 = t205 ^ 2;
t207 = t210 ^ 2;
t208 = t212 ^ 2;
t304 = qJD(4) ^ 2;
t180 = t210 * pkin(4) - t212 * pkin(8);
t243 = (t266 + (-pkin(3) * t228 - qJ(2)) * t239 + t313 * qJDD(1) + t251) * t229;
t301 = t228 * g(3);
t242 = t243 + t301;
t179 = -t229 * g(3) + t228 * t315;
t170 = -t224 * t239 * pkin(3) - pkin(7) * t268 + t179;
t278 = t237 * t170;
t106 = -t304 * pkin(4) + qJDD(4) * pkin(8) - t210 * t180 + t233 * t242 + t278;
t269 = qJD(2) * qJD(1);
t223 = 0.2e1 * t269;
t226 = qJDD(1) * qJ(2);
t257 = t238 * g(1) + t234 * g(2);
t252 = -t226 + t257;
t249 = -qJDD(3) + t252;
t117 = t223 + (-t188 + t274) * pkin(8) + (-t186 + t273) * pkin(4) - t249 + t314;
t70 = t232 * t106 - t236 * t117;
t56 = pkin(5) * t309 - pkin(9) * t133 - t70;
t171 = t205 * pkin(5) - t195 * pkin(9);
t71 = t236 * t106 + t232 * t117;
t60 = -t305 * pkin(5) - t250 * pkin(9) - t205 * t171 + t71;
t29 = t231 * t60 - t235 * t56;
t30 = t231 * t56 + t235 * t60;
t14 = t231 * t30 - t235 * t29;
t303 = pkin(5) * t14;
t84 = t151 + t102;
t51 = -t231 * t81 - t235 * t84;
t302 = pkin(5) * t51;
t136 = t233 * t170 - t237 * t242;
t137 = g(3) * t288 + t233 * t243 + t278;
t99 = -t237 * t136 + t233 * t137;
t299 = t229 * t99;
t105 = -qJDD(4) * pkin(4) - t304 * pkin(8) + t212 * t180 + t136;
t65 = pkin(5) * t250 - t305 * pkin(9) + t195 * t171 + t105;
t298 = t231 * t65;
t297 = t232 * t14;
t296 = t235 * t65;
t295 = t236 * t14;
t294 = qJDD(1) * pkin(1);
t293 = t202 * t231;
t292 = t202 * t235;
t291 = t205 * t232;
t290 = t205 * t236;
t119 = t128 + t176;
t286 = t231 * t119;
t285 = t232 * t105;
t141 = t169 + t177;
t284 = t232 * t141;
t245 = t249 - 0.2e1 * t269;
t175 = t245 - t314;
t283 = t233 * t175;
t183 = qJDD(4) + t289;
t282 = t233 * t183;
t281 = t235 * t119;
t280 = t236 * t105;
t279 = t236 * t141;
t277 = t237 * t175;
t276 = t237 * t183;
t271 = qJD(5) + t205;
t265 = t233 * t128;
t264 = t233 * t169;
t263 = t237 * t128;
t262 = t237 * t169;
t261 = -pkin(4) * t237 - pkin(3);
t15 = t231 * t29 + t235 * t30;
t41 = t232 * t70 + t236 * t71;
t100 = t233 * t136 + t237 * t137;
t197 = t300 * t239 + t245;
t258 = -t197 + t226;
t21 = t228 * (t233 * t105 + t237 * t41) + t229 * (-t237 * t105 + t233 * t41);
t40 = t232 * t71 - t236 * t70;
t143 = t229 * (t229 * t315 + t301) + t228 * t179;
t123 = -t201 - t160;
t74 = t231 * t123 + t317;
t253 = pkin(5) * t74 - t29;
t138 = -t161 - t201;
t90 = t235 * t138 - t286;
t247 = pkin(5) * t90 - t30;
t215 = t275 * qJDD(1);
t214 = t228 * t312;
t213 = t229 * t312;
t206 = -t246 + t294;
t200 = -t208 - t304;
t199 = -t208 + t304;
t198 = t207 - t304;
t187 = t209 - 0.2e1 * t274;
t185 = t308 + 0.2e1 * t273;
t181 = -t304 - t207;
t173 = -t192 + t204;
t172 = -t204 + t305;
t168 = -t207 - t208;
t167 = t192 - t305;
t159 = -t192 - t204;
t158 = -t233 * t200 - t276;
t157 = t237 * t200 - t282;
t156 = -t204 - t305;
t150 = t192 + t305;
t149 = t233 * t209 - t237 * t308;
t148 = -t237 * t209 - t233 * t308;
t147 = -t161 + t201;
t146 = t160 - t201;
t145 = t237 * t181 - t322;
t144 = t233 * t181 + t321;
t139 = (-t193 * t236 + t195 * t232) * t205;
t134 = t193 * t271 + t254;
t132 = t155 - t174;
t130 = -t195 * t271 - t259;
t127 = t161 - t160;
t126 = t236 * t155 - t195 * t291;
t125 = t193 * t290 + t232 * t250;
t124 = t229 * t157 + t228 * t158;
t122 = t236 * t172 - t284;
t121 = -t232 * t173 + t316;
t114 = -t232 * t159 - t279;
t113 = t236 * t159 - t284;
t112 = (-t162 * t235 + t164 * t231) * t202;
t111 = (-t162 * t231 - t164 * t235) * t202;
t110 = t229 * t148 + t228 * t149;
t109 = t236 * t156 - t318;
t108 = t232 * t156 + t316;
t107 = t229 * t144 + t228 * t145;
t103 = -t160 - t161;
t101 = -t164 * qJD(6) - t260;
t98 = -t129 * t236 + t232 * t133;
t97 = t236 * t130 - t232 * t132;
t96 = -t129 * t232 - t236 * t133;
t95 = t235 * t146 - t286;
t94 = -t231 * t147 + t317;
t93 = t231 * t146 + t281;
t92 = t235 * t147 + t319;
t91 = -t231 * t138 - t281;
t89 = t237 * t114 - t233 * t134;
t88 = t233 * t114 + t237 * t134;
t87 = t237 * t109 - t233 * t130;
t86 = t233 * t109 + t237 * t130;
t80 = (qJD(6) + t202) * t164 + t260;
t79 = t235 * t102 - t164 * t293;
t78 = t231 * t102 + t164 * t292;
t77 = -t231 * t101 + t162 * t292;
t76 = t235 * t101 + t162 * t293;
t75 = t235 * t123 - t319;
t73 = -t233 * t150 + t237 * t98;
t72 = t237 * t150 + t233 * t98;
t69 = -t232 * t111 + t236 * t112;
t68 = -pkin(8) * t113 + t280;
t66 = -pkin(8) * t108 + t285;
t64 = t228 * t100 + t299;
t63 = -pkin(4) * t113 + t71;
t62 = -t232 * t93 + t236 * t95;
t61 = -t232 * t92 + t236 * t94;
t59 = -pkin(4) * t108 + t70;
t58 = -t232 * t90 + t236 * t91;
t57 = t232 * t91 + t236 * t90;
t54 = t228 * t89 + t229 * t88;
t53 = t231 * t84 - t235 * t81;
t52 = -t231 * t310 - t235 * t80;
t50 = -t231 * t80 + t235 * t310;
t49 = t228 * t87 + t229 * t86;
t48 = -t232 * t78 + t236 * t79;
t47 = -t232 * t76 + t236 * t77;
t46 = -t232 * t74 + t236 * t75;
t45 = t232 * t75 + t236 * t74;
t44 = -pkin(9) * t90 + t296;
t43 = -pkin(9) * t74 + t298;
t42 = t228 * t73 + t229 * t72;
t39 = t233 * t310 + t237 * t58;
t38 = t233 * t58 - t237 * t310;
t35 = -pkin(5) * t310 + pkin(9) * t91 + t298;
t34 = t233 * t80 + t237 * t46;
t33 = t233 * t46 - t237 * t80;
t32 = -pkin(8) * t96 - t40;
t31 = -pkin(5) * t80 + pkin(9) * t75 - t296;
t27 = -t232 * t51 + t236 * t53;
t26 = -t232 * t50 + t236 * t52;
t25 = t232 * t53 + t236 * t51;
t24 = t233 * t103 + t237 * t27;
t23 = -t237 * t103 + t233 * t27;
t22 = t228 * t39 + t229 * t38;
t20 = -pkin(4) * t25 - t302;
t19 = t228 * t34 + t229 * t33;
t18 = -pkin(4) * t57 - t247;
t17 = -pkin(4) * t45 - t253;
t16 = -pkin(8) * t57 - t232 * t35 + t236 * t44;
t13 = -pkin(8) * t45 - t232 * t31 + t236 * t43;
t12 = -pkin(5) * t65 + pkin(9) * t15;
t11 = -pkin(9) * t51 - t14;
t10 = t228 * t24 + t229 * t23;
t9 = -pkin(5) * t103 + pkin(9) * t53 + t15;
t8 = t236 * t15 - t297;
t7 = t232 * t15 + t295;
t6 = t233 * t65 + t237 * t8;
t5 = t233 * t8 - t237 * t65;
t4 = -pkin(4) * t7 - t303;
t3 = -pkin(8) * t25 + t236 * t11 - t232 * t9;
t2 = -pkin(8) * t7 - pkin(9) * t295 - t232 * t12;
t1 = t228 * t6 + t229 * t5;
t28 = [0, 0, 0, 0, 0, qJDD(1), t256, t257, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t251 - 0.2e1 * t294, t223 + 0.2e1 * t226 - t257, pkin(1) * t206 + qJ(2) * (-t239 * pkin(1) + t223 - t252), t225 * qJDD(1), -0.2e1 * t228 * t267, 0, t224 * qJDD(1), 0, 0, t300 * t214 + t258 * t228, t300 * t213 + t258 * t229, -qJ(2) * t312 + t300 * t215 - t143, -qJ(2) * t197 - t300 * t143, t229 * (t237 * t188 - t233 * t273) - t228 * (t233 * t188 + t237 * t273), t229 * (-t237 * t185 - t233 * t187) - t228 * (-t233 * t185 + t237 * t187), t229 * (-t233 * t199 + t321) - t228 * (t237 * t199 + t322), t229 * (-t233 * t186 + t237 * t274) - t228 * (t237 * t186 + t233 * t274), t229 * (t237 * t198 - t282) - t228 * (t233 * t198 + t276), (t229 * (-t210 * t237 + t212 * t233) - t228 * (-t210 * t233 - t212 * t237)) * qJD(4), t229 * (-pkin(7) * t144 - t283) - t228 * (-pkin(3) * t185 + pkin(7) * t145 + t277) + qJ(2) * t185 - t300 * t107, t229 * (-pkin(7) * t157 - t277) - t228 * (-pkin(3) * t187 + pkin(7) * t158 - t283) + qJ(2) * t187 - t300 * t124, t229 * (-pkin(7) * t148 - t99) - t228 * (-pkin(3) * t168 + pkin(7) * t149 + t100) + qJ(2) * t168 - t300 * t110, -pkin(7) * t299 - t228 * (pkin(3) * t175 + pkin(7) * t100) - qJ(2) * t175 - t300 * t64, t229 * (t237 * t126 + t264) - t228 * (t233 * t126 - t262), t229 * (t233 * t167 + t237 * t97) - t228 * (-t237 * t167 + t233 * t97), t229 * (t237 * t121 + t233 * t133) - t228 * (t233 * t121 - t237 * t133), t229 * (t237 * t125 - t264) - t228 * (t233 * t125 + t262), t229 * (t237 * t122 - t233 * t129) - t228 * (t233 * t122 + t237 * t129), t229 * (t237 * t139 + t233 * t177) - t228 * (t233 * t139 - t237 * t177), t229 * (-pkin(7) * t86 - t233 * t59 + t237 * t66) - t228 * (-pkin(3) * t108 + pkin(7) * t87 + t233 * t66 + t237 * t59) + qJ(2) * t108 - t300 * t49, t229 * (-pkin(7) * t88 - t233 * t63 + t237 * t68) - t228 * (-pkin(3) * t113 + pkin(7) * t89 + t233 * t68 + t237 * t63) + qJ(2) * t113 - t300 * t54, t229 * (-pkin(7) * t72 + t237 * t32) - t228 * (pkin(7) * t73 + t233 * t32) + (pkin(4) * t287 - t228 * t261 + qJ(2)) * t96 - t300 * t42, (t229 * (pkin(4) * t233 - pkin(8) * t237) - t228 * (-pkin(8) * t233 + t261) + qJ(2)) * t40 + t313 * t21, t229 * (t237 * t48 + t265) - t228 * (t233 * t48 - t263), t229 * (t233 * t127 + t237 * t26) - t228 * (-t237 * t127 + t233 * t26), t229 * (t233 * t84 + t237 * t61) - t228 * (t233 * t61 - t237 * t84), t229 * (t237 * t47 - t265) - t228 * (t233 * t47 + t263), t229 * (-t233 * t81 + t237 * t62) - t228 * (t233 * t62 + t237 * t81), t229 * (t233 * t176 + t237 * t69) - t228 * (-t237 * t176 + t233 * t69), t229 * (-pkin(7) * t33 + t237 * t13 - t233 * t17) - t228 * (-pkin(3) * t45 + pkin(7) * t34 + t233 * t13 + t237 * t17) + qJ(2) * t45 - t300 * t19, t229 * (-pkin(7) * t38 + t237 * t16 - t233 * t18) - t228 * (-pkin(3) * t57 + pkin(7) * t39 + t233 * t16 + t237 * t18) + qJ(2) * t57 - t300 * t22, t229 * (-pkin(7) * t23 - t233 * t20 + t237 * t3) - t228 * (-pkin(3) * t25 + pkin(7) * t24 + t237 * t20 + t233 * t3) + qJ(2) * t25 - t300 * t10, t229 * (-pkin(7) * t5 + t237 * t2 - t233 * t4) - t228 * (-pkin(3) * t7 + pkin(7) * t6 + t233 * t2 + t237 * t4) + qJ(2) * t7 - t300 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t239, -t206, 0, 0, 0, 0, 0, 0, -t214, -t213, -t215, t143, 0, 0, 0, 0, 0, 0, t107, t124, t110, t64, 0, 0, 0, 0, 0, 0, t49, t54, t42, t21, 0, 0, 0, 0, 0, 0, t19, t22, t10, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t268, t267, -t312, -t197, 0, 0, 0, 0, 0, 0, t185, t187, t168, -t175, 0, 0, 0, 0, 0, 0, t108, t113, t96, t40, 0, 0, 0, 0, 0, 0, t45, t57, t25, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t289, t208 - t207, t209, -t289, -t308, qJDD(4), -t136, -t137, 0, 0, t232 * t155 + t195 * t290, t130 * t232 + t132 * t236, t173 * t236 + t318, t193 * t291 - t236 * t250, t172 * t232 + t279, (-t193 * t232 - t195 * t236) * t205, pkin(4) * t130 + pkin(8) * t109 - t280, pkin(4) * t134 + pkin(8) * t114 + t285, pkin(4) * t150 + pkin(8) * t98 + t41, -pkin(4) * t105 + pkin(8) * t41, t232 * t79 + t236 * t78, t232 * t52 + t236 * t50, t232 * t94 + t236 * t92, t232 * t77 + t236 * t76, t232 * t95 + t236 * t93, t111 * t236 + t112 * t232, -pkin(4) * t80 + pkin(8) * t46 + t232 * t43 + t236 * t31, -pkin(4) * t310 + pkin(8) * t58 + t232 * t44 + t236 * t35, -pkin(4) * t103 + pkin(8) * t27 + t11 * t232 + t236 * t9, -pkin(4) * t65 + pkin(8) * t8 - pkin(9) * t297 + t12 * t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t169, t167, t133, -t169, -t129, t177, -t70, -t71, 0, 0, t128, t127, t84, -t128, -t81, t176, t253, t247, t302, t303; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, t127, t84, -t128, -t81, t176, -t29, -t30, 0, 0;];
tauJ_reg  = t28;
