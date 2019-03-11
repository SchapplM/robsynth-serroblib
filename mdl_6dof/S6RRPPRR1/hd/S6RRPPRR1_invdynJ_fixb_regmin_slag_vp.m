% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% 
% Output:
% tau_reg [6x30]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPPRR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:48:07
% EndTime: 2019-03-09 08:48:16
% DurationCPUTime: 4.12s
% Computational Cost: add. (5076->415), mult. (11788->515), div. (0->0), fcn. (9070->12), ass. (0->218)
t183 = qJDD(2) - qJDD(5);
t184 = qJD(2) - qJD(5);
t196 = sin(qJ(6));
t200 = cos(qJ(6));
t192 = sin(pkin(10));
t193 = cos(pkin(10));
t202 = cos(qJ(2));
t271 = qJD(1) * t202;
t257 = t193 * t271;
t198 = sin(qJ(2));
t272 = qJD(1) * t198;
t131 = t192 * t272 - t257;
t147 = t192 * t202 + t193 * t198;
t134 = t147 * qJD(1);
t197 = sin(qJ(5));
t201 = cos(qJ(5));
t268 = qJD(5) * t201;
t269 = qJD(5) * t197;
t133 = t147 * qJD(2);
t263 = t202 * qJDD(1);
t264 = t198 * qJDD(1);
t230 = t192 * t264 - t193 * t263;
t89 = qJD(1) * t133 + t230;
t265 = qJD(1) * qJD(2);
t256 = t198 * t265;
t214 = qJDD(1) * t147 - t192 * t256;
t255 = t202 * t265;
t90 = t193 * t255 + t214;
t242 = -t131 * t268 + t134 * t269 - t197 * t89 - t201 * t90;
t266 = qJD(6) * t200;
t267 = qJD(6) * t196;
t314 = t131 * t197 + t134 * t201;
t10 = -t183 * t196 - t184 * t266 - t200 * t242 - t267 * t314;
t229 = t184 * t196 - t200 * t314;
t11 = -qJD(6) * t229 + t183 * t200 - t196 * t242;
t317 = -t131 * t201 + t134 * t197;
t321 = qJD(6) + t317;
t296 = t229 * t321;
t56 = t184 * t200 + t196 * t314;
t298 = t56 * t321;
t338 = t200 * (t10 - t298) + (-t11 + t296) * t196;
t295 = t229 * t314;
t28 = qJD(5) * t314 + t197 * t90 - t201 * t89;
t26 = qJDD(6) + t28;
t289 = t196 * t26;
t331 = t200 * t321;
t328 = t321 * t331 + t289;
t335 = t295 + t328;
t292 = t10 * t196;
t334 = t229 * t331 - t292;
t291 = t184 * t314;
t332 = t28 + t291;
t195 = -qJ(3) - pkin(7);
t158 = t195 * t198;
t151 = qJD(1) * t158;
t159 = t195 * t202;
t152 = qJD(1) * t159;
t280 = t192 * t152;
t95 = t151 * t193 + t280;
t276 = qJD(4) - t95;
t249 = -t200 * t26 + t267 * t321;
t288 = t196 * t321;
t240 = t288 * t317 + t249;
t297 = t56 * t314;
t330 = t240 - t297;
t290 = t184 * t317;
t329 = t242 + t290;
t327 = pkin(5) * t314;
t155 = -qJD(1) * pkin(1) - pkin(2) * t271 + qJD(3);
t61 = pkin(3) * t131 - qJ(4) * t134 + t155;
t45 = -pkin(4) * t131 - t61;
t12 = pkin(5) * t317 - pkin(9) * t314 + t45;
t143 = qJD(2) * pkin(2) + t151;
t87 = t143 * t193 + t280;
t233 = qJD(4) - t87;
t305 = pkin(8) * t134;
t307 = -pkin(3) - pkin(4);
t48 = qJD(2) * t307 + t233 - t305;
t306 = pkin(8) * t131;
t279 = t193 * t152;
t88 = t143 * t192 - t279;
t77 = qJD(2) * qJ(4) + t88;
t54 = t77 + t306;
t19 = t197 * t48 + t201 * t54;
t16 = -pkin(9) * t184 + t19;
t5 = t12 * t200 - t16 * t196;
t326 = t314 * t5;
t6 = t12 * t196 + t16 * t200;
t325 = t314 * t6;
t324 = t321 * t314;
t323 = t314 * t317;
t322 = -t305 + t276;
t320 = t314 ^ 2 - t317 ^ 2;
t203 = cos(qJ(1));
t185 = qJ(2) + pkin(10);
t179 = sin(t185);
t180 = cos(t185);
t224 = t179 * t197 + t180 * t201;
t105 = t224 * t203;
t248 = qJD(2) * t195;
t127 = -qJD(3) * t198 + t202 * t248;
t84 = qJDD(2) * pkin(2) + qJD(1) * t127 + qJDD(1) * t158;
t126 = qJD(3) * t202 + t198 * t248;
t93 = qJD(1) * t126 - qJDD(1) * t159;
t38 = -t192 * t93 + t193 * t84;
t258 = -qJDD(4) + t38;
t24 = -pkin(8) * t90 + qJDD(2) * t307 - t258;
t187 = qJD(2) * qJD(4);
t39 = t192 * t84 + t193 * t93;
t262 = qJDD(2) * qJ(4) + t39;
t34 = t187 + t262;
t25 = pkin(8) * t89 + t34;
t222 = t197 * t24 + t201 * t25 + t268 * t48 - t269 * t54;
t199 = sin(qJ(1));
t103 = t224 * t199;
t313 = -t179 * t201 + t180 * t197;
t235 = g(2) * t103 - g(3) * t313;
t319 = g(1) * t105 + t317 * t45 - t222 + t235;
t102 = t313 * t199;
t104 = t313 * t203;
t221 = g(1) * t104 + g(2) * t102 + g(3) * t224;
t228 = t197 * t25 - t201 * t24 + t268 * t54 + t269 * t48;
t318 = -t314 * t45 + t221 - t228;
t129 = t134 ^ 2;
t316 = -t131 ^ 2 - t129;
t175 = -pkin(2) * t193 - pkin(3);
t170 = -pkin(4) + t175;
t173 = pkin(2) * t192 + qJ(4);
t274 = t170 * t197 + t173 * t201;
t182 = g(2) * t203;
t315 = g(1) * t199 - t182;
t302 = g(1) * t203;
t238 = g(2) * t199 + t302;
t312 = pkin(3) * t89 - qJ(4) * t90 - qJD(4) * t134;
t107 = -pkin(9) + t274;
t4 = pkin(5) * t183 + t228;
t219 = t221 - t4;
t63 = pkin(2) * t272 + pkin(3) * t134 + qJ(4) * t131;
t49 = -pkin(4) * t134 - t63;
t311 = (-pkin(9) * t317 + qJD(6) * t107 - t327 + t49) * t321 + t219;
t310 = (pkin(9) * t321 + t327) * t321 - t219;
t309 = -g(3) * t180 - t134 * t61 + t179 * t238 + t258;
t18 = -t197 * t54 + t201 * t48;
t15 = pkin(5) * t184 - t18;
t278 = t193 * t202;
t146 = t192 * t198 - t278;
t226 = t146 * t201 - t147 * t197;
t181 = t202 * pkin(2);
t176 = t181 + pkin(1);
t86 = pkin(3) * t146 - qJ(4) * t147 - t176;
t55 = -pkin(4) * t146 - t86;
t92 = t146 * t197 + t147 * t201;
t22 = -pkin(5) * t226 - pkin(9) * t92 + t55;
t252 = -pkin(9) * t183 + qJD(6) * t12 + t222;
t100 = -t158 * t193 - t159 * t192;
t66 = -pkin(8) * t147 + t100;
t101 = t158 * t192 - t159 * t193;
t67 = pkin(8) * t146 + t101;
t33 = t197 * t66 + t201 * t67;
t270 = qJD(2) * t198;
t136 = qJD(2) * t278 - t192 * t270;
t40 = qJD(5) * t226 + t133 * t197 + t136 * t201;
t32 = t197 * t67 - t201 * t66;
t64 = t126 * t192 - t127 * t193;
t50 = -pkin(8) * t136 + t64;
t65 = t126 * t193 + t127 * t192;
t51 = pkin(8) * t133 + t65;
t8 = -qJD(5) * t32 + t197 * t50 + t201 * t51;
t308 = t15 * t40 - (qJD(6) * t22 + t8) * t321 + t252 * t226 - t33 * t26 + t4 * t92 + t302;
t304 = g(1) * t103;
t300 = g(3) * t202;
t299 = t15 * t92;
t225 = t170 * t201 - t173 * t197;
t94 = t151 * t192 - t279;
t59 = t94 + t306;
t294 = t225 * qJD(5) - t197 * t59 + t201 * t322;
t293 = t274 * qJD(5) + t197 * t322 + t201 * t59;
t286 = qJD(6) * t16;
t191 = qJDD(1) * pkin(1);
t285 = qJDD(2) * pkin(3);
t189 = t198 ^ 2;
t273 = -t202 ^ 2 + t189;
t261 = pkin(2) * t270;
t260 = t92 * t267;
t259 = t321 * t266;
t243 = t184 * t321;
t241 = t184 ^ 2;
t239 = qJDD(3) - t191 + (t256 - t263) * pkin(2);
t236 = t22 * t26 + t304;
t234 = t26 * t92 + t321 * t40;
t232 = pkin(3) * t180 + qJ(4) * t179;
t231 = -t286 - t182;
t223 = -0.2e1 * pkin(1) * t265 - pkin(7) * qJDD(2);
t53 = pkin(3) * t133 - qJ(4) * t136 - qJD(4) * t147 + t261;
t29 = t239 + t312;
t218 = t239 - t315;
t216 = t235 - t252;
t215 = -pkin(9) * t26 + (t15 + t18) * t321;
t42 = -pkin(4) * t133 - t53;
t13 = -pkin(4) * t89 - t29;
t204 = qJD(2) ^ 2;
t212 = -pkin(7) * t204 + 0.2e1 * t191 + t315;
t205 = qJD(1) ^ 2;
t211 = pkin(1) * t205 - pkin(7) * qJDD(1) + t238;
t210 = -t107 * t26 + (-t15 - t294) * t321;
t208 = t100 * t90 - t101 * t89 - t131 * t65 + t134 * t64 - t238;
t161 = t203 * t176;
t106 = pkin(5) - t225;
t97 = t105 * t200 - t196 * t199;
t96 = -t105 * t196 - t199 * t200;
t69 = -qJD(2) * pkin(3) + t233;
t41 = qJD(5) * t92 - t133 * t201 + t136 * t197;
t37 = -t258 - t285;
t9 = qJD(5) * t33 + t197 * t51 - t201 * t50;
t7 = pkin(5) * t41 - pkin(9) * t40 + t42;
t2 = pkin(5) * t28 + pkin(9) * t242 + t13;
t1 = t200 * t2;
t3 = [qJDD(1), t315, t238, qJDD(1) * t189 + 0.2e1 * t198 * t255, 0.2e1 * t198 * t263 - 0.2e1 * t265 * t273, qJDD(2) * t198 + t202 * t204, qJDD(2) * t202 - t198 * t204, 0, t198 * t223 + t202 * t212, -t198 * t212 + t202 * t223, -t133 * t88 - t136 * t87 - t146 * t39 - t147 * t38 + t208, t39 * t101 + t88 * t65 - t38 * t100 - t87 * t64 - t239 * t176 + t155 * t261 - g(1) * (-t176 * t199 - t195 * t203) - g(2) * (-t195 * t199 + t161) -qJD(2) * t64 - qJDD(2) * t100 + t131 * t53 + t133 * t61 + t146 * t29 + t180 * t315 + t86 * t89, -t133 * t77 + t136 * t69 - t146 * t34 + t147 * t37 + t208, qJD(2) * t65 + qJDD(2) * t101 - t134 * t53 - t136 * t61 - t147 * t29 + t179 * t315 - t86 * t90, -g(2) * t161 + t37 * t100 + t34 * t101 + t29 * t86 + t61 * t53 + t69 * t64 + t77 * t65 + (g(1) * t195 - g(2) * t232) * t203 + (-g(1) * (-t176 - t232) + g(2) * t195) * t199, -t242 * t92 + t314 * t40, -t226 * t242 - t28 * t92 - t314 * t41 - t317 * t40, -t183 * t92 - t184 * t40, -t183 * t226 + t184 * t41, 0, -g(2) * t105 - t13 * t226 + t183 * t32 + t184 * t9 + t28 * t55 + t317 * t42 + t41 * t45 + t304, -g(1) * t102 + g(2) * t104 + t13 * t92 + t183 * t33 + t184 * t8 - t242 * t55 + t314 * t42 + t40 * t45, t229 * t260 + (t10 * t92 - t229 * t40) * t200 (t196 * t229 - t200 * t56) * t40 + (-t292 - t11 * t200 + (t196 * t56 + t200 * t229) * qJD(6)) * t92, -t10 * t226 + t200 * t234 - t229 * t41 - t260 * t321, t11 * t226 - t196 * t234 - t259 * t92 - t41 * t56, -t226 * t26 + t321 * t41, -g(2) * t97 - t1 * t226 + t32 * t11 + t5 * t41 + t9 * t56 + (t7 * t321 + (t16 * t226 - t321 * t33 + t299) * qJD(6) + t236) * t200 + t308 * t196, -g(2) * t96 + t32 * t10 - t6 * t41 - t9 * t229 + (-(-qJD(6) * t33 + t7) * t321 + (t2 - t286) * t226 - qJD(6) * t299 - t236) * t196 + t308 * t200; 0, 0, 0, -t198 * t205 * t202, t273 * t205, t264, t263, qJDD(2), t198 * t211 - t300, g(3) * t198 + t202 * t211 (t88 - t94) * t134 + (-t87 + t95) * t131 + (-t192 * t89 - t193 * t90) * pkin(2), t87 * t94 - t88 * t95 + (-t300 + t192 * t39 + t193 * t38 + (-qJD(1) * t155 + t238) * t198) * pkin(2), qJD(2) * t94 - t131 * t63 + (pkin(3) - t175) * qJDD(2) + t309, -t173 * t89 + t175 * t90 + (t77 - t94) * t134 + (t69 - t276) * t131, -g(3) * t179 - qJD(2) * t95 + qJDD(2) * t173 - t131 * t61 + t134 * t63 - t180 * t238 + 0.2e1 * t187 + t262, t34 * t173 + t37 * t175 - t61 * t63 - t69 * t94 - g(3) * (t181 + t232) + t276 * t77 + t238 * (pkin(2) * t198 + pkin(3) * t179 - qJ(4) * t180) -t323, -t320, t329, t332, t183, -t183 * t225 + t184 * t293 - t317 * t49 - t318, t183 * t274 + t184 * t294 - t314 * t49 - t319, t334, -t338, -t335, t330, t324, t106 * t11 + t196 * t210 - t200 * t311 + t293 * t56 + t326, t106 * t10 + t196 * t311 + t200 * t210 - t229 * t293 - t325; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t316, t131 * t88 + t134 * t87 + t218, 0.2e1 * qJD(2) * t134 + t230, t316 (t131 - t257) * qJD(2) - t214, t131 * t77 - t134 * t69 + t218 + t312, 0, 0, 0, 0, 0, -t28 + t291, t242 - t290, 0, 0, 0, 0, 0, t240 + t297, -t295 + t328; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131 * t134 - qJDD(2) (t131 + t257) * qJD(2) + t214, -t129 - t204, -qJD(2) * t77 - t285 - t309, 0, 0, 0, 0, 0, -t134 * t317 - t183 * t201 - t197 * t241, -t134 * t314 + t183 * t197 - t201 * t241, 0, 0, 0, 0, 0, -t134 * t331 + (t196 * t243 - t11) * t201 + (-t184 * t56 - t259 - t289) * t197, t134 * t288 + (t200 * t243 - t10) * t201 + (t184 * t229 + t249) * t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t323, t320, -t329, -t332, -t183, -t184 * t19 + t318, -t18 * t184 + t319, -t334, t338, t335, -t330, -t324, -pkin(5) * t11 - t19 * t56 + t196 * t215 - t200 * t310 - t326, -pkin(5) * t10 + t19 * t229 + t196 * t310 + t200 * t215 + t325; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t229 * t56, t229 ^ 2 - t56 ^ 2, t10 + t298, -t11 - t296, t26, -g(1) * t96 + t15 * t229 + t196 * t216 + t200 * t231 + t321 * t6 + t1, g(1) * t97 + t15 * t56 + t5 * t321 + (-t2 - t231) * t196 + t216 * t200;];
tau_reg  = t3;
