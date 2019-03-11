% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRR8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR8_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR8_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:36:03
% EndTime: 2019-03-08 22:36:28
% DurationCPUTime: 11.84s
% Computational Cost: add. (7795->664), mult. (20946->921), div. (0->0), fcn. (15924->12), ass. (0->325)
t216 = sin(pkin(7));
t222 = sin(qJ(3));
t311 = qJD(3) * t222;
t291 = t216 * t311;
t206 = pkin(3) * t291;
t226 = cos(qJ(3));
t259 = pkin(10) * t222 - qJ(4) * t226;
t309 = qJD(4) * t222;
t235 = (qJD(3) * t259 - t309) * t216;
t117 = t206 + t235;
t334 = t216 * t222;
t210 = pkin(9) * t334;
t218 = cos(pkin(7));
t298 = -pkin(2) * t226 - pkin(3);
t120 = pkin(4) * t334 + t210 + (-pkin(10) + t298) * t218;
t228 = -pkin(3) - pkin(10);
t335 = qJ(4) * t222;
t247 = t226 * t228 - t335;
t132 = (-pkin(2) + t247) * t216;
t330 = t218 * t222;
t213 = pkin(2) * t330;
t333 = t216 * t226;
t374 = pkin(4) + pkin(9);
t138 = (t333 * t374 + t213) * qJD(3);
t223 = sin(qJ(2));
t324 = t223 * t226;
t227 = cos(qJ(2));
t325 = t222 * t227;
t245 = t218 * t324 + t325;
t217 = sin(pkin(6));
t318 = qJD(1) * t217;
t140 = t245 * t318;
t221 = sin(qJ(5));
t225 = cos(qJ(5));
t297 = t223 * t318;
t284 = t216 * t297;
t307 = qJD(5) * t225;
t308 = qJD(5) * t221;
t391 = t120 * t307 - t132 * t308 + (t117 - t284) * t225 + (t138 - t140) * t221;
t310 = qJD(3) * t226;
t289 = t218 * t310;
t207 = pkin(2) * t289;
t214 = t218 * qJD(4);
t119 = -t291 * t374 + t207 + t214;
t296 = t227 * t318;
t282 = t226 * t296;
t283 = t218 * t297;
t141 = -t222 * t283 + t282;
t403 = t119 - t141;
t290 = t216 * t310;
t402 = -pkin(11) * t290 - t391;
t173 = t218 * t221 + t225 * t333;
t288 = t221 * t311;
t126 = -qJD(5) * t173 + t216 * t288;
t303 = t221 * t333;
t127 = -qJD(5) * t303 + t218 * t307 - t225 * t291;
t401 = pkin(5) * t127 - pkin(11) * t126 + t403;
t316 = qJD(2) * t216;
t294 = t222 * t316;
t199 = qJD(5) + t294;
t209 = qJD(2) * t218 + qJD(3);
t177 = pkin(9) * t316 + t297;
t187 = qJD(2) * pkin(2) + t296;
t219 = cos(pkin(6));
t300 = t219 * t333;
t189 = qJD(1) * t300;
t329 = t218 * t226;
t320 = t187 * t329 + t189;
t86 = -t222 * (pkin(4) * t316 + t177) + t320;
t65 = t209 * t228 + qJD(4) - t86;
t317 = qJD(1) * t219;
t208 = t218 * t317;
t93 = t208 + (qJD(2) * t247 - t187) * t216;
t34 = -t221 * t93 + t225 * t65;
t24 = -pkin(5) * t199 - t34;
t220 = sin(qJ(6));
t224 = cos(qJ(6));
t35 = t221 * t65 + t225 * t93;
t25 = pkin(11) * t199 + t35;
t314 = qJD(2) * t226;
t293 = t216 * t314;
t152 = -t209 * t221 - t225 * t293;
t278 = t221 * t293;
t153 = t209 * t225 - t278;
t202 = t209 * qJ(4);
t246 = t226 * t177 + t187 * t330;
t295 = t222 * t317;
t99 = t216 * t295 + t246;
t87 = pkin(4) * t293 + t99;
t74 = t202 + t87;
t36 = -pkin(5) * t152 - pkin(11) * t153 + t74;
t11 = -t220 * t25 + t224 * t36;
t12 = t220 * t36 + t224 * t25;
t258 = t11 * t224 + t12 * t220;
t267 = mrSges(7,1) * t220 + mrSges(7,2) * t224;
t362 = t224 / 0.2e1;
t363 = -t220 / 0.2e1;
t105 = -t153 * t220 + t199 * t224;
t149 = qJD(6) - t152;
t106 = t153 * t224 + t199 * t220;
t356 = Ifges(7,4) * t106;
t38 = Ifges(7,2) * t105 + Ifges(7,6) * t149 + t356;
t104 = Ifges(7,4) * t105;
t39 = Ifges(7,1) * t106 + Ifges(7,5) * t149 + t104;
t400 = -t258 * mrSges(7,3) + t24 * t267 + t362 * t39 + t363 * t38;
t396 = -t209 / 0.2e1;
t399 = -t293 / 0.2e1;
t392 = mrSges(5,1) + mrSges(4,3);
t398 = mrSges(5,2) - mrSges(4,1);
t262 = Ifges(7,5) * t224 - Ifges(7,6) * t220;
t354 = Ifges(7,4) * t224;
t264 = -Ifges(7,2) * t220 + t354;
t355 = Ifges(7,4) * t220;
t266 = Ifges(7,1) * t224 - t355;
t366 = t149 / 0.2e1;
t370 = t106 / 0.2e1;
t372 = t105 / 0.2e1;
t397 = t262 * t366 + t264 * t372 + t266 * t370 + t400;
t395 = -t316 / 0.2e1;
t388 = t221 * t120 + t225 * t132;
t69 = pkin(11) * t334 + t388;
t176 = pkin(9) * t333 + t213;
t156 = -t218 * qJ(4) - t176;
t131 = pkin(4) * t333 - t156;
t174 = t218 * t225 - t303;
t75 = pkin(5) * t173 - pkin(11) * t174 + t131;
t26 = -t220 * t69 + t224 * t75;
t394 = qJD(6) * t26 + t401 * t220 - t224 * t402;
t27 = t220 * t75 + t224 * t69;
t393 = -qJD(6) * t27 + t220 * t402 + t401 * t224;
t315 = qJD(2) * t217;
t292 = t223 * t315;
t277 = qJD(1) * t292;
t313 = qJD(3) * t177;
t390 = t222 * (t218 * t277 + t313);
t98 = t177 * t222 - t320;
t389 = -qJD(4) - t98;
t168 = t176 * qJD(3);
t387 = t140 - t168;
t147 = Ifges(6,4) * t152;
t353 = Ifges(6,5) * t199;
t358 = Ifges(6,1) * t153;
t79 = t147 + t353 + t358;
t386 = t34 * mrSges(6,3) - t79 / 0.2e1 - t353 / 0.2e1 - t74 * mrSges(6,2) - t147 / 0.2e1;
t111 = -t209 * t308 + (-t226 * t307 + t288) * t316;
t312 = qJD(3) * t216;
t286 = qJD(2) * t312;
t276 = t222 * t286;
t112 = -qJD(5) * t278 + t209 * t307 - t225 * t276;
t299 = qJD(2) * t282 + qJD(3) * t189 + t187 * t289;
t285 = -t209 * qJD(4) - t299;
t42 = (-t313 + (-pkin(4) * t312 - t283) * qJD(2)) * t222 - t285;
t21 = pkin(5) * t112 - pkin(11) * t111 + t42;
t275 = t226 * t286;
t242 = t245 * qJD(2);
t238 = t217 * t242;
t53 = qJD(1) * t238 + ((pkin(4) * t314 + t295) * t216 + t246) * qJD(3);
t319 = pkin(3) * t276 + t216 * t277;
t90 = qJD(2) * t235 + t319;
t9 = t221 * t53 + t225 * t90 + t65 * t307 - t308 * t93;
t7 = pkin(11) * t275 + t9;
t1 = qJD(6) * t11 + t21 * t220 + t224 * t7;
t2 = -qJD(6) * t12 + t21 * t224 - t220 * t7;
t270 = t1 * t224 - t2 * t220;
t10 = -qJD(5) * t35 - t221 * t90 + t225 * t53;
t385 = t10 * mrSges(6,1) - t9 * mrSges(6,2) + Ifges(6,5) * t111 - Ifges(6,6) * t112;
t384 = Ifges(4,4) * t399 + Ifges(4,5) * t396;
t29 = -qJD(5) * t388 - t117 * t221 + t138 * t225;
t382 = Ifges(4,5) / 0.2e1;
t381 = Ifges(5,5) / 0.2e1;
t380 = Ifges(6,2) / 0.2e1;
t49 = -qJD(6) * t106 - t111 * t220 + t224 * t275;
t43 = Ifges(7,6) * t49;
t48 = qJD(6) * t105 + t111 * t224 + t220 * t275;
t44 = Ifges(7,5) * t48;
t13 = Ifges(7,3) * t112 + t43 + t44;
t379 = t13 / 0.2e1;
t15 = Ifges(7,1) * t48 + Ifges(7,4) * t49 + Ifges(7,5) * t112;
t378 = t15 / 0.2e1;
t377 = -t38 / 0.2e1;
t376 = t48 / 0.2e1;
t375 = t49 / 0.2e1;
t373 = -t105 / 0.2e1;
t371 = -t106 / 0.2e1;
t369 = t112 / 0.2e1;
t367 = -t149 / 0.2e1;
t365 = -t173 / 0.2e1;
t364 = t174 / 0.2e1;
t16 = -mrSges(7,1) * t49 + mrSges(7,2) * t48;
t94 = mrSges(6,1) * t275 - mrSges(6,3) * t111;
t359 = t94 - t16;
t357 = Ifges(6,4) * t153;
t352 = Ifges(7,5) * t106;
t351 = Ifges(6,2) * t152;
t350 = Ifges(6,6) * t199;
t349 = Ifges(7,6) * t105;
t348 = Ifges(7,3) * t149;
t347 = t111 * Ifges(6,1);
t346 = t111 * Ifges(6,4);
t345 = t112 * Ifges(6,4);
t332 = t217 * t223;
t301 = t222 * t332;
t331 = t217 * t227;
t302 = t218 * t331;
t121 = -t226 * t302 - t300 + t301;
t279 = t219 * t291;
t61 = t246 * qJD(3) + (t238 + t279) * qJD(1);
t344 = t121 * t61;
t339 = t209 * Ifges(5,5);
t338 = t221 * t24;
t204 = pkin(3) * t294;
t137 = t259 * t316 + t204;
t52 = t225 * t137 + t221 * t87;
t115 = mrSges(6,1) * t199 - mrSges(6,3) * t153;
t57 = -mrSges(7,1) * t105 + mrSges(7,2) * t106;
t337 = -t57 + t115;
t162 = -mrSges(5,1) * t293 - mrSges(5,3) * t209;
t96 = -mrSges(6,1) * t152 + mrSges(6,2) * t153;
t336 = t96 - t162;
t328 = t220 * t222;
t327 = t221 * t228;
t326 = t222 * t224;
t323 = -t209 * t398 - t294 * t392;
t161 = -mrSges(4,2) * t209 + mrSges(4,3) * t293;
t322 = t161 - t162;
t164 = (mrSges(5,2) * t226 - mrSges(5,3) * t222) * t316;
t269 = -mrSges(4,1) * t226 + mrSges(4,2) * t222;
t321 = t269 * t316 + t164;
t306 = t209 / 0.2e1 - qJD(3);
t305 = 0.3e1 / 0.2e1 * Ifges(5,6) + 0.3e1 / 0.2e1 * Ifges(4,4);
t304 = t96 + t322;
t287 = t228 * t307;
t281 = t216 * t292;
t280 = t225 * t294;
t274 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t271 = pkin(5) * t225 + pkin(11) * t221;
t268 = mrSges(7,1) * t224 - mrSges(7,2) * t220;
t265 = Ifges(7,1) * t220 + t354;
t263 = Ifges(7,2) * t224 + t355;
t261 = Ifges(7,5) * t220 + Ifges(7,6) * t224;
t260 = -pkin(3) * t226 - t335;
t257 = t11 * t220 - t12 * t224;
t32 = mrSges(7,1) * t112 - mrSges(7,3) * t48;
t33 = -mrSges(7,2) * t112 + mrSges(7,3) * t49;
t256 = -t220 * t32 + t224 * t33;
t72 = -mrSges(7,2) * t149 + mrSges(7,3) * t105;
t73 = mrSges(7,1) * t149 - mrSges(7,3) * t106;
t255 = -t220 * t72 - t224 * t73;
t254 = t221 * t34 - t225 * t35;
t244 = t218 * t325 + t324;
t122 = t217 * t244 + t219 * t334;
t172 = -t216 * t331 + t218 * t219;
t89 = t121 * t221 + t172 * t225;
t45 = t122 * t224 - t220 * t89;
t46 = t122 * t220 + t224 * t89;
t51 = -t137 * t221 + t225 * t87;
t70 = t120 * t225 - t132 * t221;
t249 = t121 * t225 - t172 * t221;
t167 = -pkin(9) * t291 + t207;
t128 = -t174 * t220 + t216 * t326;
t129 = t174 * t224 + t216 * t328;
t190 = pkin(5) * t221 - pkin(11) * t225 + qJ(4);
t151 = t190 * t220 + t224 * t327;
t150 = t190 * t224 - t220 * t327;
t241 = (-qJ(4) * t310 - t309) * t216;
t54 = t285 + t390;
t60 = t299 - t390;
t240 = -t60 * mrSges(4,2) - t54 * mrSges(5,3) + t398 * t61;
t103 = t208 + (qJD(2) * t260 - t187) * t216;
t146 = -t187 * t216 + t208;
t236 = t146 * mrSges(4,1) + Ifges(4,6) * t396 + (Ifges(4,4) * t222 + t226 * Ifges(4,2)) * t395 + t339 / 0.2e1 + (-Ifges(5,6) * t222 - t226 * Ifges(5,3)) * t316 / 0.2e1 - t103 * mrSges(5,2) - t99 * mrSges(4,3);
t234 = -t358 / 0.2e1 + t386;
t233 = t146 * mrSges(4,2) + t34 * mrSges(6,1) + t98 * mrSges(4,3) + Ifges(4,1) * t294 / 0.2e1 + Ifges(5,4) * t396 + (-Ifges(5,2) * t222 - Ifges(5,6) * t226) * t395 + t199 * Ifges(6,3) + t153 * Ifges(6,5) + t152 * Ifges(6,6) - t103 * mrSges(5,3) - t35 * mrSges(6,2) - t384;
t37 = t348 + t349 + t352;
t78 = t350 + t351 + t357;
t232 = t12 * mrSges(7,2) - t11 * mrSges(7,1) - t74 * mrSges(6,1) + t350 / 0.2e1 - t349 / 0.2e1 - t352 / 0.2e1 - t348 / 0.2e1 + t35 * mrSges(6,3) - t37 / 0.2e1 + t78 / 0.2e1 + t357 / 0.2e1;
t231 = (-t351 / 0.2e1 - t232) * t225;
t197 = Ifges(4,5) * t275;
t196 = Ifges(5,5) * t276;
t195 = Ifges(6,3) * t275;
t178 = qJD(5) * t271 + qJD(4);
t175 = pkin(2) * t329 - t210;
t166 = -qJ(4) * t293 + t204;
t158 = t218 * t298 + t210;
t157 = (-pkin(2) + t260) * t216;
t155 = (mrSges(4,1) * t222 + mrSges(4,2) * t226) * t286;
t154 = (-mrSges(5,2) * t222 - mrSges(5,3) * t226) * t286;
t148 = -t167 - t214;
t145 = (t220 * t226 + t221 * t326) * t316;
t144 = (-t221 * t328 + t224 * t226) * t316;
t143 = t209 * t220 - t224 * t280;
t142 = t209 * t224 + t220 * t280;
t139 = t206 + t241;
t114 = -mrSges(6,2) * t199 + mrSges(6,3) * t152;
t101 = -t225 * t140 + t221 * t284;
t100 = qJD(2) * t241 + t319;
t97 = pkin(5) * t153 - pkin(11) * t152;
t95 = -mrSges(6,2) * t275 - mrSges(6,3) * t112;
t92 = -qJD(6) * t151 + t178 * t224 - t220 * t287;
t91 = qJD(6) * t150 + t178 * t220 + t224 * t287;
t84 = -t202 - t99;
t83 = -pkin(3) * t209 - t389;
t81 = -t292 * t330 - qJD(3) * t301 + (t227 * t315 + (t216 * t219 + t302) * qJD(3)) * t226;
t80 = t279 + (qJD(3) * t244 + t242) * t217;
t68 = -pkin(5) * t334 - t70;
t64 = -qJD(6) * t129 - t126 * t220 + t224 * t290;
t63 = qJD(6) * t128 + t126 * t224 + t220 * t290;
t62 = mrSges(6,1) * t112 + mrSges(6,2) * t111;
t59 = (-t177 + (-pkin(4) - t271) * t316) * t222 + t320;
t56 = Ifges(6,5) * t275 - t345 + t347;
t55 = -t112 * Ifges(6,2) + Ifges(6,6) * t275 + t346;
t41 = pkin(11) * t293 + t52;
t40 = -pkin(5) * t293 - t51;
t31 = qJD(5) * t249 + t221 * t80 + t225 * t281;
t30 = qJD(5) * t89 + t221 * t281 - t80 * t225;
t23 = -pkin(5) * t290 - t29;
t20 = t220 * t59 + t224 * t41;
t19 = -t220 * t41 + t224 * t59;
t18 = t220 * t97 + t224 * t34;
t17 = -t220 * t34 + t224 * t97;
t14 = Ifges(7,4) * t48 + Ifges(7,2) * t49 + Ifges(7,6) * t112;
t8 = -pkin(5) * t275 - t10;
t6 = qJD(6) * t45 + t220 * t81 + t224 * t31;
t5 = -qJD(6) * t46 - t220 * t31 + t224 * t81;
t3 = [t31 * t114 + t122 * t62 + t45 * t32 + t46 * t33 + t5 * t73 + t6 * t72 + t89 * t95 + t359 * t249 - t323 * t80 - t337 * t30 + (-mrSges(3,1) * t223 - mrSges(3,2) * t227) * qJD(2) ^ 2 * t217 + (t154 + t155) * t172 + t304 * t81 + (t321 * t332 + t392 * qJD(3) * (t121 * t226 - t122 * t222)) * t316 + m(4) * (t344 + t122 * t60 + t80 * t98 + t81 * t99 + (qJD(1) * t172 + t146) * t281) + m(5) * (t100 * t172 + t103 * t281 - t122 * t54 + t80 * t83 - t81 * t84 + t344) + m(6) * (t10 * t249 + t122 * t42 - t30 * t34 + t31 * t35 + t74 * t81 + t89 * t9) + m(7) * (t1 * t46 + t11 * t5 + t12 * t6 + t2 * t45 + t24 * t30 - t249 * t8); (t10 * t70 + t131 * t42 + t388 * t9 + t403 * t74 + t391 * t35 + (t101 + t29) * t34) * m(6) + t388 * t95 + t393 * t73 + (t1 * t27 + t2 * t26 + t68 * t8 + (-t101 + t23) * t24 + t394 * t12 + t393 * t11) * m(7) + t394 * t72 - m(4) * (t140 * t98 + t141 * t99 + t146 * t284) + t337 * t101 + t199 * (Ifges(6,5) * t126 - Ifges(6,6) * t127) / 0.2e1 + t387 * t323 + (t100 * t157 + t156 * t54 + t158 * t61 + (t141 + t148) * t84 - t387 * t83 + (t139 - t284) * t103) * m(5) + (-pkin(2) * t155 - t321 * t297 + (-mrSges(5,1) * t54 + mrSges(5,2) * t100 + mrSges(4,3) * t60) * t226 + (t195 / 0.2e1 - t100 * mrSges(5,3) + t392 * t61 + t385) * t222) * t216 + (-t10 * t174 - t126 * t34 - t127 * t35 - t173 * t9) * mrSges(6,3) - t304 * t141 + m(4) * (t167 * t99 + t168 * t98 - t175 * t61 + t176 * t60) + (t196 / 0.2e1 + t197 / 0.2e1 + t240) * t218 + t111 * (Ifges(6,1) * t174 - Ifges(6,4) * t173) / 0.2e1 + t153 * (Ifges(6,1) * t126 - Ifges(6,4) * t127) / 0.2e1 + ((t84 * mrSges(5,1) + (-Ifges(4,6) / 0.2e1 + t381) * t209 + t236) * t222 + (t83 * mrSges(5,1) + t233 + (t382 - Ifges(5,4) / 0.2e1) * t209) * t226) * t312 + ((-m(4) * pkin(2) + t269) * t284 + ((Ifges(6,5) * t364 + Ifges(6,6) * t365 + t158 * mrSges(5,1) - t175 * mrSges(4,3) + t305 * t333 + (-Ifges(5,4) + t382) * t218) * t226 + (t156 * mrSges(5,1) - t176 * mrSges(4,3) - t305 * t334 + (t381 - Ifges(4,6)) * t218 + (0.3e1 / 0.2e1 * Ifges(5,2) - 0.3e1 / 0.2e1 * Ifges(5,3) + 0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2) + Ifges(6,3) / 0.2e1) * t333) * t222) * qJD(3)) * t316 + t391 * t114 - t112 * (Ifges(6,4) * t174 - Ifges(6,2) * t173) / 0.2e1 + t152 * (Ifges(6,4) * t126 - Ifges(6,2) * t127) / 0.2e1 + t129 * t378 + t173 * t379 + (Ifges(7,5) * t129 + Ifges(7,6) * t128 + Ifges(7,3) * t173) * t369 + (Ifges(7,1) * t63 + Ifges(7,4) * t64 + Ifges(7,5) * t127) * t370 + (Ifges(7,4) * t63 + Ifges(7,2) * t64 + Ifges(7,6) * t127) * t372 + (Ifges(7,4) * t129 + Ifges(7,2) * t128 + Ifges(7,6) * t173) * t375 + (Ifges(7,1) * t129 + Ifges(7,4) * t128 + Ifges(7,5) * t173) * t376 + t56 * t364 + t55 * t365 + (Ifges(7,5) * t63 + Ifges(7,6) * t64 + Ifges(7,3) * t127) * t366 + t26 * t32 + t27 * t33 + t23 * t57 + t63 * t39 / 0.2e1 + t64 * t38 / 0.2e1 + t24 * (-mrSges(7,1) * t64 + mrSges(7,2) * t63) + t68 * t16 + t70 * t94 + t29 * t115 + t119 * t96 + t126 * t79 / 0.2e1 + t127 * t37 / 0.2e1 + t11 * (mrSges(7,1) * t127 - mrSges(7,3) * t63) + t12 * (-mrSges(7,2) * t127 + mrSges(7,3) * t64) - t127 * t78 / 0.2e1 + t74 * (mrSges(6,1) * t127 + mrSges(6,2) * t126) + t128 * t14 / 0.2e1 + t8 * (-mrSges(7,1) * t128 + mrSges(7,2) * t129) + t131 * t62 + t157 * t154 + t148 * t162 + t139 * t164 + t167 * t161 + t1 * (-mrSges(7,2) * t173 + mrSges(7,3) * t128) + t2 * (mrSges(7,1) * t173 - mrSges(7,3) * t129) + t42 * (mrSges(6,1) * t173 + mrSges(6,2) * t174); (t11 * t145 - t12 * t144) * mrSges(7,3) + t336 * qJD(4) + (-pkin(3) * t61 - qJ(4) * t54 - t103 * t166 + t389 * t84 - t83 * t99) * m(5) + m(7) * (t1 * t151 + t11 * t92 + t12 * t91 + t150 * t2) + (-t19 + t92) * t73 + (-t10 * mrSges(6,3) + t42 * mrSges(6,2) + t347 / 0.2e1 - t345 / 0.2e1 + t8 * t267 + t266 * t376 + t264 * t375 + t262 * t369 + t15 * t362 + t14 * t363 + t56 / 0.2e1 + (-t1 * t220 - t2 * t224) * mrSges(7,3) + (m(6) * t10 - m(7) * t8 + t359) * t228 + (mrSges(7,3) * t257 + t224 * t377 + t24 * t268 + t261 * t367 + t263 * t373 + t265 * t371 + t363 * t39) * qJD(6)) * t225 + (-t9 * mrSges(6,3) + t379 - t55 / 0.2e1 + t42 * mrSges(6,1) - t346 / 0.2e1 + t44 / 0.2e1 + t43 / 0.2e1 + (m(6) * t9 + t95) * t228 + (t380 + Ifges(7,3) / 0.2e1) * t112 + t274) * t221 + t196 + t197 + (((-pkin(3) * qJD(3) - t83) * mrSges(5,1) + qJD(3) * (Ifges(6,5) * t225 - Ifges(6,6) * t221) / 0.2e1 + t306 * Ifges(5,4) - t233 + Ifges(5,6) * t399 + t384) * t226 + (-t339 / 0.2e1 + (Ifges(4,4) / 0.2e1 + Ifges(5,6) / 0.2e1) * t294 + t306 * Ifges(4,6) + (-qJ(4) * qJD(3) - t84) * mrSges(5,1) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1 - Ifges(4,1) / 0.2e1 - Ifges(5,2) / 0.2e1) * t293 + t234 * t221 + t231 - t236) * t222) * t316 + t322 * t98 + t323 * t99 + (Ifges(7,5) * t145 + Ifges(7,6) * t144) * t367 + (Ifges(7,1) * t145 + Ifges(7,4) * t144) * t371 + (Ifges(7,4) * t145 + Ifges(7,2) * t144) * t373 + t144 * t377 + t240 + (-t20 + t91) * t72 + m(6) * (qJ(4) * t42 + qJD(4) * t74) - t40 * t57 + qJ(4) * t62 - t86 * t96 - t52 * t114 - t51 * t115 - t145 * t39 / 0.2e1 - t24 * (-mrSges(7,1) * t144 + mrSges(7,2) * t145) + t150 * t32 + t151 * t33 - t166 * t164 - m(6) * (t34 * t51 + t35 * t52 + t74 * t86) - m(7) * (t11 * t19 + t12 * t20 + t24 * t40) + (t231 + (t262 * t367 + t264 * t373 + t266 * t371 + t234 - t400) * t221 + (-m(6) * t254 + m(7) * t338 + t225 * t114 - t221 * t337) * t228) * qJD(5); -t142 * t73 - t143 * t72 - t336 * t209 + (mrSges(5,1) * t310 + t164 * t222) * t316 + (t114 * t294 + (-t220 * t73 + t224 * t72 + t114) * qJD(5) + t359) * t225 + (t255 * qJD(6) - t199 * t337 + t256 + t95) * t221 + ((-qJD(5) * t257 - t8) * t225 + (qJD(5) * t24 - qJD(6) * t258 + t270) * t221 - t11 * t142 - t12 * t143 + t294 * t338) * m(7) + (t10 * t225 - t199 * t254 - t209 * t74 + t221 * t9) * m(6) + (t103 * t294 + t209 * t84 + t61) * m(5); t270 * mrSges(7,3) + t337 * t35 - t8 * t268 - pkin(5) * t16 + t385 + ((t380 - Ifges(6,1) / 0.2e1) * t153 + t386 - t397) * t152 + t195 + t232 * t153 + t220 * t378 + t261 * t369 + t263 * t375 + t265 * t376 + t14 * t362 - t18 * t72 - t17 * t73 - t34 * t114 + t397 * qJD(6) + (m(7) * t270 + t256 + (-m(7) * t258 + t255) * qJD(6)) * pkin(11) + (-pkin(5) * t8 - t11 * t17 - t12 * t18 - t24 * t35) * m(7); -t24 * (mrSges(7,1) * t106 + mrSges(7,2) * t105) + (Ifges(7,1) * t105 - t356) * t371 + t38 * t370 + (Ifges(7,5) * t105 - Ifges(7,6) * t106) * t367 - t11 * t72 + t12 * t73 + (t105 * t11 + t106 * t12) * mrSges(7,3) + t274 + t13 + (-Ifges(7,2) * t106 + t104 + t39) * t373;];
tauc  = t3(:);
