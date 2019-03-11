% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-03-09 07:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRR9_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR9_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR9_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR9_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR9_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:23:04
% EndTime: 2019-03-09 07:23:30
% DurationCPUTime: 14.38s
% Computational Cost: add. (13212->660), mult. (29040->931), div. (0->0), fcn. (19382->8), ass. (0->304)
t260 = sin(qJ(3));
t264 = cos(qJ(3));
t291 = pkin(3) * t264 + pkin(8) * t260;
t233 = t291 * qJD(1);
t265 = -pkin(1) - pkin(7);
t247 = qJD(1) * t265 + qJD(2);
t263 = cos(qJ(4));
t259 = sin(qJ(4));
t327 = t259 * t264;
t169 = t263 * t233 - t247 * t327;
t393 = -pkin(9) - pkin(8);
t303 = qJD(4) * t393;
t308 = pkin(9) * t260 * t263;
t464 = -(pkin(4) * t264 + t308) * qJD(1) - t169 + t263 * t303;
t325 = t263 * t264;
t170 = t259 * t233 + t247 * t325;
t256 = t260 * qJD(1);
t302 = t259 * t256;
t463 = pkin(9) * t302 - t259 * t303 + t170;
t330 = t247 * t264;
t215 = -qJD(3) * pkin(3) - t330;
t319 = qJD(3) * t263;
t321 = qJD(1) * t264;
t224 = -t259 * t321 + t319;
t173 = -pkin(4) * t224 + t215;
t225 = qJD(3) * t259 + t263 * t321;
t258 = sin(qJ(5));
t262 = cos(qJ(5));
t292 = t262 * t224 - t225 * t258;
t104 = -pkin(5) * t292 + t173;
t261 = cos(qJ(6));
t257 = sin(qJ(6));
t429 = pkin(10) * t292;
t237 = pkin(3) * t260 - pkin(8) * t264 + qJ(2);
t208 = t237 * qJD(1);
t236 = t260 * t247;
t214 = qJD(3) * pkin(8) + t236;
t148 = t263 * t208 - t214 * t259;
t121 = -pkin(9) * t225 + t148;
t252 = t256 + qJD(4);
t110 = pkin(4) * t252 + t121;
t149 = t208 * t259 + t214 * t263;
t122 = pkin(9) * t224 + t149;
t117 = t262 * t122;
t54 = t110 * t258 + t117;
t48 = t54 + t429;
t341 = t257 * t48;
t245 = qJD(5) + t252;
t162 = t224 * t258 + t225 * t262;
t445 = pkin(10) * t162;
t115 = t258 * t122;
t53 = t262 * t110 - t115;
t47 = t53 - t445;
t43 = pkin(5) * t245 + t47;
t14 = t261 * t43 - t341;
t339 = t261 * t48;
t15 = t257 * t43 + t339;
t433 = -t162 * t257 + t261 * t292;
t309 = qJD(3) * qJD(4);
t314 = qJD(4) * t264;
t178 = t263 * t309 + (-t259 * t314 - t260 * t319) * qJD(1);
t320 = qJD(3) * t260;
t435 = -t259 * t320 + t263 * t314;
t179 = -qJD(1) * t435 - t259 * t309;
t71 = qJD(5) * t292 + t178 * t262 + t179 * t258;
t72 = -qJD(5) * t162 - t178 * t258 + t179 * t262;
t28 = qJD(6) * t433 + t257 * t72 + t261 * t71;
t97 = t162 * t261 + t257 * t292;
t29 = -qJD(6) * t97 - t257 * t71 + t261 * t72;
t318 = qJD(3) * t264;
t295 = qJD(1) * t318;
t307 = Ifges(7,5) * t28 + Ifges(7,6) * t29 + Ifges(7,3) * t295;
t366 = Ifges(7,4) * t97;
t239 = qJD(6) + t245;
t377 = -t239 / 0.2e1;
t395 = t97 / 0.2e1;
t396 = -t97 / 0.2e1;
t398 = -t433 / 0.2e1;
t312 = qJD(5) * t262;
t313 = qJD(5) * t258;
t222 = qJD(3) * t291 + qJD(2);
t193 = t222 * qJD(1);
t271 = -qJD(4) * t149 + t263 * t193;
t61 = -pkin(9) * t178 + (pkin(4) * qJD(1) - t247 * t259) * t318 + t271;
t301 = t247 * t318;
t315 = qJD(4) * t263;
t316 = qJD(4) * t259;
t86 = t259 * t193 + t208 * t315 - t214 * t316 + t263 * t301;
t63 = pkin(9) * t179 + t86;
t12 = t110 * t312 - t122 * t313 + t258 * t61 + t262 * t63;
t10 = pkin(10) * t72 + t12;
t13 = -qJD(5) * t54 - t258 * t63 + t262 * t61;
t9 = pkin(5) * t295 - pkin(10) * t71 + t13;
t2 = qJD(6) * t14 + t10 * t261 + t257 * t9;
t3 = -qJD(6) * t15 - t10 * t257 + t261 * t9;
t428 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t45 = Ifges(7,2) * t433 + Ifges(7,6) * t239 + t366;
t90 = Ifges(7,4) * t433;
t46 = Ifges(7,1) * t97 + Ifges(7,5) * t239 + t90;
t462 = (Ifges(7,1) * t433 - t366) * t396 + (Ifges(7,5) * t433 - Ifges(7,6) * t97) * t377 + (t14 * t433 + t15 * t97) * mrSges(7,3) - t104 * (mrSges(7,1) * t97 + mrSges(7,2) * t433) + t45 * t395 + t307 + t428 + (-Ifges(7,2) * t97 + t46 + t90) * t398;
t228 = t258 * t263 + t259 * t262;
t408 = qJD(4) + qJD(5);
t168 = t408 * t228;
t206 = t228 * qJD(1);
t187 = t260 * t206;
t454 = t168 + t187;
t277 = t258 * t259 - t262 * t263;
t414 = t260 * t277;
t188 = qJD(1) * t414;
t438 = t408 * t277;
t461 = t438 + t188;
t155 = Ifges(6,4) * t292;
t306 = Ifges(6,5) * t71 + Ifges(6,6) * t72 + Ifges(6,3) * t295;
t354 = Ifges(6,4) * t162;
t375 = -t245 / 0.2e1;
t387 = -t162 / 0.2e1;
t389 = -t292 / 0.2e1;
t419 = t13 * mrSges(6,1) - t12 * mrSges(6,2);
t85 = Ifges(6,1) * t162 + Ifges(6,5) * t245 + t155;
t460 = t306 + t419 + (Ifges(6,5) * t292 - Ifges(6,6) * t162) * t375 + (t162 * t54 + t292 * t53) * mrSges(6,3) + (-Ifges(6,2) * t162 + t155 + t85) * t389 - t173 * (mrSges(6,1) * t162 + mrSges(6,2) * t292) + (Ifges(6,1) * t292 - t354) * t387 + t462;
t242 = t393 * t259;
t243 = t393 * t263;
t442 = t242 * t312 + t243 * t313 + t258 * t464 - t463 * t262;
t177 = t258 * t242 - t262 * t243;
t441 = -qJD(5) * t177 + t463 * t258 + t262 * t464;
t457 = -pkin(5) * t321 + pkin(10) * t461 + t441;
t456 = pkin(10) * t454 - t442;
t195 = t228 * t264;
t411 = t277 * qJD(1) - qJD(3) * t195 + t260 * t438;
t197 = t277 * t264;
t410 = -qJD(3) * t197 - t168 * t260 - t206;
t346 = t225 * Ifges(5,4);
t142 = t224 * Ifges(5,2) + t252 * Ifges(5,6) + t346;
t217 = Ifges(5,4) * t224;
t143 = t225 * Ifges(5,1) + t252 * Ifges(5,5) + t217;
t281 = t148 * t263 + t149 * t259;
t355 = Ifges(5,4) * t263;
t286 = -Ifges(5,2) * t259 + t355;
t356 = Ifges(5,4) * t259;
t288 = Ifges(5,1) * t263 - t356;
t289 = mrSges(5,1) * t259 + mrSges(5,2) * t263;
t352 = Ifges(5,6) * t259;
t353 = Ifges(5,5) * t263;
t369 = t263 / 0.2e1;
t371 = -t259 / 0.2e1;
t378 = t225 / 0.2e1;
t448 = -t281 * mrSges(5,3) + t142 * t371 + t143 * t369 + t215 * t289 + (-t352 + t353) * t252 / 0.2e1 + t286 * t224 / 0.2e1 + t288 * t378;
t176 = t262 * t242 + t243 * t258;
t137 = -pkin(10) * t228 + t176;
t138 = -pkin(10) * t277 + t177;
t81 = t137 * t257 + t138 * t261;
t444 = -qJD(6) * t81 + t257 * t456 + t261 * t457;
t80 = t137 * t261 - t138 * t257;
t443 = qJD(6) * t80 + t257 * t457 - t261 * t456;
t254 = pkin(4) * t262 + pkin(5);
t310 = qJD(6) * t261;
t311 = qJD(6) * t257;
t329 = t257 * t258;
t59 = -t121 * t258 - t117;
t50 = t59 - t429;
t60 = t262 * t121 - t115;
t51 = t60 - t445;
t440 = -t257 * t50 - t261 * t51 + t254 * t310 + (-t258 * t311 + (t261 * t262 - t329) * qJD(5)) * pkin(4);
t328 = t258 * t261;
t439 = t257 * t51 - t261 * t50 - t254 * t311 + (-t258 * t310 + (-t257 * t262 - t328) * qJD(5)) * pkin(4);
t189 = -pkin(4) * t302 + t236;
t437 = pkin(4) * t316 + pkin(5) * t454 - t189;
t335 = Ifges(4,5) * qJD(3);
t358 = Ifges(4,4) * t260;
t436 = t335 / 0.2e1 + (t264 * Ifges(4,1) - t358) * qJD(1) / 0.2e1 + t448;
t84 = Ifges(6,2) * t292 + Ifges(6,6) * t245 + t354;
t430 = t84 / 0.2e1;
t296 = -Ifges(4,6) * qJD(3) / 0.2e1;
t194 = t228 * t260;
t133 = -t194 * t257 - t261 * t414;
t418 = -qJD(6) * t133 - t257 * t410 + t261 * t411;
t131 = -t194 * t261 + t257 * t414;
t417 = qJD(6) * t131 + t257 * t411 + t261 * t410;
t416 = qJ(2) * (m(3) + m(4));
t101 = -mrSges(6,1) * t292 + mrSges(6,2) * t162;
t368 = m(6) * t173;
t413 = t101 + t368;
t221 = t263 * t237;
t294 = -t259 * t265 + pkin(4);
t158 = -pkin(9) * t325 + t260 * t294 + t221;
t326 = t260 * t265;
t244 = t263 * t326;
t184 = t259 * t237 + t244;
t172 = -pkin(9) * t327 + t184;
t100 = t258 * t158 + t262 * t172;
t412 = Ifges(5,5) * t178 + Ifges(5,6) * t179;
t87 = -t259 * t301 + t271;
t282 = -t259 * t87 + t263 * t86;
t409 = t87 * mrSges(5,1) - t86 * mrSges(5,2);
t181 = -mrSges(5,2) * t252 + mrSges(5,3) * t224;
t182 = mrSges(5,1) * t252 - mrSges(5,3) * t225;
t278 = -t259 * t181 - t263 * t182;
t407 = -m(5) * t281 + t278;
t405 = m(5) / 0.2e1;
t404 = t28 / 0.2e1;
t403 = t29 / 0.2e1;
t400 = t71 / 0.2e1;
t399 = t72 / 0.2e1;
t397 = t433 / 0.2e1;
t132 = -t195 * t261 + t197 * t257;
t391 = t132 / 0.2e1;
t134 = -t195 * t257 - t197 * t261;
t390 = t134 / 0.2e1;
t388 = t292 / 0.2e1;
t386 = t162 / 0.2e1;
t385 = t178 / 0.2e1;
t384 = t179 / 0.2e1;
t383 = -t195 / 0.2e1;
t382 = -t197 / 0.2e1;
t381 = -t224 / 0.2e1;
t379 = -t225 / 0.2e1;
t376 = t239 / 0.2e1;
t374 = t245 / 0.2e1;
t373 = -t252 / 0.2e1;
t367 = m(7) * t104;
t359 = mrSges(4,1) * t260;
t357 = Ifges(4,4) * t264;
t351 = qJ(2) * mrSges(4,1);
t350 = qJ(2) * mrSges(4,2);
t124 = t187 * t261 - t188 * t257;
t166 = t228 * t261 - t257 * t277;
t65 = -qJD(6) * t166 - t168 * t261 + t257 * t438;
t337 = t124 - t65;
t125 = t187 * t257 + t188 * t261;
t165 = -t228 * t257 - t261 * t277;
t64 = qJD(6) * t165 - t168 * t257 - t261 * t438;
t336 = t125 - t64;
t333 = qJD(3) * mrSges(4,2);
t322 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t224 - mrSges(5,2) * t225 - mrSges(4,3) * t321;
t317 = qJD(3) * t265;
t305 = t259 * t326;
t304 = Ifges(5,3) * t295 + t412;
t255 = -pkin(4) * t263 - pkin(3);
t226 = t247 * t320;
t299 = t264 * t317;
t297 = -t335 / 0.2e1;
t139 = -pkin(4) * t179 + t226;
t99 = t262 * t158 - t172 * t258;
t223 = pkin(4) * t327 - t264 * t265;
t290 = mrSges(5,1) * t263 - mrSges(5,2) * t259;
t287 = Ifges(5,1) * t259 + t355;
t285 = Ifges(5,2) * t263 + t356;
t283 = Ifges(5,5) * t259 + Ifges(5,6) * t263;
t68 = pkin(5) * t260 + pkin(10) * t197 + t99;
t75 = -pkin(10) * t195 + t100;
t36 = -t257 * t75 + t261 * t68;
t37 = t257 * t68 + t261 * t75;
t280 = t148 * t259 - t149 * t263;
t156 = mrSges(5,1) * t295 - mrSges(5,3) * t178;
t157 = -mrSges(5,2) * t295 + mrSges(5,3) * t179;
t279 = -t259 * t156 + t263 * t157;
t180 = pkin(4) * t435 + t260 * t317;
t119 = -qJD(4) * t305 + t259 * t222 + t237 * t315 + t263 * t299;
t102 = -pkin(9) * t435 + t119;
t201 = t263 * t222;
t91 = t201 + (-t244 + (pkin(9) * t264 - t237) * t259) * qJD(4) + (t264 * t294 + t308) * qJD(3);
t32 = t262 * t102 + t158 * t312 - t172 * t313 + t258 * t91;
t33 = -qJD(5) * t100 - t102 * t258 + t262 * t91;
t267 = t14 * mrSges(7,1) + t148 * mrSges(5,1) + t53 * mrSges(6,1) + t252 * Ifges(5,3) + t225 * Ifges(5,5) + t224 * Ifges(5,6) + t296 - (-t260 * Ifges(4,2) + t357) * qJD(1) / 0.2e1 + t239 * Ifges(7,3) + t97 * Ifges(7,5) + t433 * Ifges(7,6) + t245 * Ifges(6,3) + t162 * Ifges(6,5) + t292 * Ifges(6,6) - t149 * mrSges(5,2) - t15 * mrSges(7,2) - t54 * mrSges(6,2);
t240 = -mrSges(4,3) * t256 - t333;
t232 = (t264 * mrSges(4,2) + t359) * qJD(1);
t199 = pkin(4) * t328 + t254 * t257;
t198 = -pkin(4) * t329 + t254 * t261;
t190 = pkin(5) * t277 + t255;
t183 = t221 - t305;
t164 = pkin(5) * t195 + t223;
t129 = mrSges(6,1) * t245 - mrSges(6,3) * t162;
t128 = -mrSges(6,2) * t245 + mrSges(6,3) * t292;
t123 = pkin(4) * t225 + pkin(5) * t162;
t120 = -qJD(4) * t184 - t259 * t299 + t201;
t118 = -mrSges(5,1) * t179 + mrSges(5,2) * t178;
t112 = t178 * Ifges(5,1) + t179 * Ifges(5,4) + Ifges(5,5) * t295;
t111 = t178 * Ifges(5,4) + t179 * Ifges(5,2) + Ifges(5,6) * t295;
t109 = t228 * t320 + t264 * t438;
t107 = qJD(3) * t414 - t168 * t264;
t82 = -pkin(5) * t109 + t180;
t79 = mrSges(7,1) * t239 - mrSges(7,3) * t97;
t78 = -mrSges(7,2) * t239 + mrSges(7,3) * t433;
t67 = -mrSges(6,2) * t295 + mrSges(6,3) * t72;
t66 = mrSges(6,1) * t295 - mrSges(6,3) * t71;
t52 = -pkin(5) * t72 + t139;
t49 = -mrSges(7,1) * t433 + mrSges(7,2) * t97;
t42 = -qJD(6) * t134 - t107 * t257 + t109 * t261;
t40 = qJD(6) * t132 + t107 * t261 + t109 * t257;
t38 = -mrSges(6,1) * t72 + mrSges(6,2) * t71;
t35 = t71 * Ifges(6,1) + t72 * Ifges(6,4) + Ifges(6,5) * t295;
t34 = t71 * Ifges(6,4) + t72 * Ifges(6,2) + Ifges(6,6) * t295;
t25 = -mrSges(7,2) * t295 + mrSges(7,3) * t29;
t24 = mrSges(7,1) * t295 - mrSges(7,3) * t28;
t21 = pkin(10) * t109 + t32;
t20 = pkin(5) * t318 - pkin(10) * t107 + t33;
t17 = t261 * t47 - t341;
t16 = -t257 * t47 - t339;
t8 = -mrSges(7,1) * t29 + mrSges(7,2) * t28;
t7 = t28 * Ifges(7,1) + t29 * Ifges(7,4) + Ifges(7,5) * t295;
t6 = t28 * Ifges(7,4) + t29 * Ifges(7,2) + Ifges(7,6) * t295;
t5 = -qJD(6) * t37 + t20 * t261 - t21 * t257;
t4 = qJD(6) * t36 + t20 * t257 + t21 * t261;
t1 = [(t132 * t2 - t134 * t3 - t14 * t40 + t15 * t42) * mrSges(7,3) + m(5) * (t119 * t149 + t120 * t148 + t183 * t87 + t184 * t86) + (t297 + (-0.2e1 * t350 + 0.3e1 / 0.2e1 * t358) * qJD(1) + (m(5) * t215 - t322) * t265 - t436) * t320 + (t307 + t306 + t304 + t412) * t260 / 0.2e1 + (t111 * t371 + t112 * t369 - t265 * t118 + t286 * t384 + t288 * t385 + qJD(1) * qJD(2) * mrSges(4,2) + (-t259 * t86 - t263 * t87) * mrSges(5,3) + (t283 * t373 + t285 * t381 + t287 * t379 + t215 * t290 + t143 * t371 - t263 * t142 / 0.2e1 + t280 * mrSges(5,3)) * qJD(4) + ((0.2e1 * t351 + Ifges(6,5) * t382 + Ifges(6,6) * t383 + Ifges(7,5) * t390 + Ifges(7,6) * t391 + (-0.3e1 / 0.2e1 * Ifges(4,4) + t353 / 0.2e1 - t352 / 0.2e1) * t264) * qJD(1) + t265 * t240 + t267 + t296 + (-m(5) * t265 + t289) * t236) * qJD(3)) * t264 + m(7) * (t104 * t82 + t14 * t5 + t15 * t4 + t164 * t52 + t2 * t37 + t3 * t36) + m(6) * (t100 * t12 + t13 * t99 + t139 * t223 + t173 * t180 + t32 * t54 + t33 * t53) + (Ifges(7,1) * t134 + Ifges(7,4) * t132) * t404 + (Ifges(7,4) * t134 + Ifges(7,2) * t132) * t403 + t223 * t38 + t180 * t101 + t119 * t181 + t120 * t182 + t183 * t156 + t184 * t157 + t173 * (-mrSges(6,1) * t109 + mrSges(6,2) * t107) + t164 * t8 + t52 * (-mrSges(7,1) * t132 + mrSges(7,2) * t134) + t32 * t128 + t33 * t129 + t104 * (-mrSges(7,1) * t42 + mrSges(7,2) * t40) + t107 * t85 / 0.2e1 + t99 * t66 + t100 * t67 + t4 * t78 + t5 * t79 + t82 * t49 + t40 * t46 / 0.2e1 + t42 * t45 / 0.2e1 + t37 * t25 + t36 * t24 + (Ifges(7,4) * t40 + Ifges(7,2) * t42) * t397 + (Ifges(6,1) * t107 + Ifges(6,4) * t109) * t386 + (Ifges(6,4) * t107 + Ifges(6,2) * t109) * t388 + t7 * t390 + t6 * t391 + (Ifges(7,1) * t40 + Ifges(7,4) * t42) * t395 + (Ifges(6,5) * t107 + Ifges(6,6) * t109) * t374 + (Ifges(7,5) * t40 + Ifges(7,6) * t42) * t376 + t35 * t382 + t34 * t383 + (t232 + ((2 * mrSges(3,3)) + t359 + 0.2e1 * t416) * qJD(1)) * qJD(2) + (-Ifges(6,4) * t197 - Ifges(6,2) * t195) * t399 + (-t107 * t53 + t109 * t54 - t12 * t195 + t13 * t197) * mrSges(6,3) + t139 * (mrSges(6,1) * t195 - mrSges(6,2) * t197) + (-Ifges(6,1) * t197 - Ifges(6,4) * t195) * t400 + (Ifges(6,6) * t399 + Ifges(6,5) * t400 + Ifges(7,6) * t403 + Ifges(7,5) * t404 + (0.3e1 / 0.2e1 * Ifges(4,2) + Ifges(5,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(4,1) + Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1) * t295 + t409 + t419 + t428) * t260 + t109 * t430; t131 * t24 + t133 * t25 - t194 * t66 - t414 * t67 + t418 * t79 + t417 * t78 + t411 * t129 + t410 * t128 + (-t118 - t38 - t8 + (t181 * t263 - t182 * t259 + t240) * qJD(3)) * t264 + (t278 * qJD(4) + (t101 + t49 - t322) * qJD(3) + t279) * t260 - m(5) * t280 * t318 + 0.2e1 * ((-qJD(4) * t281 + t282) * t405 + (t367 / 0.2e1 + t368 / 0.2e1 + (t215 - t330) * t405) * qJD(3)) * t260 + (t131 * t3 + t133 * t2 + t14 * t418 + t15 * t417 - t264 * t52) * m(7) + (-t12 * t414 - t13 * t194 - t139 * t264 + t410 * t54 + t411 * t53) * m(6) + (-t232 + (-mrSges(3,3) - t416) * qJD(1) + t407) * qJD(1); (-t125 / 0.2e1 + t64 / 0.2e1) * t46 + (-t188 / 0.2e1 - t438 / 0.2e1) * t85 + (-Ifges(6,1) * t438 - Ifges(6,4) * t168) * t386 + (-t124 / 0.2e1 + t65 / 0.2e1) * t45 - t277 * t34 / 0.2e1 + t139 * (mrSges(6,1) * t277 + mrSges(6,2) * t228) + (Ifges(6,4) * t228 - Ifges(6,2) * t277) * t399 + (Ifges(6,1) * t228 - Ifges(6,4) * t277) * t400 + t437 * t49 + ((t297 + (t350 - t358 / 0.2e1) * qJD(1) + t436) * t260 + ((-t351 + t357 / 0.2e1 + (Ifges(4,1) / 0.2e1 - Ifges(4,2) / 0.2e1) * t260) * qJD(1) - t267 + t296) * t264 + (Ifges(6,5) * t228 + Ifges(7,5) * t166 - Ifges(6,6) * t277 + Ifges(7,6) * t165 + t283) * t318 / 0.2e1) * qJD(1) + (pkin(4) * t259 * t413 + pkin(8) * t407 + t448) * qJD(4) + (-pkin(3) * t226 + pkin(8) * t282 - t148 * t169 - t149 * t170 - t215 * t236) * m(5) + (-t187 / 0.2e1 - t168 / 0.2e1) * t84 + t441 * t129 + (t12 * t177 + t13 * t176 + t139 * t255 - t173 * t189 + t441 * t53 + t442 * t54) * m(6) + t442 * t128 + t282 * mrSges(5,3) + t279 * pkin(8) + t443 * t78 + t444 * t79 + (t104 * t437 + t14 * t444 + t15 * t443 + t190 * t52 + t2 * t81 + t3 * t80) * m(7) + (-t12 * t277 - t13 * t228 - t454 * t54 + t461 * t53) * mrSges(6,3) + (mrSges(6,1) * t454 - mrSges(6,2) * t461) * t173 + t259 * t112 / 0.2e1 + t255 * t38 + t228 * t35 / 0.2e1 - t189 * t101 + t190 * t8 - t170 * t181 - t169 * t182 + t177 * t67 + t176 * t66 + t166 * t7 / 0.2e1 + t165 * t6 / 0.2e1 + t52 * (-mrSges(7,1) * t165 + mrSges(7,2) * t166) - pkin(3) * t118 + t80 * t24 + t81 * t25 + (Ifges(7,4) * t64 + Ifges(7,2) * t65) * t397 + (Ifges(7,4) * t125 + Ifges(7,2) * t124) * t398 + (Ifges(7,4) * t166 + Ifges(7,2) * t165) * t403 + (Ifges(7,1) * t166 + Ifges(7,4) * t165) * t404 + t285 * t384 + t287 * t385 + (Ifges(6,1) * t188 + Ifges(6,4) * t187) * t387 + (Ifges(6,4) * t188 + Ifges(6,2) * t187) * t389 + (Ifges(7,1) * t64 + Ifges(7,4) * t65) * t395 + (Ifges(7,1) * t125 + Ifges(7,4) * t124) * t396 + (Ifges(6,5) * t188 + Ifges(6,6) * t187) * t375 + (Ifges(7,5) * t64 + Ifges(7,6) * t65) * t376 + (Ifges(7,5) * t125 + Ifges(7,6) * t124) * t377 + t111 * t369 + (-Ifges(6,4) * t438 - Ifges(6,2) * t168) * t388 + (-Ifges(6,5) * t438 - Ifges(6,6) * t168) * t374 + ((-t240 - t333) * t264 + ((-mrSges(4,1) - t290) * qJD(3) + t322) * t260) * t247 + (mrSges(7,1) * t337 - mrSges(7,2) * t336) * t104 + (t14 * t336 - t15 * t337 + t165 * t2 - t166 * t3) * mrSges(7,3); t162 * t430 + (t258 * t67 + t262 * t66 + (t128 * t262 - t129 * t258) * qJD(5) + m(6) * (t12 * t258 + t13 * t262 + t312 * t54 - t313 * t53) - t413 * t225) * pkin(4) - m(6) * (t53 * t59 + t54 * t60) + t439 * t79 + t440 * t78 + (-t104 * t123 + t14 * t439 + t15 * t440 + t198 * t3 + t199 * t2) * m(7) + (-Ifges(5,2) * t225 + t143 + t217) * t381 + (t148 * t224 + t149 * t225) * mrSges(5,3) + t304 - t215 * (mrSges(5,1) * t225 + mrSges(5,2) * t224) + t198 * t24 + t199 * t25 - t148 * t181 + t149 * t182 - t60 * t128 - t59 * t129 - t123 * t49 + (Ifges(5,5) * t224 - Ifges(5,6) * t225) * t373 + t142 * t378 + (Ifges(5,1) * t224 - t346) * t379 + t409 + t460; -m(7) * (t14 * t16 + t15 * t17) - t53 * t128 + t54 * t129 - t17 * t78 - t16 * t79 + (t261 * t24 + t257 * t25 + (-t257 * t79 + t261 * t78) * qJD(6) + m(7) * (-t14 * t311 + t15 * t310 + t2 * t257 + t261 * t3) + (-t49 - t367) * t162) * pkin(5) + t84 * t386 + t460; -t14 * t78 + t15 * t79 + t462;];
tauc  = t1(:);
