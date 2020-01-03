% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRR8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR8_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR8_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:24:23
% EndTime: 2019-12-31 22:24:40
% DurationCPUTime: 7.27s
% Computational Cost: add. (9696->530), mult. (23614->734), div. (0->0), fcn. (16635->8), ass. (0->266)
t253 = sin(qJ(4));
t345 = -t253 / 0.2e1;
t257 = cos(qJ(4));
t254 = sin(qJ(3));
t245 = pkin(2) * t254 + pkin(8);
t340 = -pkin(9) - t245;
t285 = qJD(4) * t340;
t258 = cos(qJ(3));
t333 = pkin(2) * qJD(3);
t293 = t258 * t333;
t259 = cos(qJ(2));
t299 = qJD(1) * t259;
t255 = sin(qJ(2));
t300 = qJD(1) * t255;
t214 = -t254 * t300 + t258 * t299;
t308 = t214 * t253;
t295 = pkin(9) * t308;
t226 = t254 * t259 + t258 * t255;
t215 = t226 * qJD(1);
t181 = pkin(3) * t215 - pkin(8) * t214;
t155 = pkin(2) * t300 + t181;
t365 = -pkin(7) - pkin(6);
t242 = t365 * t259;
t229 = qJD(1) * t242;
t216 = t254 * t229;
t240 = t365 * t255;
t228 = qJD(1) * t240;
t187 = t228 * t258 + t216;
t98 = t253 * t155 + t257 * t187;
t406 = t253 * t285 + t257 * t293 + t295 - t98;
t307 = t214 * t257;
t281 = t215 * pkin(4) - pkin(9) * t307;
t97 = t257 * t155 - t187 * t253;
t405 = -t253 * t293 + t257 * t285 - t281 - t97;
t251 = qJD(2) + qJD(3);
t404 = t251 * Ifges(4,6) / 0.2e1;
t219 = qJD(2) * pkin(2) + t228;
t184 = t219 * t258 + t216;
t102 = t253 * t181 + t257 * t184;
t364 = -pkin(9) - pkin(8);
t290 = qJD(4) * t364;
t403 = t253 * t290 - t102 + t295;
t101 = t257 * t181 - t184 * t253;
t402 = t257 * t290 - t101 - t281;
t196 = -t215 * t253 + t251 * t257;
t211 = qJD(4) - t214;
t197 = t215 * t257 + t251 * t253;
t337 = Ifges(5,4) * t197;
t104 = Ifges(5,2) * t196 + Ifges(5,6) * t211 + t337;
t159 = -pkin(3) * t251 - t184;
t279 = mrSges(5,1) * t253 + mrSges(5,2) * t257;
t401 = t104 * t345 + t159 * t279;
t248 = -pkin(2) * t259 - pkin(1);
t238 = qJD(1) * t248;
t320 = t251 * Ifges(4,5);
t400 = t238 * mrSges(4,2) + t320 / 0.2e1;
t321 = t215 * Ifges(4,4);
t399 = t404 + t321 / 0.2e1 + t214 * Ifges(4,2) / 0.2e1;
t117 = -pkin(4) * t196 + t159;
t252 = sin(qJ(5));
t256 = cos(qJ(5));
t127 = t196 * t252 + t197 * t256;
t150 = -t214 * pkin(3) - t215 * pkin(8) + t238;
t217 = t258 * t229;
t185 = t219 * t254 - t217;
t160 = pkin(8) * t251 + t185;
t82 = t150 * t253 + t160 * t257;
t66 = pkin(9) * t196 + t82;
t318 = t252 * t66;
t81 = t257 * t150 - t160 * t253;
t65 = -pkin(9) * t197 + t81;
t59 = pkin(4) * t211 + t65;
t16 = t256 * t59 - t318;
t316 = t256 * t66;
t17 = t252 * t59 + t316;
t192 = t251 * t226;
t176 = t192 * qJD(1);
t171 = Ifges(6,3) * t176;
t283 = t256 * t196 - t197 * t252;
t334 = Ifges(6,4) * t127;
t207 = qJD(5) + t211;
t350 = -t207 / 0.2e1;
t356 = -t127 / 0.2e1;
t398 = t171 + (Ifges(6,5) * t283 - Ifges(6,6) * t127) * t350 + (t127 * t17 + t16 * t283) * mrSges(6,3) - t117 * (mrSges(6,1) * t127 + mrSges(6,2) * t283) + (Ifges(6,1) * t283 - t334) * t356;
t221 = t340 * t253;
t250 = t257 * pkin(9);
t222 = t245 * t257 + t250;
t178 = t221 * t252 + t222 * t256;
t397 = -qJD(5) * t178 - t406 * t252 + t405 * t256;
t177 = t221 * t256 - t222 * t252;
t396 = qJD(5) * t177 + t405 * t252 + t406 * t256;
t291 = qJD(2) * t365;
t282 = qJD(1) * t291;
t220 = t255 * t282;
t269 = t259 * t282;
t108 = qJD(3) * t185 + t220 * t254 - t258 * t269;
t224 = t254 * t255 - t258 * t259;
t191 = t251 * t224;
t175 = t191 * qJD(1);
t114 = qJD(4) * t196 - t175 * t257;
t115 = -qJD(4) * t197 + t175 * t253;
t58 = -mrSges(5,1) * t115 + mrSges(5,2) * t114;
t394 = m(5) * t108 + t58;
t225 = t252 * t257 + t253 * t256;
t148 = t225 * t214;
t376 = qJD(4) + qJD(5);
t190 = t376 * t225;
t393 = t148 - t190;
t270 = t252 * t253 - t256 * t257;
t149 = t270 * t214;
t189 = t376 * t270;
t392 = t149 - t189;
t186 = t228 * t254 - t217;
t206 = pkin(4) * t308;
t297 = qJD(4) * t253;
t292 = pkin(4) * t297;
t391 = t254 * t333 - t186 - t206 + t292;
t107 = qJD(3) * t184 + t258 * t220 + t254 * t269;
t296 = qJD(4) * t257;
t298 = qJD(2) * t255;
t294 = pkin(2) * t298;
t90 = pkin(3) * t176 + pkin(8) * t175 + qJD(1) * t294;
t20 = t257 * t107 + t150 * t296 - t160 * t297 + t253 * t90;
t21 = -qJD(4) * t82 - t107 * t253 + t257 * t90;
t390 = t20 * t257 - t21 * t253;
t12 = pkin(4) * t176 - pkin(9) * t114 + t21;
t13 = pkin(9) * t115 + t20;
t3 = qJD(5) * t16 + t12 * t252 + t13 * t256;
t34 = qJD(5) * t283 + t114 * t256 + t115 * t252;
t35 = -qJD(5) * t127 - t114 * t252 + t115 * t256;
t4 = -qJD(5) * t17 + t12 * t256 - t13 * t252;
t389 = t4 * mrSges(6,1) - t3 * mrSges(6,2) + Ifges(6,5) * t34 + Ifges(6,6) * t35;
t123 = Ifges(6,4) * t283;
t388 = -Ifges(6,2) * t127 + t123;
t274 = Ifges(5,5) * t257 - Ifges(5,6) * t253;
t335 = Ifges(5,4) * t257;
t276 = -Ifges(5,2) * t253 + t335;
t336 = Ifges(5,4) * t253;
t278 = Ifges(5,1) * t257 - t336;
t195 = Ifges(5,4) * t196;
t105 = Ifges(5,1) * t197 + Ifges(5,5) * t211 + t195;
t302 = t257 * t105;
t351 = t197 / 0.2e1;
t387 = t211 * t274 / 0.2e1 + t278 * t351 + t196 * t276 / 0.2e1 + t302 / 0.2e1 + t401;
t386 = -t238 * mrSges(4,1) - t81 * mrSges(5,1) - t16 * mrSges(6,1) + t82 * mrSges(5,2) + t17 * mrSges(6,2) + t399;
t371 = t34 / 0.2e1;
t370 = t35 / 0.2e1;
t354 = t176 / 0.2e1;
t239 = t364 * t253;
t241 = pkin(8) * t257 + t250;
t200 = t239 * t252 + t241 * t256;
t380 = -qJD(5) * t200 - t403 * t252 + t402 * t256;
t198 = t239 * t256 - t241 * t252;
t379 = qJD(5) * t198 + t402 * t252 + t403 * t256;
t165 = t270 * t226;
t183 = t224 * pkin(3) - t226 * pkin(8) + t248;
t201 = t240 * t254 - t242 * t258;
t194 = t257 * t201;
t119 = t253 * t183 + t194;
t378 = t258 * t240 + t242 * t254;
t377 = -t253 * t81 + t257 * t82;
t374 = t21 * mrSges(5,1) - t20 * mrSges(5,2) + Ifges(5,5) * t114 + Ifges(5,6) * t115 + t389;
t373 = Ifges(6,4) * t371 + Ifges(6,2) * t370 + Ifges(6,6) * t354;
t372 = Ifges(6,1) * t371 + Ifges(6,4) * t370 + Ifges(6,5) * t354;
t55 = Ifges(6,2) * t283 + Ifges(6,6) * t207 + t334;
t369 = -t55 / 0.2e1;
t368 = t55 / 0.2e1;
t56 = Ifges(6,1) * t127 + Ifges(6,5) * t207 + t123;
t367 = -t56 / 0.2e1;
t366 = t56 / 0.2e1;
t363 = pkin(1) * mrSges(3,1);
t362 = pkin(1) * mrSges(3,2);
t360 = t114 / 0.2e1;
t359 = t115 / 0.2e1;
t358 = -t283 / 0.2e1;
t357 = t283 / 0.2e1;
t355 = t127 / 0.2e1;
t353 = -t196 / 0.2e1;
t352 = -t197 / 0.2e1;
t349 = t207 / 0.2e1;
t348 = -t211 / 0.2e1;
t344 = t257 / 0.2e1;
t343 = m(4) * t238;
t341 = pkin(2) * t258;
t339 = mrSges(4,3) * t214;
t338 = Ifges(3,4) * t255;
t209 = Ifges(4,4) * t214;
t332 = t283 * Ifges(6,6);
t331 = t127 * Ifges(6,5);
t330 = t196 * Ifges(5,6);
t329 = t197 * Ifges(5,5);
t327 = t207 * Ifges(6,3);
t325 = t211 * Ifges(5,3);
t323 = t215 * mrSges(4,3);
t322 = t215 * Ifges(4,1);
t314 = Ifges(3,5) * qJD(2);
t313 = Ifges(3,6) * qJD(2);
t312 = qJD(2) * mrSges(3,1);
t311 = qJD(2) * mrSges(3,2);
t310 = t108 * t378;
t305 = t226 * t253;
t301 = mrSges(4,1) * t251 + mrSges(5,1) * t196 - mrSges(5,2) * t197 - t323;
t247 = -pkin(4) * t257 - pkin(3);
t289 = t314 / 0.2e1;
t288 = -t313 / 0.2e1;
t113 = pkin(3) * t192 + pkin(8) * t191 + t294;
t233 = t255 * t291;
t234 = t259 * t291;
t130 = qJD(3) * t378 + t233 * t258 + t234 * t254;
t284 = t257 * t113 - t130 * t253;
t118 = t257 * t183 - t201 * t253;
t280 = mrSges(5,1) * t257 - mrSges(5,2) * t253;
t277 = Ifges(5,1) * t253 + t335;
t275 = Ifges(5,2) * t257 + t336;
t273 = Ifges(5,5) * t253 + Ifges(5,6) * t257;
t73 = pkin(4) * t224 - t226 * t250 + t118;
t88 = -pkin(9) * t305 + t119;
t40 = -t252 * t88 + t256 * t73;
t41 = t252 * t73 + t256 * t88;
t70 = mrSges(5,1) * t176 - mrSges(5,3) * t114;
t71 = -mrSges(5,2) * t176 + mrSges(5,3) * t115;
t272 = -t253 * t70 + t257 * t71;
t271 = -t253 * t82 - t257 * t81;
t268 = t271 * mrSges(5,3);
t267 = -t191 * t253 + t226 * t296;
t36 = t253 * t113 + t257 * t130 + t183 * t296 - t201 * t297;
t131 = qJD(3) * t201 + t233 * t254 - t258 * t234;
t261 = m(5) * (qJD(4) * t271 + t390);
t103 = t325 + t329 + t330;
t157 = t209 + t320 + t322;
t44 = t114 * Ifges(5,4) + t115 * Ifges(5,2) + t176 * Ifges(5,6);
t45 = t114 * Ifges(5,1) + t115 * Ifges(5,4) + t176 * Ifges(5,5);
t54 = t327 + t331 + t332;
t57 = -pkin(4) * t115 + t108;
t260 = t387 * qJD(4) + (t81 * t307 + t82 * t308 + t390) * mrSges(5,3) + (-t392 * t16 + t393 * t17 - t4 * t225 - t3 * t270) * mrSges(6,3) + (-t393 * mrSges(6,1) + t392 * mrSges(6,2)) * t117 - (Ifges(4,1) * t214 + t103 - t321 + t54) * t215 / 0.2e1 - (-Ifges(4,2) * t215 + t157 + t209 + t302) * t214 / 0.2e1 + t57 * (mrSges(6,1) * t270 + mrSges(6,2) * t225) + (Ifges(6,4) * t225 - Ifges(6,2) * t270) * t370 + (Ifges(6,1) * t225 - Ifges(6,4) * t270) * t371 + (Ifges(6,5) * t225 - Ifges(6,6) * t270 + t273) * t354 - t270 * t373 + (-Ifges(6,5) * t189 - Ifges(6,6) * t190) * t349 + (-Ifges(6,1) * t189 - Ifges(6,4) * t190) * t355 + (-Ifges(6,4) * t189 - Ifges(6,2) * t190) * t357 + (-Ifges(6,1) * t149 - Ifges(6,4) * t148) * t356 + t184 * t339 + (Ifges(5,5) * t352 + Ifges(6,5) * t356 + Ifges(5,6) * t353 + Ifges(6,6) * t358 + Ifges(5,3) * t348 + Ifges(6,3) * t350 + t386 + t404) * t215 + (-t280 - mrSges(4,1)) * t108 + (t274 * t348 + t276 * t353 + t278 * t352 - t400 - t401) * t214 + (-Ifges(6,4) * t149 - Ifges(6,2) * t148) * t358 + t253 * t45 / 0.2e1 + t44 * t344 + t275 * t359 + t277 * t360 - t189 * t366 - t149 * t367 - t190 * t368 - t148 * t369 + t225 * t372 + (-Ifges(6,5) * t149 - Ifges(6,6) * t148) * t350 - t107 * mrSges(4,2) - Ifges(4,5) * t175 - Ifges(4,6) * t176;
t249 = Ifges(3,4) * t299;
t237 = mrSges(3,3) * t299 - t311;
t236 = -mrSges(3,3) * t300 + t312;
t235 = t247 - t341;
t213 = Ifges(3,1) * t300 + t249 + t314;
t212 = t313 + (Ifges(3,2) * t259 + t338) * qJD(1);
t204 = -mrSges(4,2) * t251 + t339;
t180 = -mrSges(4,1) * t214 + mrSges(4,2) * t215;
t172 = Ifges(5,3) * t176;
t164 = t225 * t226;
t151 = pkin(4) * t305 - t378;
t143 = mrSges(5,1) * t211 - mrSges(5,3) * t197;
t142 = -mrSges(5,2) * t211 + mrSges(5,3) * t196;
t135 = t185 + t206;
t92 = mrSges(6,1) * t207 - mrSges(6,3) * t127;
t91 = -mrSges(6,2) * t207 + mrSges(6,3) * t283;
t69 = pkin(4) * t267 + t131;
t64 = -mrSges(6,1) * t283 + mrSges(6,2) * t127;
t53 = t165 * t376 + t225 * t191;
t52 = -t190 * t226 + t191 * t270;
t37 = -qJD(4) * t119 + t284;
t25 = t256 * t65 - t318;
t24 = -t252 * t65 - t316;
t23 = -mrSges(6,2) * t176 + mrSges(6,3) * t35;
t22 = mrSges(6,1) * t176 - mrSges(6,3) * t34;
t19 = -pkin(9) * t267 + t36;
t15 = t191 * t250 + pkin(4) * t192 + (-t194 + (pkin(9) * t226 - t183) * t253) * qJD(4) + t284;
t11 = -mrSges(6,1) * t35 + mrSges(6,2) * t34;
t6 = -qJD(5) * t41 + t15 * t256 - t19 * t252;
t5 = qJD(5) * t40 + t15 * t252 + t19 * t256;
t1 = [(-t16 * t52 - t164 * t3 + t165 * t4 + t17 * t53) * mrSges(6,3) + (-Ifges(6,5) * t165 - Ifges(6,6) * t164) * t354 + (-Ifges(6,4) * t165 - Ifges(6,2) * t164) * t370 + (-Ifges(6,1) * t165 - Ifges(6,4) * t164) * t371 + t57 * (mrSges(6,1) * t164 - mrSges(6,2) * t165) + m(4) * (t107 * t201 + t130 * t185 - t131 * t184 - t310) + m(5) * (t118 * t21 + t119 * t20 + t131 * t159 + t36 * t82 + t37 * t81 - t310) - t378 * t58 + (t175 * t378 - t176 * t201 + t184 * t191 - t185 * t192) * mrSges(4,3) + (Ifges(4,4) * t175 + t171 / 0.2e1 + t172 / 0.2e1 - t107 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) + Ifges(6,3) / 0.2e1) * t176 + t374) * t224 + (t276 * t359 + t274 * t354 + t278 * t360 - Ifges(4,1) * t175 - Ifges(4,4) * t176 + t44 * t345 + t45 * t344 + (mrSges(4,3) + t279) * t108 + (-t20 * t253 - t21 * t257) * mrSges(5,3) + (-t257 * t104 / 0.2e1 + t105 * t345 + t273 * t348 + t159 * t280 + t275 * t353 + t277 * t352 - t377 * mrSges(5,3)) * qJD(4)) * t226 + t248 * (mrSges(4,1) * t176 - mrSges(4,2) * t175) + (-pkin(6) * t237 - t212 / 0.2e1 + t288 + (-0.2e1 * t363 - 0.3e1 / 0.2e1 * t338 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t259) * qJD(1) + (t180 + 0.2e1 * t343 + qJD(1) * (mrSges(4,1) * t224 + mrSges(4,2) * t226)) * pkin(2)) * t298 + m(6) * (t117 * t69 + t151 * t57 + t16 * t6 + t17 * t5 + t3 * t41 + t4 * t40) - t301 * t131 - (t209 / 0.2e1 + t322 / 0.2e1 + t157 / 0.2e1 + t268 + t387 + t400) * t191 + (-pkin(6) * t236 + t213 / 0.2e1 + t289 + (-0.2e1 * t362 + 0.3e1 / 0.2e1 * Ifges(3,4) * t259) * qJD(1)) * t259 * qJD(2) + (t325 / 0.2e1 + t330 / 0.2e1 + t329 / 0.2e1 + t327 / 0.2e1 + t54 / 0.2e1 + t103 / 0.2e1 + t332 / 0.2e1 + t331 / 0.2e1 - t386 - t399) * t192 + t130 * t204 + (Ifges(6,5) * t52 + Ifges(6,6) * t53) * t349 + (Ifges(6,1) * t52 + Ifges(6,4) * t53) * t355 + (Ifges(6,4) * t52 + Ifges(6,2) * t53) * t357 + t52 * t366 + t53 * t368 - t165 * t372 - t164 * t373 + t40 * t22 + t41 * t23 + t69 * t64 + t5 * t91 + t6 * t92 + t117 * (-mrSges(6,1) * t53 + mrSges(6,2) * t52) + t118 * t70 + t119 * t71 + t36 * t142 + t37 * t143 + t151 * t11; -m(4) * (-t184 * t186 + t185 * t187) + (m(4) * (t107 * t254 - t108 * t258) + (t258 * t175 - t254 * t176) * mrSges(4,3) + ((-m(4) * t184 + m(5) * t159 - t301) * t254 + (m(4) * t185 + m(5) * t377 + t257 * t142 - t253 * t143 + t204) * t258) * qJD(3)) * pkin(2) + ((t289 - t213 / 0.2e1 - t249 / 0.2e1 + qJD(1) * t362 + (t236 - t312) * pkin(6)) * t259 + (t288 + t212 / 0.2e1 + (t363 + t338 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t259) * qJD(1) + (t237 + t311) * pkin(6) + (-t180 - t343) * pkin(2)) * t255) * qJD(1) + t268 * qJD(4) - m(5) * (t159 * t186 + t81 * t97 + t82 * t98) + t397 * t92 + t301 * t186 + t396 * t91 + t391 * t64 + t260 + t235 * t11 - t187 * t204 + t185 * t323 - t98 * t142 - t97 * t143 + t177 * t22 + t178 * t23 + t394 * (-pkin(3) - t341) + (t261 + t272 + (-t142 * t253 - t143 * t257) * qJD(4)) * t245 + (t391 * t117 + t397 * t16 + t396 * t17 + t177 * t4 + t178 * t3 + t235 * t57) * m(6); -m(5) * (t101 * t81 + t102 * t82 + t159 * t185) + t380 * t92 + t379 * t91 + (t301 + t323) * t185 + t260 + t247 * t11 + t198 * t22 + t200 * t23 - t184 * t204 + pkin(8) * t261 + t272 * pkin(8) - t135 * t64 - t102 * t142 - t101 * t143 + ((-t81 * mrSges(5,3) - pkin(8) * t143) * t257 + (-t82 * mrSges(5,3) + pkin(4) * t64 - pkin(8) * t142) * t253) * qJD(4) - t394 * pkin(3) + (t198 * t4 + t200 * t3 + t247 * t57 + t379 * t17 + t380 * t16 + (-t135 + t292) * t117) * m(6); t283 * t367 - t127 * t369 + t374 - m(6) * (t16 * t24 + t17 * t25) + t172 + (-t197 * t64 + t256 * t22 + t252 * t23 + (-t252 * t92 + t256 * t91) * qJD(5) + (-t117 * t197 + t252 * t3 + t256 * t4 + (-t16 * t252 + t17 * t256) * qJD(5)) * m(6)) * pkin(4) + (-Ifges(5,2) * t197 + t105 + t195) * t353 + (t196 * t81 + t197 * t82) * mrSges(5,3) - t159 * (mrSges(5,1) * t197 + mrSges(5,2) * t196) + (Ifges(5,5) * t196 - Ifges(5,6) * t197) * t348 + t104 * t351 + (Ifges(5,1) * t196 - t337) * t352 + t388 * t358 - t25 * t91 - t24 * t92 - t81 * t142 + t82 * t143 + t398; t55 * t355 - t16 * t91 + t17 * t92 + (t388 + t56) * t358 + t389 + t398;];
tauc = t1(:);
