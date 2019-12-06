% Calculate vector of inverse dynamics joint torques for
% S5RRPRR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR6_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR6_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:35:42
% EndTime: 2019-12-05 18:35:53
% DurationCPUTime: 5.13s
% Computational Cost: add. (5392->428), mult. (8267->587), div. (0->0), fcn. (5206->14), ass. (0->244)
t230 = sin(pkin(9));
t231 = cos(pkin(9));
t191 = -pkin(3) * t231 - pkin(7) * t230 - pkin(2);
t233 = sin(qJ(4));
t237 = cos(qJ(4));
t317 = t231 * t237;
t142 = qJ(3) * t317 + t233 * t191;
t234 = sin(qJ(2));
t308 = qJD(3) * t231;
t238 = cos(qJ(2));
t316 = t231 * t238;
t343 = qJD(1) * pkin(1);
t374 = -qJD(4) * t142 - t233 * t308 - (-t233 * t316 + t234 * t237) * t343;
t305 = qJD(4) * t237;
t394 = (t233 * t234 + t237 * t316) * t343 - t191 * t305 - t237 * t308;
t296 = t238 * t343;
t271 = qJD(3) - t296;
t393 = mrSges(3,2) - mrSges(4,3);
t225 = t230 ^ 2;
t226 = t231 ^ 2;
t387 = mrSges(4,3) * (t225 + t226);
t306 = qJD(4) * t233;
t285 = t230 * t306;
t200 = pkin(8) * t285;
t392 = t200 + t374;
t318 = t231 * t233;
t291 = qJ(3) * t318;
t320 = t230 * t237;
t302 = pkin(8) * t320;
t391 = -(-t291 - t302) * qJD(4) + t394;
t351 = mrSges(4,2) * t230;
t390 = -m(5) * t191 + mrSges(3,1) - m(6) * (-(pkin(4) * t237 + pkin(3)) * t231 - pkin(2)) - t351 + (mrSges(5,3) + mrSges(6,3)) * t230;
t232 = sin(qJ(5));
t236 = cos(qJ(5));
t184 = t232 * t237 + t233 * t236;
t162 = t184 * t230;
t227 = qJD(1) + qJD(2);
t116 = t227 * t162;
t224 = qJDD(1) + qJDD(2);
t143 = (t224 * t237 - t227 * t306) * t230;
t144 = (-t224 * t233 - t227 * t305) * t230;
t45 = -qJD(5) * t116 + t143 * t236 + t144 * t232;
t389 = Ifges(6,5) * t45;
t262 = t232 * t233 - t236 * t237;
t257 = t262 * qJD(5);
t323 = t227 * t230;
t46 = -t143 * t232 + t144 * t236 + t257 * t323;
t388 = Ifges(6,6) * t46;
t386 = Ifges(5,5) * t143;
t385 = Ifges(5,6) * t144;
t324 = t224 * t231;
t195 = qJDD(4) - t324;
t384 = Ifges(5,3) * t195;
t189 = qJDD(5) + t195;
t383 = Ifges(6,3) * t189;
t382 = t271 * t230;
t322 = t227 * t231;
t196 = qJD(4) - t322;
t297 = t234 * t343;
t187 = qJ(3) * t227 + t297;
t333 = t187 * t225;
t381 = (mrSges(5,1) * t237 - mrSges(5,2) * t233) * t333 + t196 * t230 * (-Ifges(5,5) * t233 - Ifges(5,6) * t237) / 0.2e1;
t347 = Ifges(5,4) * t237;
t348 = Ifges(5,4) * t233;
t380 = (-t233 * (-Ifges(5,2) * t237 - t348) / 0.2e1 + t237 * (-Ifges(5,1) * t233 - t347) / 0.2e1) * t225;
t292 = qJD(2) * t343;
t336 = qJDD(1) * pkin(1);
t179 = t234 * t336 + t238 * t292;
t310 = qJD(3) * t227;
t128 = qJ(3) * t224 + t179 + t310;
t113 = t191 * t227 + t271;
t68 = t113 * t233 + t187 * t317;
t178 = -t234 * t292 + t238 * t336;
t260 = qJDD(3) - t178;
t87 = t191 * t224 + t260;
t27 = -qJD(4) * t68 - t128 * t318 + t237 * t87;
t16 = pkin(4) * t195 - pkin(8) * t143 + t27;
t283 = t231 * t306;
t26 = t113 * t305 + t128 * t317 - t187 * t283 + t233 * t87;
t17 = pkin(8) * t144 + t26;
t321 = t230 * t233;
t289 = t227 * t321;
t59 = -pkin(8) * t289 + t68;
t341 = t232 * t59;
t288 = t227 * t320;
t67 = t237 * t113 - t187 * t318;
t58 = -pkin(8) * t288 + t67;
t47 = pkin(4) * t196 + t58;
t18 = t236 * t47 - t341;
t4 = qJD(5) * t18 + t16 * t232 + t17 * t236;
t340 = t236 * t59;
t19 = t232 * t47 + t340;
t5 = -qJD(5) * t19 + t16 * t236 - t17 * t232;
t379 = t5 * mrSges(6,1) - t4 * mrSges(6,2);
t378 = -t27 * mrSges(5,1) + t26 * mrSges(5,2);
t303 = pkin(8) * t321;
t110 = -t303 + t142;
t177 = t237 * t191;
t91 = -t302 + t177 + (-qJ(3) * t233 - pkin(4)) * t231;
t48 = -t110 * t232 + t236 * t91;
t377 = qJD(5) * t48 + t392 * t232 - t391 * t236;
t49 = t110 * t236 + t232 * t91;
t376 = -qJD(5) * t49 + t391 * t232 + t392 * t236;
t361 = m(6) * pkin(4);
t375 = mrSges(5,1) + t361;
t373 = -qJ(3) * t283 - t394;
t365 = -qJD(4) - qJD(5);
t105 = t365 * t184;
t372 = t184 * t322 + t105;
t371 = (t322 + t365) * t262;
t284 = t230 * t305;
t201 = pkin(4) * t284;
t370 = t201 + t382;
t355 = pkin(1) * t238;
t173 = t191 - t355;
t357 = pkin(1) * t234;
t214 = qJ(3) + t357;
t101 = t233 * t173 + t214 * t317;
t229 = qJ(1) + qJ(2);
t220 = sin(t229);
t319 = t230 * (-pkin(8) - pkin(7));
t222 = cos(t229);
t326 = t222 * t233;
t369 = pkin(4) * t326 + t220 * t319;
t311 = qJD(3) * t187;
t335 = t128 * t225;
t366 = qJ(3) * t335 + t225 * t311 - t296 * t333;
t228 = qJ(4) + qJ(5);
t219 = sin(t228);
t221 = cos(t228);
t327 = t222 * t231;
t133 = t219 * t327 - t220 * t221;
t134 = -t219 * t220 - t221 * t327;
t157 = -t220 * t237 + t222 * t318;
t328 = t220 * t233;
t158 = -t222 * t317 - t328;
t364 = mrSges(4,1) * t327 - t158 * mrSges(5,1) - t134 * mrSges(6,1) - t157 * mrSges(5,2) - t133 * mrSges(6,2) - t393 * t220 + t390 * t222;
t329 = t220 * t231;
t131 = t219 * t329 + t221 * t222;
t132 = -t219 * t222 + t221 * t329;
t155 = t220 * t318 + t222 * t237;
t156 = t220 * t317 - t326;
t363 = mrSges(4,1) * t329 + t156 * mrSges(5,1) + t132 * mrSges(6,1) - t155 * mrSges(5,2) - t131 * mrSges(6,2) + t390 * t220 + t393 * t222;
t362 = t227 ^ 2;
t118 = -t232 * t289 + t236 * t288;
t359 = t118 / 0.2e1;
t235 = sin(qJ(1));
t356 = pkin(1) * t235;
t239 = cos(qJ(1));
t354 = pkin(1) * t239;
t353 = pkin(4) * t233;
t352 = g(1) * t230;
t350 = mrSges(4,3) * t226;
t349 = mrSges(6,3) * t116;
t346 = Ifges(6,4) * t118;
t345 = pkin(1) * qJD(2);
t342 = t118 * mrSges(6,3);
t294 = t238 * t345;
t205 = qJD(3) + t294;
t339 = t205 * t333 + t214 * t335;
t338 = qJ(3) * t220;
t337 = qJ(3) * t230;
t334 = t128 * t230;
t332 = t187 * t227;
t331 = t205 * t230;
t330 = t214 * t230;
t325 = t224 * t230;
t314 = t131 * mrSges(6,1) + t132 * mrSges(6,2);
t313 = -t133 * mrSges(6,1) + t134 * mrSges(6,2);
t301 = t383 + t388 + t389;
t299 = mrSges(4,3) * t224 * t225;
t298 = t224 * t350;
t295 = t234 * t345;
t290 = t214 * t318;
t287 = t384 + t385 + t386;
t286 = t173 * t305 + t205 * t317 + t233 * t295;
t282 = t320 / 0.2e1;
t211 = t222 * qJ(3);
t280 = t211 - t356;
t279 = -pkin(2) * t220 + t211;
t277 = t227 * mrSges(3,1) * t357;
t276 = pkin(4) * t288;
t275 = mrSges(5,3) * t289;
t274 = mrSges(5,3) * t288;
t273 = t227 * t296;
t166 = -mrSges(4,1) * t324 + mrSges(4,2) * t325;
t268 = -t231 * mrSges(4,1) + t351;
t267 = mrSges(5,1) * t233 + mrSges(5,2) * t237;
t266 = -mrSges(6,1) * t219 - mrSges(6,2) * t221;
t265 = -pkin(2) * t222 - t338;
t165 = t237 * t173;
t75 = -t302 + t165 + (-t214 * t233 - pkin(4)) * t231;
t86 = -t303 + t101;
t39 = -t232 * t86 + t236 * t75;
t40 = t232 * t75 + t236 * t86;
t138 = -mrSges(5,2) * t196 - t275;
t139 = mrSges(5,1) * t196 - t274;
t263 = t138 * t237 - t139 * t233;
t259 = (t237 * Ifges(5,1) - t348) * t230;
t258 = (-t233 * Ifges(5,2) + t347) * t230;
t256 = -pkin(4) * t328 + t222 * t319 - t338;
t63 = -qJD(4) * t101 - t205 * t318 + t237 * t295;
t114 = Ifges(6,4) * t116;
t115 = (t227 * t353 + t187) * t230;
t190 = qJD(5) + t196;
t51 = -Ifges(6,2) * t116 + Ifges(6,6) * t190 + t346;
t52 = Ifges(6,1) * t118 + Ifges(6,5) * t190 - t114;
t244 = -t115 * (mrSges(6,1) * t118 - mrSges(6,2) * t116) - t18 * t349 + t51 * t359 + t301 - t118 * (-Ifges(6,1) * t116 - t346) / 0.2e1 - t190 * (-Ifges(6,5) * t116 - Ifges(6,6) * t118) / 0.2e1 + (-Ifges(6,2) * t118 - t114 + t52) * t116 / 0.2e1 + t379;
t150 = -pkin(2) * t224 + t260;
t170 = t267 * t230;
t66 = -pkin(4) * t144 + t334;
t78 = t105 * t230;
t79 = (qJD(4) * t262 + t257) * t230;
t94 = Ifges(5,6) * t196 + t227 * t258;
t95 = Ifges(5,5) * t196 + t227 * t259;
t241 = -(t301 + t287) * t231 / 0.2e1 + (Ifges(5,1) * t143 + Ifges(5,4) * t144 + Ifges(5,5) * t195) * t282 + t143 * t259 / 0.2e1 + (-t26 * t321 - t27 * t320 - t68 * t284 + t67 * t285) * mrSges(5,3) + (mrSges(6,1) * t66 - mrSges(6,3) * t4 - Ifges(6,4) * t45 - Ifges(6,2) * t46 - Ifges(6,6) * t189) * t162 + (Ifges(6,1) * t78 + Ifges(6,4) * t79) * t359 - (Ifges(5,4) * t143 + Ifges(5,2) * t144 + Ifges(5,6) * t195) * t321 / 0.2e1 + t115 * (-mrSges(6,1) * t79 + mrSges(6,2) * t78) - t116 * (Ifges(6,4) * t78 + Ifges(6,2) * t79) / 0.2e1 + t78 * t52 / 0.2e1 + t79 * t51 / 0.2e1 + t178 * mrSges(3,1) - t179 * mrSges(3,2) + t190 * (Ifges(6,5) * t78 + Ifges(6,6) * t79) / 0.2e1 + Ifges(3,3) * t224 + t144 * t258 / 0.2e1 + t170 * t334 + mrSges(4,3) * t335 + t128 * t350 + t150 * t268 + (-t18 * t78 + t19 * t79) * mrSges(6,3) + (-t383 / 0.2e1 - t388 / 0.2e1 - t389 / 0.2e1 - t384 / 0.2e1 + Ifges(4,2) * t324 + Ifges(4,4) * t325 - t385 / 0.2e1 - t386 / 0.2e1 + t378 - t379) * t231 + ((-mrSges(6,2) * t66 + mrSges(6,3) * t5 - Ifges(6,1) * t45 - Ifges(6,4) * t46 - Ifges(6,5) * t189) * t262 + t195 * (Ifges(5,5) * t237 - Ifges(5,6) * t233) / 0.2e1 + Ifges(4,4) * t324 + Ifges(4,1) * t325) * t230 + (-(t233 * t95 + t237 * t94) * t230 / 0.2e1 + t381 + t380 * t227) * qJD(4);
t216 = -pkin(2) - t355;
t208 = pkin(4) * t321;
t182 = t208 + t337;
t181 = -pkin(2) * t227 + t271;
t168 = t268 * t227;
t167 = t208 + t330;
t161 = t225 * t332;
t149 = t201 + t331;
t141 = t177 - t291;
t140 = t267 * t323;
t100 = t165 - t290;
t89 = -mrSges(5,2) * t195 + mrSges(5,3) * t144;
t88 = mrSges(5,1) * t195 - mrSges(5,3) * t143;
t85 = mrSges(6,1) * t190 - t342;
t84 = -mrSges(6,2) * t190 - t349;
t69 = -mrSges(5,1) * t144 + mrSges(5,2) * t143;
t65 = mrSges(6,1) * t116 + mrSges(6,2) * t118;
t62 = -t214 * t283 + t286;
t56 = t200 + t63;
t55 = (-t290 - t302) * qJD(4) + t286;
t38 = -mrSges(6,2) * t189 + mrSges(6,3) * t46;
t37 = mrSges(6,1) * t189 - mrSges(6,3) * t45;
t23 = t236 * t58 - t341;
t22 = -t232 * t58 - t340;
t15 = -mrSges(6,1) * t46 + mrSges(6,2) * t45;
t7 = -qJD(5) * t40 - t232 * t55 + t236 * t56;
t6 = qJD(5) * t39 + t232 * t56 + t236 * t55;
t1 = [t39 * t37 + t40 * t38 + t224 * mrSges(3,1) * t355 + t241 + t214 * t298 + t214 * t299 + t168 * t295 + (-m(6) * (t280 + t369) - m(5) * t280 + m(3) * t356 - m(4) * (t279 - t356) + t235 * mrSges(2,1) + t239 * mrSges(2,2) + t363) * g(3) + (m(3) * t354 - m(5) * (-t338 - t354) - m(6) * (t256 - t354) - m(4) * (t265 - t354) + t239 * mrSges(2,1) - t235 * mrSges(2,2) + t364) * g(2) + (-t224 * t357 - t227 * t294) * mrSges(3,2) + m(6) * (t115 * t149 + t167 * t66 + t18 * t7 + t19 * t6 + t39 * t5 + t4 * t40) - qJD(2) * t277 + t149 * t65 + t62 * t138 + t63 * t139 + t100 * t88 + t101 * t89 + t6 * t84 + t7 * t85 + t167 * t15 + t216 * t166 + m(4) * (t181 * t295 + t150 * t216 + (t128 * t214 + t187 * t205) * t226 + t339) + m(5) * (t100 * t27 + t101 * t26 + t62 * t68 + t63 * t67 + t339) + t205 * t227 * t387 + t69 * t330 + t140 * t331 + m(3) * (t178 * t238 + t179 * t234) * pkin(1) + Ifges(2,3) * qJDD(1); t382 * t140 + t48 * t37 + t49 * t38 - t168 * t297 + t241 + mrSges(3,2) * t273 + qJD(1) * t277 + t373 * t138 + t374 * t139 + (t141 * t27 + t142 * t26 + t373 * t68 + t374 * t67 + t366) * m(5) + (-m(4) * t279 - m(5) * t211 - m(6) * (t211 + t369) + t363) * g(3) + t370 * t65 + (-pkin(2) * t150 + (qJ(3) * t128 + t311) * t226 - (t187 * t226 * t238 + t181 * t234) * t343 + t366) * m(4) + (-m(4) * t265 + m(5) * t338 - m(6) * t256 + t364) * g(2) + (t298 + t299) * qJ(3) - pkin(2) * t166 + t141 * t88 + t142 * t89 + t376 * t85 + (t370 * t115 + t376 * t18 + t182 * t66 + t377 * t19 + t4 * t49 + t48 * t5) * m(6) + t377 * t84 + t182 * t15 + t69 * t337 + (-t273 + t310) * t387; -t262 * t37 + t184 * t38 + t233 * t89 + t237 * t88 + t372 * t85 + t371 * t84 + t263 * qJD(4) - t362 * t387 + (-t263 * t231 + (-t140 - t65) * t230) * t227 + t166 + (g(2) * t222 + g(3) * t220) * (-m(4) - m(5) - m(6)) + (-t115 * t323 + t18 * t372 + t184 * t4 + t371 * t19 - t262 * t5) * m(6) + (t233 * t26 + t237 * t27 - t161 + t196 * (-t233 * t67 + t237 * t68)) * m(5) + (-t226 * t332 + t150 - t161) * m(4); t95 * t289 / 0.2e1 + t244 - (-m(6) * t353 + t266) * t352 + (t232 * t4 + t236 * t5 + (-t18 * t232 + t19 * t236) * qJD(5)) * t361 - t23 * t84 - t22 * t85 + g(1) * t170 + t19 * t342 + t287 - t65 * t276 - m(6) * (t115 * t276 + t18 * t22 + t19 * t23) + (t274 + t139) * t68 + (-t138 - t275) * t67 - t380 * t362 + (-mrSges(5,2) * t158 + t157 * t375 - t313) * g(3) + (-mrSges(5,2) * t156 - t155 * t375 - t314) * g(2) + (t94 * t282 - t381) * t227 + ((-t232 * t85 + t236 * t84) * qJD(5) + t232 * t38 + t236 * t37) * pkin(4) - t378; t244 - g(2) * t314 - t18 * t84 - t266 * t352 + (t85 + t342) * t19 - g(3) * t313;];
tau = t1;
