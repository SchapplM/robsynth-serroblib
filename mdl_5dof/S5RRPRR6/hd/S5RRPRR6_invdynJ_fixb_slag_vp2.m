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
% m [6x1]
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
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 11:16:51
% EndTime: 2022-01-20 11:17:03
% DurationCPUTime: 5.34s
% Computational Cost: add. (5392->433), mult. (8267->592), div. (0->0), fcn. (5206->14), ass. (0->246)
t239 = sin(pkin(9));
t240 = cos(pkin(9));
t192 = -pkin(3) * t240 - pkin(7) * t239 - pkin(2);
t242 = sin(qJ(4));
t246 = cos(qJ(4));
t325 = t240 * t246;
t142 = qJ(3) * t325 + t242 * t192;
t243 = sin(qJ(2));
t316 = qJD(3) * t240;
t247 = cos(qJ(2));
t324 = t240 * t247;
t354 = pkin(1) * qJD(1);
t381 = -qJD(4) * t142 - t242 * t316 - (-t242 * t324 + t243 * t246) * t354;
t313 = qJD(4) * t246;
t401 = (t242 * t243 + t246 * t324) * t354 - t192 * t313 - t246 * t316;
t303 = t247 * t354;
t275 = qJD(3) - t303;
t400 = -mrSges(4,3) + mrSges(3,2);
t399 = mrSges(5,3) + mrSges(6,3);
t234 = t239 ^ 2;
t235 = t240 ^ 2;
t394 = mrSges(4,3) * (t234 + t235);
t314 = qJD(4) * t242;
t290 = t239 * t314;
t205 = pkin(8) * t290;
t398 = t205 + t381;
t326 = t240 * t242;
t298 = qJ(3) * t326;
t328 = t239 * t246;
t309 = pkin(8) * t328;
t397 = -(-t298 - t309) * qJD(4) + t401;
t241 = sin(qJ(5));
t245 = cos(qJ(5));
t185 = t241 * t246 + t242 * t245;
t162 = t185 * t239;
t236 = qJD(1) + qJD(2);
t116 = t236 * t162;
t233 = qJDD(1) + qJDD(2);
t143 = (t233 * t246 - t236 * t314) * t239;
t144 = (-t233 * t242 - t236 * t313) * t239;
t45 = -qJD(5) * t116 + t143 * t245 + t144 * t241;
t396 = Ifges(6,5) * t45;
t267 = t241 * t242 - t245 * t246;
t261 = t267 * qJD(5);
t331 = t236 * t239;
t46 = -t143 * t241 + t144 * t245 + t261 * t331;
t395 = Ifges(6,6) * t46;
t393 = Ifges(5,5) * t143;
t392 = Ifges(5,6) * t144;
t332 = t233 * t240;
t195 = qJDD(4) - t332;
t391 = Ifges(5,3) * t195;
t190 = qJDD(5) + t195;
t390 = Ifges(6,3) * t190;
t389 = t275 * t239;
t330 = t236 * t240;
t196 = qJD(4) - t330;
t304 = t243 * t354;
t188 = qJ(3) * t236 + t304;
t342 = t188 * t234;
t388 = t196 * t239 * (-Ifges(5,5) * t242 - Ifges(5,6) * t246) / 0.2e1 + (mrSges(5,1) * t246 - mrSges(5,2) * t242) * t342;
t356 = Ifges(5,4) * t246;
t357 = Ifges(5,4) * t242;
t387 = (-t242 * (-Ifges(5,2) * t246 - t357) / 0.2e1 + t246 * (-Ifges(5,1) * t242 - t356) / 0.2e1) * t234;
t353 = pkin(1) * qJD(2);
t299 = qJD(1) * t353;
t346 = pkin(1) * qJDD(1);
t180 = t243 * t346 + t247 * t299;
t317 = qJD(3) * t236;
t128 = qJ(3) * t233 + t180 + t317;
t113 = t192 * t236 + t275;
t68 = t113 * t242 + t188 * t325;
t179 = -t243 * t299 + t247 * t346;
t264 = qJDD(3) - t179;
t87 = t192 * t233 + t264;
t27 = -qJD(4) * t68 - t128 * t326 + t246 * t87;
t16 = pkin(4) * t195 - pkin(8) * t143 + t27;
t288 = t240 * t314;
t26 = t113 * t313 + t128 * t325 - t188 * t288 + t242 * t87;
t17 = pkin(8) * t144 + t26;
t329 = t239 * t242;
t296 = t236 * t329;
t59 = -pkin(8) * t296 + t68;
t349 = t241 * t59;
t295 = t236 * t328;
t67 = t246 * t113 - t188 * t326;
t58 = -pkin(8) * t295 + t67;
t47 = pkin(4) * t196 + t58;
t18 = t245 * t47 - t349;
t4 = qJD(5) * t18 + t16 * t241 + t17 * t245;
t348 = t245 * t59;
t19 = t241 * t47 + t348;
t5 = -qJD(5) * t19 + t16 * t245 - t17 * t241;
t386 = t5 * mrSges(6,1) - t4 * mrSges(6,2);
t385 = -t27 * mrSges(5,1) + t26 * mrSges(5,2);
t310 = pkin(8) * t329;
t110 = -t310 + t142;
t178 = t246 * t192;
t91 = -t309 + t178 + (-qJ(3) * t242 - pkin(4)) * t240;
t49 = t110 * t245 + t241 * t91;
t384 = -qJD(5) * t49 + t397 * t241 + t398 * t245;
t48 = -t110 * t241 + t245 * t91;
t383 = qJD(5) * t48 + t398 * t241 - t397 * t245;
t368 = m(6) * pkin(4);
t382 = mrSges(5,1) + t368;
t380 = -qJ(3) * t288 - t401;
t372 = -qJD(4) - qJD(5);
t105 = t372 * t185;
t379 = t185 * t330 + t105;
t378 = (t330 + t372) * t267;
t289 = t239 * t313;
t206 = pkin(4) * t289;
t377 = t206 + t389;
t362 = pkin(1) * t247;
t173 = t192 - t362;
t364 = pkin(1) * t243;
t221 = qJ(3) + t364;
t101 = t242 * t173 + t221 * t325;
t238 = qJ(1) + qJ(2);
t228 = sin(t238);
t327 = t239 * (-pkin(8) - pkin(7));
t230 = cos(t238);
t334 = t230 * t242;
t376 = pkin(4) * t334 + t228 * t327;
t318 = qJD(3) * t188;
t344 = t128 * t234;
t373 = qJ(3) * t344 + t234 * t318 - t303 * t342;
t237 = qJ(4) + qJ(5);
t227 = sin(t237);
t229 = cos(t237);
t338 = t228 * t240;
t131 = t227 * t338 + t229 * t230;
t132 = t227 * t230 - t229 * t338;
t155 = t228 * t326 + t230 * t246;
t156 = -t228 * t325 + t334;
t223 = pkin(4) * t246 + pkin(3);
t350 = t239 * mrSges(4,2);
t371 = mrSges(4,1) * t338 - t156 * mrSges(5,1) - t132 * mrSges(6,1) - t155 * mrSges(5,2) - t131 * mrSges(6,2) + t400 * t230 + (-m(5) * t192 - m(6) * (-t223 * t240 - pkin(2)) - t350 + mrSges(3,1) + t399 * t239) * t228;
t335 = t230 * t240;
t133 = -t227 * t335 + t228 * t229;
t134 = t227 * t228 + t229 * t335;
t157 = t228 * t246 - t230 * t326;
t337 = t228 * t242;
t158 = t230 * t325 + t337;
t336 = t230 * t239;
t370 = -t230 * mrSges(3,1) - mrSges(4,1) * t335 - t158 * mrSges(5,1) - t134 * mrSges(6,1) - t157 * mrSges(5,2) - t133 * mrSges(6,2) + t400 * t228 + (mrSges(4,2) - t399) * t336;
t369 = t236 ^ 2;
t118 = -t241 * t296 + t245 * t295;
t366 = t118 / 0.2e1;
t244 = sin(qJ(1));
t363 = pkin(1) * t244;
t361 = pkin(4) * t242;
t360 = g(3) * t239;
t248 = cos(qJ(1));
t231 = t248 * pkin(1);
t359 = mrSges(4,3) * t235;
t358 = mrSges(6,3) * t116;
t355 = Ifges(6,4) * t118;
t351 = t118 * mrSges(6,3);
t301 = t247 * t353;
t210 = qJD(3) + t301;
t347 = t210 * t342 + t221 * t344;
t345 = qJ(3) * t239;
t343 = t128 * t239;
t341 = t188 * t236;
t340 = t210 * t239;
t339 = t221 * t239;
t333 = t233 * t239;
t322 = -t131 * mrSges(6,1) + t132 * mrSges(6,2);
t321 = t133 * mrSges(6,1) - t134 * mrSges(6,2);
t319 = t230 * pkin(2) + t228 * qJ(3);
t308 = t390 + t395 + t396;
t306 = mrSges(4,3) * t233 * t234;
t305 = t233 * t359;
t302 = t243 * t353;
t297 = t221 * t326;
t293 = t391 + t392 + t393;
t291 = t173 * t313 + t210 * t325 + t242 * t302;
t287 = t328 / 0.2e1;
t217 = t230 * qJ(3);
t285 = t217 - t363;
t284 = -pkin(2) * t228 + t217;
t282 = t236 * mrSges(3,1) * t364;
t281 = pkin(4) * t295;
t280 = mrSges(5,3) * t296;
t279 = mrSges(5,3) * t295;
t278 = t236 * t303;
t276 = pkin(3) * t335 + pkin(7) * t336 + t319;
t166 = -mrSges(4,1) * t332 + mrSges(4,2) * t333;
t272 = -t240 * mrSges(4,1) + t350;
t271 = mrSges(5,1) * t242 + mrSges(5,2) * t246;
t270 = -mrSges(6,1) * t227 - mrSges(6,2) * t229;
t165 = t246 * t173;
t75 = -t309 + t165 + (-t221 * t242 - pkin(4)) * t240;
t86 = -t310 + t101;
t39 = -t241 * t86 + t245 * t75;
t40 = t241 * t75 + t245 * t86;
t138 = -mrSges(5,2) * t196 - t280;
t139 = mrSges(5,1) * t196 - t279;
t268 = t138 * t246 - t139 * t242;
t263 = (t246 * Ifges(5,1) - t357) * t239;
t262 = (-t242 * Ifges(5,2) + t356) * t239;
t256 = pkin(4) * t337 + t223 * t335 - t230 * t327 + t319;
t63 = -qJD(4) * t101 - t210 * t326 + t246 * t302;
t114 = Ifges(6,4) * t116;
t115 = (t236 * t361 + t188) * t239;
t191 = qJD(5) + t196;
t51 = -Ifges(6,2) * t116 + Ifges(6,6) * t191 + t355;
t52 = Ifges(6,1) * t118 + Ifges(6,5) * t191 - t114;
t252 = -t115 * (mrSges(6,1) * t118 - mrSges(6,2) * t116) - t18 * t358 + t51 * t366 + t308 - t118 * (-Ifges(6,1) * t116 - t355) / 0.2e1 - t191 * (-Ifges(6,5) * t116 - Ifges(6,6) * t118) / 0.2e1 + (-Ifges(6,2) * t118 - t114 + t52) * t116 / 0.2e1 + t386;
t150 = -pkin(2) * t233 + t264;
t170 = t271 * t239;
t66 = -pkin(4) * t144 + t343;
t78 = t105 * t239;
t79 = (qJD(4) * t267 + t261) * t239;
t94 = Ifges(5,6) * t196 + t236 * t262;
t95 = Ifges(5,5) * t196 + t236 * t263;
t250 = Ifges(3,3) * t233 + t191 * (Ifges(6,5) * t78 + Ifges(6,6) * t79) / 0.2e1 - t180 * mrSges(3,2) + t179 * mrSges(3,1) + t115 * (-mrSges(6,1) * t79 + mrSges(6,2) * t78) - t116 * (Ifges(6,4) * t78 + Ifges(6,2) * t79) / 0.2e1 + t78 * t52 / 0.2e1 + t79 * t51 / 0.2e1 + t150 * t272 + (-t26 * t329 - t27 * t328 - t289 * t68 + t290 * t67) * mrSges(5,3) + (mrSges(6,1) * t66 - mrSges(6,3) * t4 - Ifges(6,4) * t45 - Ifges(6,2) * t46 - Ifges(6,6) * t190) * t162 + t144 * t262 / 0.2e1 + (-t18 * t78 + t19 * t79) * mrSges(6,3) + t128 * t359 + (Ifges(6,1) * t78 + Ifges(6,4) * t79) * t366 + t170 * t343 + mrSges(4,3) * t344 + (Ifges(5,1) * t143 + Ifges(5,4) * t144 + Ifges(5,5) * t195) * t287 + t143 * t263 / 0.2e1 + (-t393 / 0.2e1 - t392 / 0.2e1 - t391 / 0.2e1 - t390 / 0.2e1 - t395 / 0.2e1 - t396 / 0.2e1 + Ifges(4,2) * t332 + Ifges(4,4) * t333 + t385 - t386) * t240 - (t308 + t293) * t240 / 0.2e1 - (Ifges(5,4) * t143 + Ifges(5,2) * t144 + Ifges(5,6) * t195) * t329 / 0.2e1 + (t195 * (Ifges(5,5) * t246 - Ifges(5,6) * t242) / 0.2e1 + Ifges(4,4) * t332 + Ifges(4,1) * t333 + (-mrSges(6,2) * t66 + mrSges(6,3) * t5 - Ifges(6,1) * t45 - Ifges(6,4) * t46 - Ifges(6,5) * t190) * t267) * t239 + (t387 * t236 - (t242 * t95 + t246 * t94) * t239 / 0.2e1 + t388) * qJD(4);
t224 = -pkin(2) - t362;
t213 = pkin(4) * t329;
t183 = t213 + t345;
t182 = -pkin(2) * t236 + t275;
t168 = t272 * t236;
t167 = t213 + t339;
t161 = t234 * t341;
t149 = t206 + t340;
t141 = t178 - t298;
t140 = t271 * t331;
t100 = t165 - t297;
t89 = -mrSges(5,2) * t195 + mrSges(5,3) * t144;
t88 = mrSges(5,1) * t195 - mrSges(5,3) * t143;
t85 = mrSges(6,1) * t191 - t351;
t84 = -mrSges(6,2) * t191 - t358;
t69 = -mrSges(5,1) * t144 + mrSges(5,2) * t143;
t65 = mrSges(6,1) * t116 + mrSges(6,2) * t118;
t62 = -t221 * t288 + t291;
t56 = t205 + t63;
t55 = (-t297 - t309) * qJD(4) + t291;
t38 = -mrSges(6,2) * t190 + mrSges(6,3) * t46;
t37 = mrSges(6,1) * t190 - mrSges(6,3) * t45;
t23 = t245 * t58 - t349;
t22 = -t241 * t58 - t348;
t15 = -mrSges(6,1) * t46 + mrSges(6,2) * t45;
t7 = -qJD(5) * t40 - t241 * t55 + t245 * t56;
t6 = qJD(5) * t39 + t241 * t56 + t245 * t55;
t1 = [t224 * t166 + t167 * t15 + t149 * t65 + t62 * t138 + t63 * t139 + t100 * t88 + t101 * t89 + t6 * t84 + t7 * t85 + t39 * t37 + t40 * t38 + (t244 * mrSges(2,1) + t248 * mrSges(2,2) - m(6) * (t285 + t376) - m(5) * t285 + m(3) * t363 - m(4) * (t284 - t363) + t371) * g(1) + (-t248 * mrSges(2,1) + t244 * mrSges(2,2) - m(5) * (t231 + t276) - m(6) * (t231 + t256) - m(4) * (t231 + t319) - m(3) * t231 + t370) * g(2) + (-t233 * t364 - t236 * t301) * mrSges(3,2) + m(3) * (t179 * t247 + t180 * t243) * pkin(1) + t221 * t305 + t221 * t306 + t168 * t302 + t250 + t69 * t339 + t140 * t340 + t233 * mrSges(3,1) * t362 + t210 * t236 * t394 - qJD(2) * t282 + m(6) * (t115 * t149 + t167 * t66 + t18 * t7 + t19 * t6 + t39 * t5 + t4 * t40) + Ifges(2,3) * qJDD(1) + m(4) * (t182 * t302 + t150 * t224 + (t128 * t221 + t188 * t210) * t235 + t347) + m(5) * (t100 * t27 + t101 * t26 + t62 * t68 + t63 * t67 + t347); t183 * t15 - pkin(2) * t166 + t141 * t88 + t142 * t89 + t48 * t37 + t49 * t38 + t383 * t84 + t380 * t138 + (t141 * t27 + t142 * t26 + t380 * t68 + t381 * t67 + t373) * m(5) + t381 * t139 + (-m(6) * (t217 + t376) - m(5) * t217 - m(4) * t284 + t371) * g(1) + t377 * t65 + (-pkin(2) * t150 + (qJ(3) * t128 + t318) * t235 - (t188 * t235 * t247 + t182 * t243) * t354 + t373) * m(4) + (-m(4) * t319 - m(5) * t276 - m(6) * t256 + t370) * g(2) + (t305 + t306) * qJ(3) + t389 * t140 + t250 + qJD(1) * t282 + mrSges(3,2) * t278 + t69 * t345 + (t377 * t115 + t384 * t18 + t183 * t66 + t383 * t19 + t4 * t49 + t48 * t5) * m(6) + t384 * t85 - t168 * t304 + (-t278 + t317) * t394; -t267 * t37 + t185 * t38 + t242 * t89 + t246 * t88 + t379 * t85 + t378 * t84 + t268 * qJD(4) - t369 * t394 + (-t268 * t240 + (-t140 - t65) * t239) * t236 + t166 + (-g(1) * t228 + g(2) * t230) * (m(4) + m(5) + m(6)) + (-t115 * t331 + t18 * t379 + t185 * t4 + t19 * t378 - t267 * t5) * m(6) + (t242 * t26 + t246 * t27 - t161 + t196 * (-t242 * t67 + t246 * t68)) * m(5) + (-t235 * t341 + t150 - t161) * m(4); g(3) * t170 - t23 * t84 - t22 * t85 + t293 + (t241 * t4 + t245 * t5 + (-t18 * t241 + t19 * t245) * qJD(5)) * t368 + t19 * t351 + t252 - t65 * t281 - m(6) * (t115 * t281 + t18 * t22 + t19 * t23) + t95 * t296 / 0.2e1 - (-m(6) * t361 + t270) * t360 + (t139 + t279) * t68 + (-t138 - t280) * t67 - t387 * t369 + (-t156 * mrSges(5,2) + t155 * t382 - t322) * g(2) + (t158 * mrSges(5,2) - t157 * t382 - t321) * g(1) + (t94 * t287 - t388) * t236 + ((-t241 * t85 + t245 * t84) * qJD(5) + t241 * t38 + t245 * t37) * pkin(4) - t385; -t18 * t84 - g(1) * t321 - g(2) * t322 + t252 - t270 * t360 + (t85 + t351) * t19;];
tau = t1;
