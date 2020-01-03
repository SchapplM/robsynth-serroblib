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
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:05:30
% EndTime: 2020-01-03 12:05:40
% DurationCPUTime: 4.72s
% Computational Cost: add. (5392->436), mult. (8267->598), div. (0->0), fcn. (5206->14), ass. (0->238)
t244 = sin(pkin(9));
t245 = cos(pkin(9));
t193 = -pkin(3) * t245 - pkin(7) * t244 - pkin(2);
t247 = sin(qJ(4));
t251 = cos(qJ(4));
t319 = t245 * t251;
t142 = qJ(3) * t319 + t247 * t193;
t248 = sin(qJ(2));
t309 = qJD(3) * t245;
t252 = cos(qJ(2));
t318 = t245 * t252;
t342 = qJD(1) * pkin(1);
t362 = -qJD(4) * t142 - t247 * t309 - (-t247 * t318 + t248 * t251) * t342;
t306 = qJD(4) * t251;
t380 = (t247 * t248 + t251 * t318) * t342 - t193 * t306 - t251 * t309;
t307 = qJD(4) * t247;
t291 = t244 * t307;
t208 = pkin(8) * t291;
t379 = t208 + t362;
t320 = t245 * t247;
t297 = qJ(3) * t320;
t322 = t244 * t251;
t302 = pkin(8) * t322;
t378 = -(-t297 - t302) * qJD(4) + t380;
t246 = sin(qJ(5));
t250 = cos(qJ(5));
t186 = t246 * t251 + t247 * t250;
t162 = t186 * t244;
t241 = qJD(1) + qJD(2);
t116 = t241 * t162;
t238 = qJDD(1) + qJDD(2);
t143 = (t238 * t251 - t241 * t307) * t244;
t144 = (-t238 * t247 - t241 * t306) * t244;
t45 = -qJD(5) * t116 + t143 * t250 + t144 * t246;
t377 = Ifges(6,5) * t45;
t271 = t246 * t247 - t250 * t251;
t266 = t271 * qJD(5);
t325 = t241 * t244;
t46 = -t143 * t246 + t144 * t250 + t266 * t325;
t376 = Ifges(6,6) * t46;
t375 = Ifges(5,5) * t143;
t374 = Ifges(5,6) * t144;
t326 = t238 * t245;
t195 = qJDD(4) - t326;
t373 = Ifges(5,3) * t195;
t191 = qJDD(5) + t195;
t372 = Ifges(6,3) * t191;
t324 = t241 * t245;
t196 = qJD(4) - t324;
t189 = qJ(3) * t241 + t248 * t342;
t239 = t244 ^ 2;
t334 = t189 * t239;
t371 = t196 * t244 * (-Ifges(5,5) * t247 - Ifges(5,6) * t251) / 0.2e1 + (mrSges(5,1) * t251 - mrSges(5,2) * t247) * t334;
t370 = mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t346 = Ifges(5,4) * t251;
t347 = Ifges(5,4) * t247;
t369 = (-t247 * (-Ifges(5,2) * t251 - t347) / 0.2e1 + t251 * (-Ifges(5,1) * t247 - t346) / 0.2e1) * t239;
t298 = qJD(2) * t342;
t337 = qJDD(1) * pkin(1);
t181 = t248 * t337 + t252 * t298;
t128 = qJ(3) * t238 + qJD(3) * t241 + t181;
t300 = t252 * t342;
t277 = qJD(3) - t300;
t113 = t193 * t241 + t277;
t68 = t113 * t247 + t189 * t319;
t180 = -t248 * t298 + t252 * t337;
t269 = qJDD(3) - t180;
t87 = t193 * t238 + t269;
t27 = -qJD(4) * t68 - t128 * t320 + t251 * t87;
t16 = pkin(4) * t195 - pkin(8) * t143 + t27;
t290 = t245 * t307;
t26 = t113 * t306 + t128 * t319 - t189 * t290 + t247 * t87;
t17 = pkin(8) * t144 + t26;
t323 = t244 * t247;
t295 = t241 * t323;
t59 = -pkin(8) * t295 + t68;
t340 = t246 * t59;
t294 = t241 * t322;
t67 = t251 * t113 - t189 * t320;
t58 = -pkin(8) * t294 + t67;
t47 = pkin(4) * t196 + t58;
t18 = t250 * t47 - t340;
t4 = t18 * qJD(5) + t16 * t246 + t17 * t250;
t339 = t250 * t59;
t19 = t246 * t47 + t339;
t5 = -t19 * qJD(5) + t16 * t250 - t17 * t246;
t368 = t5 * mrSges(6,1) - t4 * mrSges(6,2);
t367 = -t27 * mrSges(5,1) + t26 * mrSges(5,2);
t303 = pkin(8) * t323;
t110 = -t303 + t142;
t179 = t251 * t193;
t91 = -t302 + t179 + (-qJ(3) * t247 - pkin(4)) * t245;
t48 = -t110 * t246 + t250 * t91;
t366 = qJD(5) * t48 + t246 * t379 - t378 * t250;
t49 = t110 * t250 + t246 * t91;
t365 = -qJD(5) * t49 + t378 * t246 + t250 * t379;
t355 = m(6) * pkin(4);
t364 = -mrSges(5,1) - t355;
t363 = -qJ(3) * t290 - t380;
t357 = -qJD(4) - qJD(5);
t105 = t357 * t186;
t361 = t186 * t324 + t105;
t360 = (t324 + t357) * t271;
t351 = pkin(1) * t252;
t173 = t193 - t351;
t224 = pkin(1) * t248 + qJ(3);
t101 = t247 * t173 + t224 * t319;
t243 = qJ(1) + qJ(2);
t232 = sin(t243);
t331 = t232 * t245;
t332 = t232 * t244;
t359 = -pkin(3) * t331 - pkin(7) * t332;
t310 = qJD(3) * t189;
t336 = t128 * t239;
t358 = qJ(3) * t336 + t239 * t310 - t300 * t334;
t356 = t241 ^ 2;
t118 = -t246 * t295 + t250 * t294;
t353 = t118 / 0.2e1;
t350 = pkin(4) * t247;
t349 = g(1) * t244;
t348 = mrSges(6,3) * t116;
t345 = Ifges(6,4) * t118;
t344 = pkin(1) * qJD(2);
t341 = t118 * mrSges(6,3);
t213 = t252 * t344 + qJD(3);
t338 = t213 * t334 + t224 * t336;
t335 = t128 * t244;
t333 = t189 * t241;
t330 = t232 * t247;
t234 = cos(t243);
t329 = t234 * t244;
t328 = t234 * t245;
t327 = t238 * t244;
t321 = t244 * (-pkin(8) - pkin(7));
t242 = qJ(4) + qJ(5);
t231 = sin(t242);
t233 = cos(t242);
t131 = -t231 * t331 - t233 * t234;
t132 = -t231 * t234 + t233 * t331;
t316 = t131 * mrSges(6,1) - t132 * mrSges(6,2);
t133 = t231 * t328 - t232 * t233;
t134 = t231 * t232 + t233 * t328;
t315 = t133 * mrSges(6,1) + t134 * mrSges(6,2);
t225 = t232 * pkin(2);
t249 = sin(qJ(1));
t235 = t249 * pkin(1);
t313 = t225 + t235;
t312 = t234 * pkin(2) + t232 * qJ(3);
t240 = t245 ^ 2;
t311 = t239 + t240;
t305 = m(4) + m(5) + m(6);
t304 = m(6) * t350;
t301 = t372 + t376 + t377;
t299 = t248 * t344;
t296 = t224 * t320;
t293 = t373 + t374 + t375;
t292 = t173 * t306 + t213 * t319 + t247 * t299;
t289 = t322 / 0.2e1;
t275 = mrSges(5,1) * t247 + mrSges(5,2) * t251;
t140 = t275 * t325;
t65 = mrSges(6,1) * t116 + mrSges(6,2) * t118;
t287 = (-t140 - t65) * t244;
t286 = t311 * mrSges(4,3);
t285 = t311 * t238;
t284 = t311 * t241;
t283 = -t234 * qJ(3) + t225;
t282 = pkin(4) * t294;
t281 = mrSges(5,3) * t295;
t280 = mrSges(5,3) * t294;
t279 = t68 * mrSges(5,3) * t306;
t278 = pkin(3) * t328 + pkin(7) * t329 + t312;
t166 = -mrSges(4,1) * t326 + mrSges(4,2) * t327;
t276 = -mrSges(4,1) * t245 + t244 * mrSges(4,2);
t274 = -mrSges(6,1) * t231 - mrSges(6,2) * t233;
t165 = t251 * t173;
t75 = -t302 + t165 + (-t224 * t247 - pkin(4)) * t245;
t86 = -t303 + t101;
t39 = -t246 * t86 + t250 * t75;
t40 = t246 * t75 + t250 * t86;
t138 = -mrSges(5,2) * t196 - t281;
t139 = mrSges(5,1) * t196 - t280;
t272 = t138 * t251 - t139 * t247;
t227 = pkin(4) * t251 + pkin(3);
t270 = t227 * t331 - t232 * t321 + t225;
t157 = -t232 * t251 + t234 * t320;
t155 = -t232 * t320 - t234 * t251;
t268 = (t251 * Ifges(5,1) - t347) * t244;
t267 = (-t247 * Ifges(5,2) + t346) * t244;
t261 = pkin(4) * t330 + t227 * t328 - t234 * t321 + t312;
t156 = t232 * t319 - t234 * t247;
t258 = -t232 * mrSges(3,1) - mrSges(4,1) * t331 - t156 * mrSges(5,1) - t132 * mrSges(6,1) - t234 * mrSges(3,2) - t155 * mrSges(5,2) - t131 * mrSges(6,2) + t370 * t332;
t63 = -qJD(4) * t101 - t213 * t320 + t251 * t299;
t114 = Ifges(6,4) * t116;
t115 = (t241 * t350 + t189) * t244;
t192 = qJD(5) + t196;
t51 = -Ifges(6,2) * t116 + Ifges(6,6) * t192 + t345;
t52 = Ifges(6,1) * t118 + Ifges(6,5) * t192 - t114;
t257 = -t115 * (mrSges(6,1) * t118 - mrSges(6,2) * t116) - t18 * t348 + t51 * t353 + t301 - t118 * (-Ifges(6,1) * t116 - t345) / 0.2e1 - t192 * (-Ifges(6,5) * t116 - Ifges(6,6) * t118) / 0.2e1 + (-Ifges(6,2) * t118 - t114 + t52) * t116 / 0.2e1 + t368;
t158 = t234 * t319 + t330;
t256 = -t234 * mrSges(3,1) - mrSges(4,1) * t328 - t158 * mrSges(5,1) - t134 * mrSges(6,1) + t157 * mrSges(5,2) + t133 * mrSges(6,2) + (mrSges(3,2) - mrSges(4,3)) * t232 + t370 * t329;
t150 = -pkin(2) * t238 + t269;
t170 = t275 * t244;
t66 = -pkin(4) * t144 + t335;
t78 = t105 * t244;
t79 = (qJD(4) * t271 + t266) * t244;
t94 = Ifges(5,6) * t196 + t241 * t267;
t95 = Ifges(5,5) * t196 + t241 * t268;
t255 = (Ifges(6,1) * t78 + Ifges(6,4) * t79) * t353 + t170 * t335 + t144 * t267 / 0.2e1 + (-t18 * t78 + t19 * t79) * mrSges(6,3) - (Ifges(5,4) * t143 + Ifges(5,2) * t144 + Ifges(5,6) * t195) * t323 / 0.2e1 + t150 * t276 + (mrSges(6,1) * t66 - mrSges(6,3) * t4 - Ifges(6,4) * t45 - Ifges(6,2) * t46 - Ifges(6,6) * t191) * t162 - t116 * (Ifges(6,4) * t78 + Ifges(6,2) * t79) / 0.2e1 + t115 * (-mrSges(6,1) * t79 + mrSges(6,2) * t78) + t78 * t52 / 0.2e1 + t79 * t51 / 0.2e1 + t180 * mrSges(3,1) - t181 * mrSges(3,2) + t192 * (Ifges(6,5) * t78 + Ifges(6,6) * t79) / 0.2e1 + Ifges(3,3) * t238 + (Ifges(5,1) * t143 + Ifges(5,4) * t144 + Ifges(5,5) * t195) * t289 + (-t373 / 0.2e1 - t372 / 0.2e1 - t376 / 0.2e1 - t377 / 0.2e1 + Ifges(4,2) * t326 + Ifges(4,4) * t327 - t375 / 0.2e1 - t374 / 0.2e1 + t367 - t368) * t245 + (t128 * t240 + t336) * mrSges(4,3) + (-t26 * t323 - t27 * t322 + t291 * t67) * mrSges(5,3) - (t301 + t293) * t245 / 0.2e1 + t143 * t268 / 0.2e1 + (t195 * (Ifges(5,5) * t251 - Ifges(5,6) * t247) / 0.2e1 + Ifges(4,4) * t326 + Ifges(4,1) * t327 + (-mrSges(6,2) * t66 + mrSges(6,3) * t5 - Ifges(6,1) * t45 - Ifges(6,4) * t46 - Ifges(6,5) * t191) * t271) * t244 + (t369 * t241 + t371 - (t247 * t95 + t251 * t94) * t244 / 0.2e1) * qJD(4);
t253 = cos(qJ(1));
t236 = t253 * pkin(1);
t228 = -pkin(2) - t351;
t216 = pkin(4) * t323;
t209 = t244 * pkin(4) * t306;
t184 = qJ(3) * t244 + t216;
t183 = -pkin(2) * t241 + t277;
t174 = qJD(3) * t244 + t209;
t168 = t276 * t241;
t167 = t224 * t244 + t216;
t161 = t239 * t333;
t149 = t213 * t244 + t209;
t141 = t179 - t297;
t100 = t165 - t296;
t89 = -mrSges(5,2) * t195 + mrSges(5,3) * t144;
t88 = mrSges(5,1) * t195 - mrSges(5,3) * t143;
t85 = mrSges(6,1) * t192 - t341;
t84 = -mrSges(6,2) * t192 - t348;
t69 = -mrSges(5,1) * t144 + mrSges(5,2) * t143;
t62 = -t224 * t290 + t292;
t56 = t208 + t63;
t55 = (-t296 - t302) * qJD(4) + t292;
t38 = -mrSges(6,2) * t191 + mrSges(6,3) * t46;
t37 = mrSges(6,1) * t191 - mrSges(6,3) * t45;
t23 = t250 * t58 - t340;
t22 = -t246 * t58 - t339;
t15 = -mrSges(6,1) * t46 + mrSges(6,2) * t45;
t7 = -qJD(5) * t40 - t246 * t55 + t250 * t56;
t6 = qJD(5) * t39 + t246 * t56 + t250 * t55;
t1 = [(-m(5) * (t313 - t359) + t258 + (qJ(3) * t305 + mrSges(4,3) + t304) * t234 - t253 * mrSges(2,2) - t249 * mrSges(2,1) - m(4) * t313 - m(6) * (t235 + t270)) * g(3) + (t213 * t284 + t224 * t285) * mrSges(4,3) + (t140 * t213 + t224 * t69 - t279) * t244 + t40 * t38 + t39 * t37 + m(4) * (t183 * t299 + t150 * t228 + (t128 * t224 + t189 * t213) * t240 + t338) + (t256 - t253 * mrSges(2,1) - m(6) * (t236 + t261) - m(5) * (t236 + t278) - m(4) * (t236 + t312) + t249 * mrSges(2,2)) * g(2) + t255 + t62 * t138 + t63 * t139 + t101 * t89 + t6 * t84 + t7 * t85 + t100 * t88 + t149 * t65 + t167 * t15 + t228 * t166 + m(6) * (t115 * t149 + t167 * t66 + t18 * t7 + t19 * t6 + t39 * t5 + t4 * t40) + m(5) * (t100 * t27 + t101 * t26 + t62 * t68 + t63 * t67 + t338) + Ifges(2,3) * qJDD(1) + ((mrSges(3,1) * t252 - mrSges(3,2) * t248) * t238 + (t168 * t248 + (-mrSges(3,1) * t248 - mrSges(3,2) * t252) * t241) * qJD(2) + (-g(2) * t253 - g(3) * t249 + t180 * t252 + t181 * t248) * m(3)) * pkin(1); t365 * t85 + t366 * t84 + (g(3) * t234 + qJ(3) * t285 + qJD(3) * t284) * mrSges(4,3) + t48 * t37 + t49 * t38 + t256 * g(2) + t258 * g(3) + (qJ(3) * t69 + qJD(3) * t140 - t279) * t244 + t362 * t139 + t363 * t138 + t255 + t141 * t88 + t142 * t89 - pkin(2) * t166 + t174 * t65 + t184 * t15 + ((mrSges(3,1) * t241 - t168) * t248 + (t287 + (mrSges(3,2) - t286) * t241) * t252) * t342 + (t184 * t66 + t4 * t49 + t48 * t5 - t261 * g(2) + (-(-qJ(3) - t350) * t234 - t270) * g(3) + t366 * t19 + t365 * t18 + (-t244 * t300 + t174) * t115) * m(6) + (t141 * t27 + t142 * t26 - t278 * g(2) + (-t283 + t359) * g(3) + t363 * t68 + t362 * t67 + t358) * m(5) + (-t312 * g(2) - t283 * g(3) - pkin(2) * t150 + (qJ(3) * t128 + t310) * t240 - (t189 * t240 * t252 + t183 * t248) * t342 + t358) * m(4); -t271 * t37 + t186 * t38 + t247 * t89 + t251 * t88 + t361 * t85 + t360 * t84 + t272 * qJD(4) - t356 * t286 + (-t245 * t272 + t287) * t241 + t166 + (g(2) * t234 + g(3) * t232) * t305 + (-t115 * t325 + t18 * t361 + t4 * t186 + t19 * t360 - t271 * t5) * m(6) + (t247 * t26 + t251 * t27 - t161 + t196 * (-t247 * t67 + t251 * t68)) * m(5) + (-t240 * t333 + t150 - t161) * m(4); (t246 * t4 + t250 * t5 + (-t18 * t246 + t19 * t250) * qJD(5)) * t355 + t19 * t341 + t257 - (t274 - t304) * t349 + t95 * t295 / 0.2e1 - t65 * t282 - m(6) * (t115 * t282 + t18 * t22 + t19 * t23) + t293 - t23 * t84 - t22 * t85 + g(1) * t170 + (t280 + t139) * t68 + (-t281 - t138) * t67 - t369 * t356 + (-t158 * mrSges(5,2) + t157 * t364 - t315) * g(3) + (t156 * mrSges(5,2) + t155 * t364 - t316) * g(2) + (t94 * t289 - t371) * t241 + ((-t246 * t85 + t250 * t84) * qJD(5) + t246 * t38 + t250 * t37) * pkin(4) - t367; t257 - g(2) * t316 - g(3) * t315 - t18 * t84 - t274 * t349 + (t85 + t341) * t19;];
tau = t1;
