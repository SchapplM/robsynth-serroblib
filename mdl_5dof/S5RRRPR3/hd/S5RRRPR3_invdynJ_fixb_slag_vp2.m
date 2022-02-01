% Calculate vector of inverse dynamics joint torques for
% S5RRRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2022-01-20 11:44
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR3_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR3_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:42:28
% EndTime: 2022-01-20 11:42:40
% DurationCPUTime: 5.05s
% Computational Cost: add. (6393->444), mult. (9682->590), div. (0->0), fcn. (6446->16), ass. (0->226)
t262 = sin(qJ(3));
t342 = t262 / 0.2e1;
t266 = cos(qJ(3));
t244 = t266 * qJD(4);
t260 = -qJ(4) - pkin(7);
t288 = qJD(3) * t260;
t174 = t262 * t288 + t244;
t175 = -qJD(4) * t262 + t266 * t288;
t258 = sin(pkin(9));
t259 = cos(pkin(9));
t189 = t258 * t266 + t259 * t262;
t267 = cos(qJ(2));
t323 = pkin(1) * qJD(1);
t294 = t267 * t323;
t355 = -t174 * t258 + t259 * t175 + t189 * t294;
t188 = -t258 * t262 + t259 * t266;
t354 = t259 * t174 + t258 * t175 - t188 * t294;
t254 = qJ(3) + pkin(9);
t243 = qJ(5) + t254;
t227 = sin(t243);
t228 = cos(t243);
t378 = t228 * mrSges(6,1) - t227 * mrSges(6,2);
t210 = -t266 * mrSges(4,1) + t262 * mrSges(4,2);
t365 = -mrSges(3,1) + t210;
t178 = t188 * qJD(3);
t333 = pkin(8) * t178;
t377 = -t333 + t355;
t177 = t189 * qJD(3);
t170 = t177 * pkin(8);
t376 = t170 - t354;
t241 = sin(t254);
t242 = cos(t254);
t375 = -t242 * mrSges(5,1) + t241 * mrSges(5,2) - t378;
t374 = t365 + t375;
t261 = sin(qJ(5));
t265 = cos(qJ(5));
t253 = qJD(1) + qJD(2);
t160 = t189 * t253;
t334 = pkin(8) * t160;
t263 = sin(qJ(2));
t295 = t263 * t323;
t200 = pkin(7) * t253 + t295;
t287 = qJ(4) * t253 + t200;
t145 = t287 * t266;
t129 = t258 * t145;
t144 = t287 * t262;
t321 = qJD(3) * pkin(3);
t135 = -t144 + t321;
t75 = t259 * t135 - t129;
t50 = qJD(3) * pkin(4) - t334 + t75;
t159 = t188 * t253;
t335 = pkin(8) * t159;
t303 = t259 * t145;
t76 = t258 * t135 + t303;
t52 = t76 + t335;
t22 = -t261 * t52 + t265 * t50;
t23 = t261 * t50 + t265 * t52;
t249 = qJDD(3) + qJDD(5);
t252 = qJD(3) + qJD(5);
t286 = t265 * t159 - t160 * t261;
t250 = qJDD(1) + qJDD(2);
t297 = qJD(3) * t262;
t179 = t250 * t266 - t253 * t297;
t296 = qJD(3) * t266;
t180 = t250 * t262 + t253 * t296;
t106 = t179 * t258 + t180 * t259;
t322 = pkin(1) * qJD(2);
t292 = qJD(1) * t322;
t313 = pkin(1) * qJDD(1);
t187 = t263 * t313 + t267 * t292;
t166 = pkin(7) * t250 + t187;
t291 = t200 * t296;
t65 = -t291 + qJDD(3) * pkin(3) - qJ(4) * t180 + (-qJD(4) * t253 - t166) * t262;
t108 = t266 * t166 - t200 * t297;
t70 = qJ(4) * t179 + t244 * t253 + t108;
t31 = -t258 * t70 + t259 * t65;
t18 = qJDD(3) * pkin(4) - pkin(8) * t106 + t31;
t105 = t179 * t259 - t180 * t258;
t32 = t258 * t65 + t259 * t70;
t19 = pkin(8) * t105 + t32;
t3 = qJD(5) * t22 + t18 * t261 + t19 * t265;
t92 = t159 * t261 + t160 * t265;
t340 = Ifges(6,4) * t92;
t36 = qJD(5) * t286 + t105 * t261 + t106 * t265;
t37 = -qJD(5) * t92 + t105 * t265 - t106 * t261;
t4 = -qJD(5) * t23 + t18 * t265 - t19 * t261;
t86 = Ifges(6,4) * t286;
t41 = Ifges(6,1) * t92 + Ifges(6,5) * t252 + t86;
t336 = pkin(3) * t266;
t234 = pkin(2) + t336;
t156 = -t234 * t253 + qJD(4) - t294;
t99 = -pkin(4) * t159 + t156;
t373 = t4 * mrSges(6,1) - t3 * mrSges(6,2) + Ifges(6,5) * t36 + Ifges(6,6) * t37 + Ifges(6,3) * t249 - (Ifges(6,5) * t286 - Ifges(6,6) * t92) * t252 / 0.2e1 + (t22 * t286 + t23 * t92) * mrSges(6,3) - (-Ifges(6,2) * t92 + t41 + t86) * t286 / 0.2e1 - t99 * (mrSges(6,1) * t92 + mrSges(6,2) * t286) - (Ifges(6,1) * t286 - t340) * t92 / 0.2e1;
t328 = qJD(3) / 0.2e1;
t229 = pkin(3) * t259 + pkin(4);
t338 = pkin(3) * t258;
t173 = t229 * t261 + t265 * t338;
t80 = t144 * t258 - t303;
t54 = t80 - t335;
t81 = -t259 * t144 - t129;
t55 = t81 - t334;
t369 = -t173 * qJD(5) + t261 * t55 - t265 * t54;
t172 = t229 * t265 - t261 * t338;
t368 = t172 * qJD(5) - t261 * t54 - t265 * t55;
t306 = t253 * t262;
t198 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t306;
t305 = t253 * t266;
t199 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t305;
t367 = (t253 * mrSges(3,2) + t262 * t198 - t266 * t199) * t267;
t257 = qJ(1) + qJ(2);
t245 = sin(t257);
t246 = cos(t257);
t366 = g(1) * t246 + g(2) * t245;
t364 = mrSges(3,2) - mrSges(6,3) - mrSges(5,3) - mrSges(4,3);
t40 = Ifges(6,2) * t286 + Ifges(6,6) * t252 + t340;
t362 = t40 / 0.2e1;
t209 = t260 * t262;
t247 = t266 * qJ(4);
t211 = pkin(7) * t266 + t247;
t133 = t259 * t209 - t211 * t258;
t332 = pkin(8) * t189;
t103 = t133 - t332;
t134 = t258 * t209 + t259 * t211;
t183 = t188 * pkin(8);
t104 = t183 + t134;
t45 = t103 * t265 - t104 * t261;
t361 = qJD(5) * t45 + t377 * t261 - t376 * t265;
t46 = t103 * t261 + t104 * t265;
t360 = -qJD(5) * t46 + t376 * t261 + t377 * t265;
t352 = t365 * t253;
t201 = -pkin(2) * t253 - t294;
t350 = t201 * t263 + (t262 ^ 2 + t266 ^ 2) * t200 * t267;
t349 = m(5) * pkin(3);
t346 = t92 / 0.2e1;
t344 = t160 / 0.2e1;
t339 = pkin(1) * t267;
t337 = pkin(3) * t262;
t329 = -qJD(3) / 0.2e1;
t327 = mrSges(4,2) * t266;
t326 = Ifges(4,4) * t262;
t325 = Ifges(4,4) * t266;
t324 = Ifges(5,4) * t160;
t315 = t266 * Ifges(4,2);
t312 = t108 * t266;
t109 = -t166 * t262 - t291;
t311 = t109 * t262;
t233 = pkin(1) * t263 + pkin(7);
t307 = t233 * t266;
t299 = -qJ(4) - t233;
t283 = qJD(3) * t299;
t293 = t267 * t322;
t121 = t262 * t283 + t266 * t293 + t244;
t122 = (-qJD(4) - t293) * t262 + t266 * t283;
t67 = t259 * t121 + t258 * t122;
t184 = t299 * t262;
t185 = t247 + t307;
t112 = t258 * t184 + t259 * t185;
t298 = t246 * pkin(2) + t245 * pkin(7);
t237 = pkin(3) * t297;
t9 = -t37 * mrSges(6,1) + t36 * mrSges(6,2);
t290 = t296 / 0.2e1;
t289 = pkin(4) * t242 + t336;
t142 = pkin(4) * t177 + t237;
t49 = -t105 * mrSges(5,1) + t106 * mrSges(5,2);
t66 = -t121 * t258 + t259 * t122;
t111 = t259 * t184 - t185 * t258;
t196 = pkin(2) + t289;
t251 = pkin(8) - t260;
t285 = t246 * t196 + t245 * t251;
t284 = t246 * t234 - t245 * t260;
t186 = -t263 * t292 + t267 * t313;
t282 = -g(1) * t245 + g(2) * t246;
t281 = mrSges(6,1) * t227 + mrSges(6,2) * t228;
t279 = t315 + t326;
t278 = Ifges(4,5) * t266 - Ifges(4,6) * t262;
t84 = t111 - t332;
t85 = t183 + t112;
t38 = -t261 * t85 + t265 * t84;
t39 = t261 * t84 + t265 * t85;
t115 = t188 * t265 - t189 * t261;
t116 = t188 * t261 + t189 * t265;
t277 = t198 * t266 + t199 * t262;
t152 = -pkin(4) * t188 - t234;
t275 = t201 * (mrSges(4,1) * t262 + t327);
t274 = t262 * (Ifges(4,1) * t266 - t326);
t165 = -pkin(2) * t250 - t186;
t110 = -pkin(3) * t179 + qJDD(4) + t165;
t273 = -qJD(3) * t277 + t266 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t179) - t262 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t180);
t272 = t364 * t245 + t374 * t246;
t270 = (m(4) * pkin(2) + m(5) * t234 + m(6) * t196 - t374) * t245 + (-m(4) * pkin(7) + m(5) * t260 - m(6) * t251 + t364) * t246;
t163 = Ifges(4,6) * qJD(3) + t253 * t279;
t215 = Ifges(4,4) * t305;
t164 = Ifges(4,1) * t306 + Ifges(4,5) * qJD(3) + t215;
t51 = -pkin(4) * t105 + t110;
t57 = qJD(5) * t115 - t177 * t261 + t178 * t265;
t58 = -qJD(5) * t116 - t177 * t265 - t178 * t261;
t87 = Ifges(5,2) * t159 + Ifges(5,6) * qJD(3) + t324;
t149 = Ifges(5,4) * t159;
t88 = Ifges(5,1) * t160 + Ifges(5,5) * qJD(3) + t149;
t269 = (0.2e1 * Ifges(4,5) * t342 + Ifges(5,5) * t189 + Ifges(4,6) * t266 + Ifges(5,6) * t188) * qJDD(3) + (t274 * t328 + (-Ifges(4,2) * t262 + t325) * t290) * t253 + t179 * t279 / 0.2e1 + t164 * t290 + t266 * (Ifges(4,4) * t180 + Ifges(4,2) * t179) / 0.2e1 + (t278 * t328 + t275) * qJD(3) + (Ifges(4,1) * t180 + Ifges(4,4) * t179) * t342 + (-t22 * t57 + t23 * t58) * mrSges(6,3) + t58 * t362 + t286 * (Ifges(6,4) * t57 + Ifges(6,2) * t58) / 0.2e1 + (Ifges(5,1) * t178 - Ifges(5,4) * t177) * t344 + (Ifges(5,5) * t178 - Ifges(5,6) * t177) * t328 + t156 * (mrSges(5,1) * t177 + mrSges(5,2) * t178) + t159 * (Ifges(5,4) * t178 - Ifges(5,2) * t177) / 0.2e1 + (-mrSges(6,1) * t51 + mrSges(6,3) * t3 + Ifges(6,4) * t36 + Ifges(6,2) * t37 + Ifges(6,6) * t249) * t115 - t163 * t297 / 0.2e1 + (Ifges(6,1) * t57 + Ifges(6,4) * t58) * t346 + mrSges(4,3) * t312 + (mrSges(6,2) * t51 - mrSges(6,3) * t4 + Ifges(6,1) * t36 + Ifges(6,4) * t37 + Ifges(6,5) * t249) * t116 + t180 * (t262 * Ifges(4,1) + t325) / 0.2e1 + (-t177 * t76 - t75 * t178 + t188 * t32 - t31 * t189) * mrSges(5,3) + t57 * t41 / 0.2e1 + t99 * (-mrSges(6,1) * t58 + mrSges(6,2) * t57) - t177 * t87 / 0.2e1 + t178 * t88 / 0.2e1 + t186 * mrSges(3,1) - t187 * mrSges(3,2) + t110 * (-mrSges(5,1) * t188 + mrSges(5,2) * t189) + t106 * (Ifges(5,1) * t189 + Ifges(5,4) * t188) + t105 * (Ifges(5,4) * t189 + Ifges(5,2) * t188) + t165 * t210 + Ifges(3,3) * t250 + t252 * (Ifges(6,5) * t57 + Ifges(6,6) * t58) / 0.2e1;
t268 = cos(qJ(1));
t264 = sin(qJ(1));
t248 = t268 * pkin(1);
t238 = t263 * t322;
t235 = -pkin(2) - t339;
t206 = -t234 - t339;
t197 = t238 + t237;
t140 = t152 - t339;
t137 = qJD(3) * mrSges(5,1) - t160 * mrSges(5,3);
t136 = -qJD(3) * mrSges(5,2) + t159 * mrSges(5,3);
t125 = t142 + t238;
t120 = pkin(3) * t306 + pkin(4) * t160;
t113 = -mrSges(4,1) * t179 + mrSges(4,2) * t180;
t95 = -mrSges(5,1) * t159 + mrSges(5,2) * t160;
t94 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t106;
t93 = -qJDD(3) * mrSges(5,2) + mrSges(5,3) * t105;
t79 = mrSges(6,1) * t252 - mrSges(6,3) * t92;
t78 = -mrSges(6,2) * t252 + mrSges(6,3) * t286;
t48 = -t170 + t67;
t47 = t66 - t333;
t44 = -mrSges(6,1) * t286 + mrSges(6,2) * t92;
t30 = -mrSges(6,2) * t249 + mrSges(6,3) * t37;
t29 = mrSges(6,1) * t249 - mrSges(6,3) * t36;
t8 = -qJD(5) * t39 - t261 * t48 + t265 * t47;
t7 = qJD(5) * t38 + t261 * t47 + t265 * t48;
t1 = [t273 * t233 + (t264 * mrSges(2,1) + mrSges(2,2) * t268 + t270) * g(1) + m(5) * (t110 * t206 + t111 * t31 + t112 * t32 + t156 * t197 + t66 * t75 + t67 * t76) + m(6) * (t125 * t99 + t140 * t51 + t22 * t8 + t23 * t7 + t3 * t39 + t38 * t4) + (-m(6) * (t248 + t285) - m(5) * (t248 + t284) - m(4) * (t248 + t298) - mrSges(2,1) * t268 + t264 * mrSges(2,2) + t272) * g(2) + t269 + m(4) * (t108 * t307 + t165 * t235 - t233 * t311) + ((mrSges(3,1) * t267 - mrSges(3,2) * t263) * t250 + (-g(2) * t268 + t186 * t267 + t187 * t263) * m(3) + (m(3) + m(4) + m(5) + m(6)) * t264 * g(1) + (m(4) * t350 + t263 * t352 - t367) * qJD(2)) * pkin(1) - mrSges(4,3) * t311 + t38 * t29 + t39 * t30 + t7 * t78 + t8 * t79 + t111 * t94 + t112 * t93 + t125 * t44 + t67 * t136 + t66 * t137 + t140 * t9 + t197 * t95 + t206 * t49 + t235 * t113 + Ifges(2,3) * qJDD(1); t360 * t79 + t361 * t78 + (t367 + (-t352 - t44 - t95) * t263) * t323 + t273 * pkin(7) + t270 * g(1) + t272 * g(2) + t269 + t354 * t136 + t355 * t137 + t45 * t29 + t46 * t30 - pkin(2) * t113 + t133 * t94 + t134 * t93 + t142 * t44 + t152 * t9 - t234 * t49 + (-t109 * mrSges(4,3) + t321 * t95) * t262 + (-t285 * g(2) + t152 * t51 + t3 * t46 + t4 * t45 + (t142 - t295) * t99 + t361 * t23 + t360 * t22) * m(6) + (-g(2) * t284 - t110 * t234 + t133 * t31 + t134 * t32 + t354 * t76 + t355 * t75 + (t237 - t295) * t156) * m(5) + (-t298 * g(2) - pkin(2) * t165 + (-t311 + t312) * pkin(7) - t350 * t323) * m(4); t373 + (-m(5) * t336 + t210 + t375) * g(3) + (Ifges(5,3) + Ifges(4,3)) * qJDD(3) + t368 * t78 + (-t289 * g(3) - t120 * t99 + t172 * t4 + t173 * t3 + t369 * t22 + t368 * t23) * m(6) + t369 * t79 + t366 * (-m(6) * (-pkin(4) * t241 - t337) + mrSges(5,1) * t241 + t327 + mrSges(5,2) * t242 + (mrSges(4,1) + t349) * t262 + t281) + (t163 * t342 - t275 + t278 * t329 + (t315 * t342 - t274 / 0.2e1) * t253 + (-m(5) * t156 - t95) * t337 - (t215 + t164) * t266 / 0.2e1) * t253 - m(5) * (t75 * t80 + t76 * t81) + t87 * t344 + (t258 * t32 + t259 * t31) * t349 + t277 * t200 + (t258 * t93 + t259 * t94) * pkin(3) + t31 * mrSges(5,1) - t32 * mrSges(5,2) + Ifges(5,6) * t105 + Ifges(5,5) * t106 - t108 * mrSges(4,2) + t109 * mrSges(4,1) - t120 * t44 - t81 * t136 - t80 * t137 + t172 * t29 + t173 * t30 + Ifges(4,6) * t179 + Ifges(4,5) * t180 + t92 * t362 + (t159 * t75 + t160 * t76) * mrSges(5,3) + (Ifges(5,5) * t159 - Ifges(5,6) * t160) * t329 - t156 * (t160 * mrSges(5,1) + t159 * mrSges(5,2)) - (-Ifges(5,2) * t160 + t149 + t88) * t159 / 0.2e1 - t160 * (Ifges(5,1) * t159 - t324) / 0.2e1; -t159 * t136 + t160 * t137 - t286 * t78 + t92 * t79 + t49 + t9 + (t22 * t92 - t23 * t286 + t282 + t51) * m(6) + (-t76 * t159 + t75 * t160 + t110 + t282) * m(5); -g(3) * t378 - t22 * t78 + t23 * t79 + t366 * t281 + t40 * t346 + t373;];
tau = t1;
