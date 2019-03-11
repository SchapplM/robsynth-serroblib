% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRRPR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:00:19
% EndTime: 2019-03-08 23:00:39
% DurationCPUTime: 9.47s
% Computational Cost: add. (10524->551), mult. (27560->778), div. (0->0), fcn. (21021->12), ass. (0->273)
t245 = sin(qJ(4));
t246 = sin(qJ(3));
t249 = cos(qJ(4));
t250 = cos(qJ(3));
t215 = t245 * t250 + t246 * t249;
t251 = cos(qJ(2));
t242 = sin(pkin(6));
t314 = qJD(1) * t242;
t296 = t251 * t314;
t180 = t215 * t296;
t370 = -pkin(9) - pkin(8);
t227 = t370 * t246;
t228 = t370 * t250;
t183 = t245 * t227 - t249 * t228;
t299 = qJD(3) * t370;
t220 = t246 * t299;
t221 = t250 * t299;
t376 = -qJD(4) * t183 - t220 * t245 + t249 * t221 + t180;
t214 = -t245 * t246 + t249 * t250;
t305 = qJD(4) * t249;
t306 = qJD(4) * t245;
t375 = -t214 * t296 + t249 * t220 + t245 * t221 + t227 * t305 + t228 * t306;
t240 = qJD(3) + qJD(4);
t175 = t240 * t215;
t393 = qJ(5) * t175 - qJD(5) * t214 - t375;
t174 = t240 * t214;
t392 = -qJ(5) * t174 - qJD(5) * t215 + t376;
t241 = sin(pkin(12));
t326 = cos(pkin(12));
t391 = -t392 * t241 + t326 * t393;
t115 = t174 * t241 + t175 * t326;
t116 = t174 * t326 - t241 * t175;
t307 = qJD(3) * t246;
t239 = pkin(3) * t307;
t154 = pkin(4) * t175 + t239;
t247 = sin(qJ(2));
t297 = t247 * t314;
t390 = pkin(5) * t115 - pkin(10) * t116 + t154 - t297;
t244 = sin(qJ(6));
t248 = cos(qJ(6));
t207 = t214 * qJD(2);
t208 = t215 * qJD(2);
t260 = t241 * t207 + t208 * t326;
t134 = t240 * t248 - t244 * t260;
t287 = t326 * t207 - t208 * t241;
t143 = qJD(6) - t287;
t282 = mrSges(7,1) * t244 + mrSges(7,2) * t248;
t222 = qJD(2) * pkin(8) + t297;
t288 = pkin(9) * qJD(2) + t222;
t243 = cos(pkin(6));
t313 = qJD(1) * t246;
t295 = t243 * t313;
t179 = t250 * t288 + t295;
t171 = t249 * t179;
t312 = qJD(1) * t250;
t230 = t243 * t312;
t178 = -t288 * t246 + t230;
t172 = qJD(3) * pkin(3) + t178;
t118 = t172 * t245 + t171;
t323 = qJ(5) * t207;
t99 = t118 + t323;
t329 = t241 * t99;
t169 = t245 * t179;
t117 = t249 * t172 - t169;
t198 = t208 * qJ(5);
t98 = t117 - t198;
t85 = pkin(4) * t240 + t98;
t41 = t326 * t85 - t329;
t36 = -t240 * pkin(5) - t41;
t264 = t36 * t282;
t277 = Ifges(7,5) * t248 - Ifges(7,6) * t244;
t346 = Ifges(7,4) * t248;
t279 = -Ifges(7,2) * t244 + t346;
t347 = Ifges(7,4) * t244;
t281 = Ifges(7,1) * t248 - t347;
t357 = t248 / 0.2e1;
t358 = -t244 / 0.2e1;
t135 = t240 * t244 + t248 * t260;
t364 = t135 / 0.2e1;
t342 = t135 * Ifges(7,4);
t60 = t134 * Ifges(7,2) + t143 * Ifges(7,6) + t342;
t133 = Ifges(7,4) * t134;
t61 = t135 * Ifges(7,1) + t143 * Ifges(7,5) + t133;
t389 = t61 * t357 + t60 * t358 + t143 * t277 / 0.2e1 + t281 * t364 + t134 * t279 / 0.2e1 + t264;
t388 = -Ifges(4,1) / 0.2e1;
t308 = qJD(2) * t250;
t387 = -Ifges(4,4) * t308 / 0.2e1;
t167 = -t214 * t326 + t215 * t241;
t168 = t241 * t214 + t215 * t326;
t236 = -pkin(3) * t250 - pkin(2);
t189 = -pkin(4) * t214 + t236;
t86 = pkin(5) * t167 - pkin(10) * t168 + t189;
t146 = qJ(5) * t214 + t183;
t182 = t249 * t227 + t228 * t245;
t265 = -qJ(5) * t215 + t182;
t90 = t146 * t326 + t241 * t265;
t40 = t244 * t86 + t248 * t90;
t384 = -qJD(6) * t40 + t244 * t391 + t390 * t248;
t39 = -t244 * t90 + t248 * t86;
t383 = qJD(6) * t39 + t390 * t244 - t248 * t391;
t382 = mrSges(6,3) * t287;
t335 = t260 * Ifges(6,4);
t338 = t287 * Ifges(6,4);
t381 = t241 * t393 + t392 * t326;
t328 = mrSges(6,1) * t240 + mrSges(7,1) * t134 - mrSges(7,2) * t135 - mrSges(6,3) * t260;
t120 = t249 * t178 - t169;
t103 = -t198 + t120;
t119 = -t178 * t245 - t171;
t266 = t119 - t323;
t289 = t326 * t245;
t345 = pkin(3) * qJD(4);
t379 = -t103 * t241 + t266 * t326 + (t241 * t249 + t289) * t345;
t317 = t241 * t245;
t195 = (t249 * t326 - t317) * t345;
t49 = t103 * t326 + t241 * t266;
t378 = t195 - t49;
t164 = -mrSges(5,1) * t207 + mrSges(5,2) * t208;
t200 = qJD(2) * t236 - t296;
t377 = -m(5) * t200 - t164;
t185 = -t222 * t246 + t230;
t87 = t326 * t99;
t42 = t241 * t85 + t87;
t37 = pkin(10) * t240 + t42;
t162 = -pkin(4) * t207 + qJD(5) + t200;
t70 = -pkin(5) * t287 - pkin(10) * t260 + t162;
t16 = -t244 * t37 + t248 * t70;
t17 = t244 * t70 + t248 * t37;
t374 = -t16 * t244 + t17 * t248;
t160 = t174 * qJD(2);
t161 = t175 * qJD(2);
t110 = t160 * t241 + t161 * t326;
t111 = t160 * t326 - t241 * t161;
t310 = qJD(2) * t246;
t238 = pkin(3) * t310;
t311 = qJD(2) * t242;
t290 = qJD(1) * t311;
t201 = qJD(3) * t238 + t247 * t290;
t130 = pkin(4) * t161 + t201;
t32 = pkin(5) * t110 - pkin(10) * t111 + t130;
t286 = t251 * t290;
t269 = t250 * t286;
t270 = t246 * t286;
t52 = -t179 * t306 + t245 * (-qJD(3) * t179 - t270) + t249 * (qJD(3) * t178 + t269) + t172 * t305;
t25 = -qJ(5) * t161 + qJD(5) * t207 + t52;
t318 = t222 * t250;
t186 = t295 + t318;
t254 = (-t245 * (-pkin(9) * t310 + t185) + t249 * (-pkin(9) * t308 - t186)) * qJD(3);
t268 = -t160 * qJ(5) - t208 * qJD(5);
t9 = t326 * t25 + (-t172 * t306 - t179 * t305 - t245 * t269 - t249 * t270 + t254 + t268) * t241;
t2 = qJD(6) * t16 + t244 * t32 + t248 * t9;
t322 = qJD(6) * t17;
t3 = -t244 * t9 + t248 * t32 - t322;
t64 = qJD(6) * t134 + t111 * t248;
t65 = -qJD(6) * t135 - t111 * t244;
t373 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t64 + Ifges(7,6) * t65;
t372 = t64 / 0.2e1;
t371 = t65 / 0.2e1;
t53 = -qJD(2) * t180 - qJD(4) * t118 + t254;
t8 = t241 * t25 - t326 * (t268 + t53);
t89 = t146 * t241 - t265 * t326;
t369 = t8 * t89;
t316 = t242 * t247;
t203 = t243 * t250 - t246 * t316;
t204 = t243 * t246 + t250 * t316;
t148 = t203 * t249 - t204 * t245;
t149 = t203 * t245 + t204 * t249;
t95 = -t148 * t326 + t149 * t241;
t368 = t8 * t95;
t367 = t110 / 0.2e1;
t366 = -t134 / 0.2e1;
t365 = -t135 / 0.2e1;
t363 = -t143 / 0.2e1;
t361 = t207 / 0.2e1;
t360 = t208 / 0.2e1;
t355 = pkin(4) * t208;
t354 = pkin(4) * t241;
t353 = t16 * mrSges(7,3);
t352 = t2 * t248;
t351 = t3 * t244;
t350 = mrSges(5,3) * t207;
t349 = mrSges(5,3) * t208;
t348 = Ifges(4,4) * t246;
t344 = qJD(2) * pkin(2);
t343 = t134 * Ifges(7,6);
t341 = t135 * Ifges(7,5);
t340 = t143 * Ifges(7,3);
t339 = t260 * t42;
t337 = t287 * Ifges(6,2);
t336 = t260 * Ifges(6,1);
t332 = t208 * Ifges(5,4);
t331 = t240 * Ifges(6,5);
t330 = t240 * Ifges(6,6);
t325 = Ifges(4,5) * qJD(3);
t324 = Ifges(4,6) * qJD(3);
t321 = t287 * t244;
t320 = t287 * t248;
t315 = t242 * t251;
t235 = pkin(3) * t249 + pkin(4);
t197 = pkin(3) * t289 + t241 * t235;
t309 = qJD(2) * t247;
t304 = qJD(6) * t244;
t303 = qJD(6) * t248;
t302 = qJD(2) * qJD(3);
t298 = t326 * pkin(4);
t294 = t242 * t309;
t293 = t251 * t311;
t292 = t325 / 0.2e1;
t291 = -t324 / 0.2e1;
t54 = t110 * mrSges(6,1) + t111 * mrSges(6,2);
t284 = -t2 * t244 - t248 * t3;
t283 = mrSges(7,1) * t248 - mrSges(7,2) * t244;
t280 = Ifges(7,1) * t244 + t346;
t278 = Ifges(7,2) * t248 + t347;
t276 = Ifges(7,5) * t244 + Ifges(7,6) * t248;
t275 = -t16 * t248 - t17 * t244;
t81 = -mrSges(7,2) * t143 + mrSges(7,3) * t134;
t82 = mrSges(7,1) * t143 - mrSges(7,3) * t135;
t273 = -t244 * t82 + t248 * t81;
t261 = qJD(3) * t243 + t293;
t155 = -t222 * t307 + t261 * t312;
t156 = -qJD(3) * t318 - t261 * t313;
t271 = t155 * t250 - t156 * t246;
t96 = t241 * t148 + t149 * t326;
t78 = -t244 * t96 - t248 * t315;
t267 = t244 * t315 - t248 * t96;
t223 = -t296 - t344;
t263 = t185 * mrSges(4,3) + t310 * t388 + t387 - t325 / 0.2e1 - t223 * mrSges(4,2);
t262 = t186 * mrSges(4,3) + t324 / 0.2e1 + (t250 * Ifges(4,2) + t348) * qJD(2) / 0.2e1 - t223 * mrSges(4,1);
t80 = pkin(5) * t260 - pkin(10) * t287 + t355;
t196 = -pkin(3) * t317 + t235 * t326;
t256 = qJD(6) * t275 - t351 + t352;
t12 = t64 * Ifges(7,4) + t65 * Ifges(7,2) + t110 * Ifges(7,6);
t13 = t64 * Ifges(7,1) + t65 * Ifges(7,4) + t110 * Ifges(7,5);
t144 = t207 * Ifges(5,2) + t240 * Ifges(5,6) + t332;
t199 = Ifges(5,4) * t207;
t145 = t208 * Ifges(5,1) + t240 * Ifges(5,5) + t199;
t59 = t340 + t341 + t343;
t92 = t330 + t335 + t337;
t93 = t331 + t336 + t338;
t253 = t41 * t382 - t200 * (mrSges(5,1) * t208 + mrSges(5,2) * t207) - t17 * (-mrSges(7,2) * t260 - mrSges(7,3) * t321) - t16 * (mrSges(7,1) * t260 - mrSges(7,3) * t320) + t260 * t92 / 0.2e1 - t162 * (mrSges(6,1) * t260 + mrSges(6,2) * t287) + (Ifges(7,3) * t260 + t277 * t287) * t363 + (Ifges(7,5) * t260 + t281 * t287) * t365 + (Ifges(7,6) * t260 + t279 * t287) * t366 - (Ifges(5,5) * t207 + Ifges(6,5) * t287 - Ifges(5,6) * t208 - Ifges(6,6) * t260) * t240 / 0.2e1 - (-Ifges(5,2) * t208 + t145 + t199) * t207 / 0.2e1 - (-Ifges(6,2) * t260 + t338 + t93) * t287 / 0.2e1 - (Ifges(6,1) * t287 - t335 + t59) * t260 / 0.2e1 + t244 * t13 / 0.2e1 + Ifges(5,5) * t160 - Ifges(5,6) * t161 - Ifges(6,6) * t110 + Ifges(6,5) * t111 - t52 * mrSges(5,2) + t53 * mrSges(5,1) - t208 * (Ifges(5,1) * t207 - t332) / 0.2e1 - t9 * mrSges(6,2) - t61 * t320 / 0.2e1 + t60 * t321 / 0.2e1 + t389 * qJD(6) + (-mrSges(6,1) - t283) * t8 + t118 * t349 + t117 * t350 + mrSges(7,3) * t352 + t12 * t357 + t144 * t360 + t276 * t367 + t278 * t371 + t280 * t372 - t287 * t264;
t252 = qJD(2) ^ 2;
t234 = -t298 - pkin(5);
t225 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t308;
t224 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t310;
t212 = (mrSges(4,1) * t246 + mrSges(4,2) * t250) * t302;
t193 = pkin(10) + t197;
t192 = -pkin(5) - t196;
t188 = mrSges(5,1) * t240 - t349;
t187 = -mrSges(5,2) * t240 + t350;
t184 = t238 + t355;
t177 = qJD(3) * t203 + t250 * t293;
t176 = -qJD(3) * t204 - t246 * t293;
t138 = -mrSges(6,2) * t240 + t382;
t114 = mrSges(5,1) * t161 + mrSges(5,2) * t160;
t107 = Ifges(7,3) * t110;
t105 = -mrSges(6,1) * t287 + mrSges(6,2) * t260;
t76 = t238 + t80;
t74 = -qJD(4) * t149 + t176 * t249 - t177 * t245;
t73 = qJD(4) * t148 + t176 * t245 + t177 * t249;
t46 = t326 * t98 - t329;
t45 = t241 * t98 + t87;
t31 = -mrSges(7,2) * t110 + mrSges(7,3) * t65;
t30 = mrSges(7,1) * t110 - mrSges(7,3) * t64;
t27 = t241 * t74 + t326 * t73;
t26 = t241 * t73 - t326 * t74;
t22 = -mrSges(7,1) * t65 + mrSges(7,2) * t64;
t21 = t244 * t76 + t248 * t49;
t20 = -t244 * t49 + t248 * t76;
t19 = t244 * t80 + t248 * t46;
t18 = -t244 * t46 + t248 * t80;
t15 = qJD(6) * t267 - t244 * t27 + t248 * t294;
t14 = qJD(6) * t78 + t244 * t294 + t248 * t27;
t1 = [t27 * t138 + t14 * t81 + t15 * t82 + t176 * t224 + t177 * t225 + t73 * t187 + t74 * t188 + t95 * t22 + t78 * t30 - t267 * t31 - t328 * t26 + (-t110 * t96 + t111 * t95) * mrSges(6,3) + (-t148 * t160 - t149 * t161) * mrSges(5,3) + (-t203 * t250 - t204 * t246) * mrSges(4,3) * t302 + ((-mrSges(3,2) * t252 - t114 - t212 - t54) * t251 + (-mrSges(3,1) * t252 + (t105 + t164 + qJD(2) * (-mrSges(4,1) * t250 + mrSges(4,2) * t246)) * qJD(2)) * t247) * t242 + m(7) * (t14 * t17 + t15 * t16 - t2 * t267 + t26 * t36 + t3 * t78 + t368) + m(6) * (-t26 * t41 + t27 * t42 + t368 + t9 * t96 + (-t130 * t251 + t162 * t309) * t242) + m(5) * (t117 * t74 + t118 * t73 + t148 * t53 + t149 * t52 + (t200 * t309 - t201 * t251) * t242) + m(4) * (t155 * t204 + t156 * t203 + t176 * t185 + t177 * t186 + (t223 - t296) * t294); m(4) * ((-t185 * t250 - t186 * t246) * qJD(3) + t271) * pkin(8) + (-t110 * t90 + t111 * t89 - t115 * t42 - t116 * t41) * mrSges(6,3) + (-0.3e1 / 0.2e1 * t246 ^ 2 + 0.3e1 / 0.2e1 * t250 ^ 2) * Ifges(4,4) * t302 + (t130 * t189 + t154 * t162 + t381 * t41 - t391 * t42 + t9 * t90 + t369) * m(6) - t391 * t138 + t240 * (Ifges(5,5) * t174 - Ifges(5,6) * t175) / 0.2e1 + t200 * (mrSges(5,1) * t175 + mrSges(5,2) * t174) + t236 * t114 + t201 * (-mrSges(5,1) * t214 + mrSges(5,2) * t215) - pkin(2) * t212 + t189 * t54 + t174 * t145 / 0.2e1 - t175 * t144 / 0.2e1 + t154 * t105 + t89 * t22 + ((t224 * t246 - t225 * t250) * t251 + (-m(6) * t162 - t105 + t377) * t247 + ((t185 * t246 - t186 * t250) * t251 + (-t344 - t223) * t247) * m(4)) * t314 + (t160 * t215 + t174 * t360) * Ifges(5,1) + (-t117 * t174 - t118 * t175 - t160 * t182 - t161 * t183 + t214 * t52 - t215 * t53) * mrSges(5,3) + (-t161 * t214 - t175 * t361) * Ifges(5,2) + (t160 * t214 - t161 * t215 + t174 * t361 - t175 * t360) * Ifges(5,4) + t39 * t30 + t40 * t31 + (t279 * t371 + t277 * t367 + t281 * t372 + t130 * mrSges(6,2) - Ifges(6,4) * t110 + Ifges(6,1) * t111 + t12 * t358 + t13 * t357 + (mrSges(6,3) + t282) * t8 + t284 * mrSges(7,3) + (t276 * t363 + t278 * t366 + t280 * t365 + t36 * t283 + t61 * t358 - t248 * t60 / 0.2e1 - t374 * mrSges(7,3)) * qJD(6)) * t168 + (t343 / 0.2e1 + t341 / 0.2e1 + t16 * mrSges(7,1) - t17 * mrSges(7,2) + t340 / 0.2e1 - t330 / 0.2e1 + t162 * mrSges(6,1) - t337 / 0.2e1 - t335 / 0.2e1 + t59 / 0.2e1 - t92 / 0.2e1) * t115 + t383 * t81 + t384 * t82 + (t16 * t384 + t17 * t383 + t2 * t40 + t3 * t39 - t36 * t381 + t369) * m(7) + (t331 / 0.2e1 + t162 * mrSges(6,2) + t338 / 0.2e1 + t336 / 0.2e1 + t93 / 0.2e1 + t275 * mrSges(7,3) + t389) * t116 + t271 * mrSges(4,3) + (t130 * mrSges(6,1) - Ifges(6,4) * t111 + t107 / 0.2e1 - t9 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2)) * t110 + t373) * t167 + t375 * t187 + t376 * t188 + (t117 * t376 + t118 * t375 + t182 * t53 + t183 * t52 + t200 * t239 + t201 * t236) * m(5) + t381 * t328 + ((-pkin(8) * t224 - t263 + t292) * t250 + (-pkin(8) * t225 + pkin(3) * t164 + t291 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2)) * t308 - t262) * t246) * qJD(3); t253 + t378 * t138 + t186 * t224 - t185 * t225 + t192 * t22 - t120 * t187 - t119 * t188 - t184 * t105 - t155 * mrSges(4,2) + t156 * mrSges(4,1) - t21 * t81 - t20 * t82 + ((t292 + t387 + t263) * t250 + (t291 + (t348 / 0.2e1 + (t388 + Ifges(4,2) / 0.2e1) * t250) * qJD(2) + t377 * pkin(3) + t262) * t246) * qJD(2) + (t193 * t31 + t195 * t81 + (-t193 * t82 - t353) * qJD(6)) * t248 + (-t110 * t197 - t111 * t196 + t339) * mrSges(6,3) + (-t195 * t82 + (-qJD(6) * t81 - t30) * t193 + (-t3 - t322) * mrSges(7,3)) * t244 + ((t187 * t249 - t188 * t245) * qJD(4) + (-t160 * t249 - t161 * t245) * mrSges(5,3)) * pkin(3) - t328 * t379 + (-t162 * t184 - t196 * t8 + t197 * t9 + t378 * t42 - t379 * t41) * m(6) + (-t117 * t119 - t118 * t120 + (t245 * t52 + t249 * t53 + (-t117 * t245 + t118 * t249) * qJD(4)) * pkin(3)) * m(5) + (-t16 * t20 - t17 * t21 + t192 * t8 + t256 * t193 + t195 * t374 + t36 * t379) * m(7); -t105 * t355 - t117 * t187 + t118 * t188 - t46 * t138 - t18 * t82 - t19 * t81 + t234 * t22 - t303 * t353 + t253 + t328 * t45 + (-t17 * t304 - t351) * mrSges(7,3) + (-t16 * t18 - t17 * t19 + t234 * t8 - t36 * t45) * m(7) + ((t241 * t9 - t326 * t8) * pkin(4) - t162 * t355 + t41 * t45 - t42 * t46) * m(6) + (-t110 * t354 - t111 * t298 + t339) * mrSges(6,3) + (m(7) * t256 - t244 * t30 + t248 * t31 - t303 * t82 - t304 * t81) * (pkin(10) + t354); t244 * t31 + t248 * t30 + t328 * t260 + t273 * qJD(6) + (-t138 - t273) * t287 + t54 + (t143 * t374 - t260 * t36 - t284) * m(7) + (t260 * t41 - t287 * t42 + t130) * m(6); t107 - t36 * (mrSges(7,1) * t135 + mrSges(7,2) * t134) + (Ifges(7,1) * t134 - t342) * t365 + t60 * t364 + (Ifges(7,5) * t134 - Ifges(7,6) * t135) * t363 - t16 * t81 + t17 * t82 + (t134 * t16 + t135 * t17) * mrSges(7,3) + (-Ifges(7,2) * t135 + t133 + t61) * t366 + t373;];
tauc  = t1(:);
