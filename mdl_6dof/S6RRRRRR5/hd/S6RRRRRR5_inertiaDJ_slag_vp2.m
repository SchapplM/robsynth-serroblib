% Calculate time derivative of joint inertia matrix for
% S6RRRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:56:44
% EndTime: 2019-03-10 03:57:11
% DurationCPUTime: 12.02s
% Computational Cost: add. (26400->695), mult. (66206->1032), div. (0->0), fcn. (67966->12), ass. (0->299)
t405 = qJD(3) + qJD(4);
t278 = sin(qJ(6));
t274 = t278 ^ 2;
t283 = cos(qJ(6));
t275 = t283 ^ 2;
t412 = t274 + t275;
t280 = sin(qJ(4));
t285 = cos(qJ(4));
t281 = sin(qJ(3));
t392 = -pkin(10) - pkin(9);
t326 = t281 * t392;
t286 = cos(qJ(3));
t408 = t392 * t286;
t206 = t280 * t408 + t285 * t326;
t239 = -t280 * t281 + t285 * t286;
t197 = t405 * t239;
t240 = t280 * t286 + t281 * t285;
t198 = t405 * t240;
t279 = sin(qJ(5));
t284 = cos(qJ(5));
t304 = t284 * t239 - t240 * t279;
t120 = qJD(5) * t304 + t197 * t284 - t198 * t279;
t191 = t239 * t279 + t240 * t284;
t332 = qJD(6) * t283;
t301 = t278 * t120 + t191 * t332;
t333 = qJD(6) * t278;
t345 = t283 * t120;
t300 = t191 * t333 - t345;
t409 = -t280 * t326 + t285 * t408;
t181 = pkin(11) * t239 - t409;
t161 = t405 * t409;
t294 = -t197 * pkin(11) + t161;
t335 = qJD(5) * t279;
t160 = t405 * t206;
t297 = -t240 * pkin(11) + t206;
t407 = -t198 * pkin(11) + qJD(5) * t297 + t160;
t66 = -t181 * t335 + t279 * t294 + t284 * t407;
t121 = qJD(5) * t191 + t197 * t279 + t284 * t198;
t339 = qJD(3) * t281;
t273 = pkin(3) * t339;
t184 = pkin(4) * t198 + t273;
t71 = pkin(5) * t121 - pkin(12) * t120 + t184;
t270 = -pkin(3) * t286 - pkin(2);
t215 = -pkin(4) * t239 + t270;
t134 = -pkin(5) * t304 - pkin(12) * t191 + t215;
t136 = t284 * t181 + t279 * t297;
t88 = t134 * t283 - t136 * t278;
t20 = qJD(6) * t88 + t278 * t71 + t283 * t66;
t89 = t134 * t278 + t136 * t283;
t21 = -qJD(6) * t89 - t278 * t66 + t283 * t71;
t411 = t20 * t283 - t21 * t278;
t277 = cos(pkin(6));
t276 = sin(pkin(6));
t282 = sin(qJ(2));
t351 = t276 * t282;
t227 = t277 * t286 - t281 * t351;
t228 = t277 * t281 + t286 * t351;
t182 = t227 * t285 - t228 * t280;
t287 = cos(qJ(2));
t341 = qJD(2) * t276;
t321 = t287 * t341;
t202 = -qJD(3) * t228 - t281 * t321;
t203 = qJD(3) * t227 + t286 * t321;
t114 = qJD(4) * t182 + t202 * t280 + t203 * t285;
t183 = t227 * t280 + t228 * t285;
t115 = -qJD(4) * t183 + t202 * t285 - t203 * t280;
t305 = t284 * t182 - t183 * t279;
t56 = qJD(5) * t305 + t114 * t284 + t115 * t279;
t138 = t182 * t279 + t183 * t284;
t57 = qJD(5) * t138 + t114 * t279 - t284 * t115;
t350 = t276 * t287;
t235 = t277 * t282 * pkin(1) + pkin(8) * t350;
t222 = t235 * qJD(2);
t166 = -t202 * pkin(3) + t222;
t91 = -t115 * pkin(4) + t166;
t16 = t57 * pkin(5) - t56 * pkin(12) + t91;
t217 = pkin(9) * t277 + t235;
t218 = (-pkin(2) * t287 - pkin(9) * t282 - pkin(1)) * t276;
t170 = -t281 * t217 + t286 * t218;
t147 = -pkin(3) * t350 - t228 * pkin(10) + t170;
t171 = t286 * t217 + t281 * t218;
t157 = pkin(10) * t227 + t171;
t92 = t285 * t147 - t280 * t157;
t81 = -pkin(4) * t350 - t183 * pkin(11) + t92;
t93 = t280 * t147 + t285 * t157;
t85 = pkin(11) * t182 + t93;
t375 = t279 * t81 + t284 * t85;
t47 = -pkin(12) * t350 + t375;
t263 = pkin(8) * t351;
t379 = pkin(1) * t287;
t216 = t263 + (-pkin(2) - t379) * t277;
t186 = -t227 * pkin(3) + t216;
t140 = -t182 * pkin(4) + t186;
t68 = -pkin(5) * t305 - t138 * pkin(12) + t140;
t22 = -t278 * t47 + t283 * t68;
t340 = qJD(2) * t282;
t322 = t276 * t340;
t220 = (pkin(2) * t282 - pkin(9) * t287) * t341;
t234 = t277 * t379 - t263;
t221 = t234 * qJD(2);
t130 = -qJD(3) * t171 + t286 * t220 - t221 * t281;
t96 = pkin(3) * t322 - pkin(10) * t203 + t130;
t338 = qJD(3) * t286;
t129 = -t217 * t339 + t218 * t338 + t281 * t220 + t286 * t221;
t98 = pkin(10) * t202 + t129;
t39 = -qJD(4) * t93 - t280 * t98 + t285 * t96;
t28 = pkin(4) * t322 - pkin(11) * t114 + t39;
t336 = qJD(4) * t285;
t337 = qJD(4) * t280;
t38 = t147 * t336 - t157 * t337 + t280 * t96 + t285 * t98;
t30 = pkin(11) * t115 + t38;
t334 = qJD(5) * t284;
t8 = t279 * t28 + t284 * t30 + t81 * t334 - t335 * t85;
t5 = pkin(12) * t322 + t8;
t2 = qJD(6) * t22 + t16 * t278 + t283 * t5;
t23 = t278 * t68 + t283 * t47;
t3 = -qJD(6) * t23 + t16 * t283 - t278 * t5;
t410 = t2 * t283 - t3 * t278;
t382 = t278 / 0.2e1;
t381 = t283 / 0.2e1;
t318 = -t333 / 0.2e1;
t406 = t412 * t284;
t404 = Ifges(5,5) * t114 + Ifges(5,6) * t115 + Ifges(5,3) * t322;
t403 = Ifges(4,5) * t203 + Ifges(4,6) * t202 + Ifges(4,3) * t322;
t9 = -qJD(5) * t375 - t279 * t30 + t28 * t284;
t311 = mrSges(7,1) * t278 + mrSges(7,2) * t283;
t243 = t311 * qJD(6);
t268 = -pkin(4) * t284 - pkin(5);
t219 = t268 * t243;
t250 = -mrSges(7,1) * t283 + mrSges(7,2) * t278;
t236 = pkin(4) * t250 * t335;
t328 = pkin(4) * t334;
t315 = mrSges(7,3) * t328;
t256 = t274 * t315;
t257 = t275 * t315;
t402 = t219 + t236 + t256 + t257;
t401 = 2 * m(5);
t400 = 2 * m(6);
t399 = 2 * m(7);
t398 = -2 * mrSges(3,3);
t397 = -2 * mrSges(6,3);
t67 = t181 * t334 + t279 * t407 - t284 * t294;
t396 = 0.2e1 * t67;
t135 = t279 * t181 - t284 * t297;
t395 = 0.2e1 * t135;
t394 = 0.2e1 * t222;
t302 = -t283 * t138 + t278 * t350;
t36 = qJD(6) * t302 - t278 * t56 + t283 * t322;
t393 = t36 / 0.2e1;
t124 = -t278 * t138 - t283 * t350;
t391 = t124 / 0.2e1;
t390 = t239 / 0.2e1;
t389 = t240 / 0.2e1;
t271 = Ifges(7,5) * t332;
t388 = Ifges(7,6) * t318 + t271 / 0.2e1;
t372 = Ifges(7,4) * t278;
t310 = Ifges(7,1) * t283 - t372;
t248 = t310 * qJD(6);
t387 = t248 / 0.2e1;
t386 = Ifges(7,5) * t382 + Ifges(7,6) * t381;
t371 = Ifges(7,4) * t283;
t254 = Ifges(7,1) * t278 + t371;
t384 = t254 / 0.2e1;
t383 = -t278 / 0.2e1;
t378 = pkin(5) * t243;
t374 = Ifges(4,4) * t281;
t373 = Ifges(4,4) * t286;
t370 = Ifges(7,6) * t278;
t369 = pkin(4) * qJD(5);
t368 = t135 * t67;
t269 = pkin(3) * t285 + pkin(4);
t347 = t279 * t280;
t187 = t269 * t334 + (-t280 * t335 + (t284 * t285 - t347) * qJD(4)) * pkin(3);
t367 = t187 * mrSges(6,2);
t364 = t221 * mrSges(3,2);
t363 = t222 * mrSges(3,1);
t127 = -mrSges(6,1) * t350 - t138 * mrSges(6,3);
t78 = -mrSges(7,1) * t124 - mrSges(7,2) * t302;
t362 = -t127 + t78;
t346 = t280 * t284;
t188 = t269 * t335 + (t280 * t334 + (t279 * t285 + t346) * qJD(4)) * pkin(3);
t361 = t135 * t188;
t360 = t135 * t279;
t359 = t187 * t274;
t358 = t187 * t275;
t357 = t187 * t278;
t356 = t187 * t283;
t355 = t191 * t278;
t354 = t191 * t283;
t348 = t278 * t284;
t344 = t283 * t284;
t343 = Ifges(6,5) * t120 - Ifges(6,6) * t121;
t342 = Ifges(5,5) * t197 - Ifges(5,6) * t198;
t226 = pkin(3) * t346 + t279 * t269;
t331 = 0.2e1 * t276;
t35 = qJD(6) * t124 + t278 * t322 + t283 * t56;
t12 = Ifges(7,5) * t35 + Ifges(7,6) * t36 + Ifges(7,3) * t57;
t329 = Ifges(6,5) * t56 - Ifges(6,6) * t57 + Ifges(6,3) * t322;
t15 = -mrSges(7,1) * t36 + mrSges(7,2) * t35;
t6 = -pkin(5) * t322 - t9;
t327 = m(7) * t6 + t15;
t60 = mrSges(7,1) * t301 - mrSges(7,2) * t300;
t323 = m(7) * t67 + t60;
t317 = t332 / 0.2e1;
t316 = t412 * t187;
t312 = -t279 * mrSges(6,1) - t284 * mrSges(6,2);
t309 = -Ifges(7,2) * t278 + t371;
t48 = -t279 * t85 + t284 * t81;
t145 = mrSges(7,2) * t304 - mrSges(7,3) * t355;
t146 = -mrSges(7,1) * t304 - mrSges(7,3) * t354;
t306 = t283 * t145 - t278 * t146;
t225 = -pkin(3) * t347 + t269 * t284;
t126 = mrSges(6,2) * t350 + mrSges(6,3) * t305;
t82 = mrSges(7,2) * t305 + mrSges(7,3) * t124;
t83 = -mrSges(7,1) * t305 + mrSges(7,3) * t302;
t303 = -t278 * t83 + t283 * t82 + t126;
t246 = t309 * qJD(6);
t252 = Ifges(7,2) * t283 + t372;
t299 = t283 * t246 + t278 * t248 - t252 * t333 + t254 * t332;
t298 = (-mrSges(5,1) * t280 - mrSges(5,2) * t285) * qJD(4) * pkin(3);
t172 = t188 * t250;
t179 = mrSges(7,3) * t359;
t180 = mrSges(7,3) * t358;
t185 = t188 * mrSges(6,1);
t223 = -pkin(5) - t225;
t196 = t223 * t243;
t295 = t172 + t179 + t180 - t185 + t196 + t299 - t367;
t43 = -Ifges(7,5) * t300 - Ifges(7,6) * t301 + Ifges(7,3) * t121;
t17 = mrSges(7,1) * t57 - mrSges(7,3) * t35;
t18 = -mrSges(7,2) * t57 + mrSges(7,3) * t36;
t293 = -t83 * t332 + m(7) * (-t22 * t332 - t23 * t333 + t410) + t283 * t18 - t278 * t17 - t82 * t333;
t69 = mrSges(7,1) * t121 + mrSges(7,3) * t300;
t70 = -mrSges(7,2) * t121 - mrSges(7,3) * t301;
t292 = -t146 * t332 + m(7) * (-t332 * t88 - t333 * t89 + t411) + t283 * t70 - t278 * t69 - t145 * t333;
t13 = Ifges(7,4) * t35 + Ifges(7,2) * t36 + Ifges(7,6) * t57;
t14 = Ifges(7,1) * t35 + Ifges(7,4) * t36 + Ifges(7,5) * t57;
t46 = pkin(5) * t350 - t48;
t62 = -Ifges(7,4) * t302 + Ifges(7,2) * t124 - Ifges(7,6) * t305;
t63 = -Ifges(7,1) * t302 + Ifges(7,4) * t124 - Ifges(7,5) * t305;
t291 = t9 * mrSges(6,1) - t8 * mrSges(6,2) - t302 * t387 + t13 * t381 - t305 * t388 + t14 * t382 + t46 * t243 + t246 * t391 + t6 * t250 + t252 * t393 + t63 * t317 + t62 * t318 + t35 * t384 + t57 * t386 + t329 + ((-t22 * t283 - t23 * t278) * qJD(6) + t410) * mrSges(7,3);
t110 = -Ifges(7,6) * t304 + t191 * t309;
t111 = -Ifges(7,5) * t304 + t191 * t310;
t44 = -Ifges(7,4) * t300 - Ifges(7,2) * t301 + Ifges(7,6) * t121;
t45 = -Ifges(7,1) * t300 - Ifges(7,4) * t301 + Ifges(7,5) * t121;
t290 = t345 * t384 + t111 * t317 + t121 * t386 + t135 * t243 - t246 * t355 / 0.2e1 + t354 * t387 - t304 * t388 + t45 * t382 + t44 * t381 - t66 * mrSges(6,2) + t343 + (t250 - mrSges(6,1)) * t67 - t301 * t252 / 0.2e1 + (t191 * t254 + t110) * t318 + ((-t278 * t89 - t283 * t88) * qJD(6) + t411) * mrSges(7,3);
t289 = t39 * mrSges(5,1) - t38 * mrSges(5,2) + t291 + t404;
t288 = t161 * mrSges(5,1) - t160 * mrSges(5,2) + t290 + t342;
t272 = Ifges(4,5) * t338;
t267 = pkin(4) * t279 + pkin(12);
t262 = Ifges(3,5) * t321;
t255 = Ifges(4,1) * t281 + t373;
t253 = Ifges(4,2) * t286 + t374;
t249 = (Ifges(4,1) * t286 - t374) * qJD(3);
t247 = (-Ifges(4,2) * t281 + t373) * qJD(3);
t244 = (mrSges(4,1) * t281 + mrSges(4,2) * t286) * qJD(3);
t224 = pkin(12) + t226;
t205 = -mrSges(4,1) * t350 - t228 * mrSges(4,3);
t204 = mrSges(4,2) * t350 + t227 * mrSges(4,3);
t201 = Ifges(5,1) * t240 + Ifges(5,4) * t239;
t200 = Ifges(5,4) * t240 + Ifges(5,2) * t239;
t199 = -mrSges(5,1) * t239 + mrSges(5,2) * t240;
t178 = mrSges(4,1) * t322 - mrSges(4,3) * t203;
t177 = -mrSges(4,2) * t322 + mrSges(4,3) * t202;
t174 = Ifges(4,1) * t228 + Ifges(4,4) * t227 - Ifges(4,5) * t350;
t173 = Ifges(4,4) * t228 + Ifges(4,2) * t227 - Ifges(4,6) * t350;
t165 = -mrSges(5,1) * t350 - t183 * mrSges(5,3);
t164 = mrSges(5,2) * t350 + t182 * mrSges(5,3);
t158 = -mrSges(4,1) * t202 + mrSges(4,2) * t203;
t156 = Ifges(5,1) * t197 - Ifges(5,4) * t198;
t155 = Ifges(5,4) * t197 - Ifges(5,2) * t198;
t154 = mrSges(5,1) * t198 + mrSges(5,2) * t197;
t152 = Ifges(6,1) * t191 + Ifges(6,4) * t304;
t151 = Ifges(6,4) * t191 + Ifges(6,2) * t304;
t150 = -mrSges(6,1) * t304 + mrSges(6,2) * t191;
t149 = Ifges(4,1) * t203 + Ifges(4,4) * t202 + Ifges(4,5) * t322;
t148 = Ifges(4,4) * t203 + Ifges(4,2) * t202 + Ifges(4,6) * t322;
t142 = t311 * t191;
t139 = -mrSges(5,1) * t182 + mrSges(5,2) * t183;
t133 = Ifges(5,1) * t183 + Ifges(5,4) * t182 - Ifges(5,5) * t350;
t132 = Ifges(5,4) * t183 + Ifges(5,2) * t182 - Ifges(5,6) * t350;
t109 = -Ifges(7,3) * t304 + (Ifges(7,5) * t283 - t370) * t191;
t100 = -mrSges(5,2) * t322 + mrSges(5,3) * t115;
t99 = mrSges(5,1) * t322 - mrSges(5,3) * t114;
t90 = -mrSges(6,1) * t305 + mrSges(6,2) * t138;
t87 = Ifges(6,1) * t138 + Ifges(6,4) * t305 - Ifges(6,5) * t350;
t86 = Ifges(6,4) * t138 + Ifges(6,2) * t305 - Ifges(6,6) * t350;
t77 = Ifges(6,1) * t120 - Ifges(6,4) * t121;
t76 = Ifges(6,4) * t120 - Ifges(6,2) * t121;
t75 = mrSges(6,1) * t121 + mrSges(6,2) * t120;
t74 = -mrSges(5,1) * t115 + mrSges(5,2) * t114;
t73 = Ifges(5,1) * t114 + Ifges(5,4) * t115 + Ifges(5,5) * t322;
t72 = Ifges(5,4) * t114 + Ifges(5,2) * t115 + Ifges(5,6) * t322;
t61 = -Ifges(7,5) * t302 + Ifges(7,6) * t124 - Ifges(7,3) * t305;
t52 = -mrSges(6,2) * t322 - mrSges(6,3) * t57;
t51 = mrSges(6,1) * t322 - mrSges(6,3) * t56;
t26 = mrSges(6,1) * t57 + mrSges(6,2) * t56;
t25 = Ifges(6,1) * t56 - Ifges(6,4) * t57 + Ifges(6,5) * t322;
t24 = Ifges(6,4) * t56 - Ifges(6,2) * t57 + Ifges(6,6) * t322;
t1 = [(-t86 + t61) * t57 + t228 * t149 + t227 * t148 + t202 * t173 + t203 * t174 + 0.2e1 * t129 * t204 + 0.2e1 * t130 * t205 + 0.2e1 * t216 * t158 + t182 * t72 + t183 * t73 + 0.2e1 * t186 * t74 + 0.2e1 * t166 * t139 + 0.2e1 * t171 * t177 + 0.2e1 * t170 * t178 + 0.2e1 * t38 * t164 + 0.2e1 * t39 * t165 + t138 * t25 + 0.2e1 * t140 * t26 + t115 * t132 + t114 * t133 + t124 * t13 + 0.2e1 * t8 * t126 + 0.2e1 * t9 * t127 + 0.2e1 * t92 * t99 + 0.2e1 * t93 * t100 + 0.2e1 * t2 * t82 + 0.2e1 * t3 * t83 + t56 * t87 + 0.2e1 * t91 * t90 + 0.2e1 * t6 * t78 + t36 * t62 + t35 * t63 + 0.2e1 * t48 * t51 + 0.2e1 * t46 * t15 + 0.2e1 * t23 * t18 + 0.2e1 * t22 * t17 + 0.2e1 * t375 * t52 + (t140 * t91 + t375 * t8 + t48 * t9) * t400 + 0.2e1 * m(3) * (t221 * t235 - t222 * t234) + 0.2e1 * m(4) * (t129 * t171 + t130 * t170 + t216 * t222) - t302 * t14 - (t12 - t24) * t305 + (mrSges(3,3) * t282 * t394 + (0.2e1 * mrSges(3,3) * t221 - t329 - t403 - t404) * t287 + ((t234 * t398 + Ifges(3,5) * t277 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t287) * t331) * t287 + (t235 * t398 + Ifges(4,5) * t228 + Ifges(5,5) * t183 + Ifges(6,5) * t138 - 0.2e1 * Ifges(3,6) * t277 + Ifges(4,6) * t227 + Ifges(5,6) * t182 + Ifges(6,6) * t305 + (-mrSges(3,1) * pkin(1) - Ifges(3,4) * t282) * t331 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3) - Ifges(5,3) - Ifges(6,3)) * t350) * t282) * qJD(2)) * t276 + (t262 - 0.2e1 * t363 - 0.2e1 * t364) * t277 + (-mrSges(4,1) * t227 + mrSges(4,2) * t228) * t394 + (t2 * t23 + t22 * t3 + t46 * t6) * t399 + (t166 * t186 + t38 * t93 + t39 * t92) * t401; (t15 - t51) * t135 + m(7) * (t135 * t6 + t2 * t89 + t20 * t23 + t21 * t22 + t3 * t88 + t46 * t67) + ((t129 * t286 - t130 * t281 + (-t170 * t286 - t171 * t281) * qJD(3)) * pkin(9) - pkin(2) * t222) * m(4) + (-t151 / 0.2e1 + t109 / 0.2e1) * t57 - t363 - t364 + m(5) * (t160 * t93 + t161 * t92 + t166 * t270 + t186 * t273 + t206 * t39 - t38 * t409) - t409 * t100 + t270 * t74 + t216 * t244 + t227 * t247 / 0.2e1 + t228 * t249 / 0.2e1 + t202 * t253 / 0.2e1 + t203 * t255 / 0.2e1 + t206 * t99 + t215 * t26 + t197 * t133 / 0.2e1 - t198 * t132 / 0.2e1 + t166 * t199 + t115 * t200 / 0.2e1 + t114 * t201 / 0.2e1 + t182 * t155 / 0.2e1 + t183 * t156 / 0.2e1 + t184 * t90 + t186 * t154 - pkin(2) * t158 + t160 * t164 + t161 * t165 + t2 * t145 + t3 * t146 + t91 * t150 + t56 * t152 / 0.2e1 + t136 * t52 + t138 * t77 / 0.2e1 + t140 * t75 + t6 * t142 + t66 * t126 + t35 * t111 / 0.2e1 + t20 * t82 + t21 * t83 + t88 * t17 + t89 * t18 + t22 * t69 + t23 * t70 + t46 * t60 + t262 + m(6) * (-t135 * t9 + t136 * t8 + t140 * t184 + t215 * t91 + t375 * t66 - t48 * t67) + (-t86 / 0.2e1 + t61 / 0.2e1 - t375 * mrSges(6,3)) * t121 + ((-pkin(9) * t205 - t170 * mrSges(4,3) + t174 / 0.2e1) * t286 + (-t173 / 0.2e1 - pkin(9) * t204 + pkin(3) * t139 - t171 * mrSges(4,3)) * t281) * qJD(3) + (-t197 * t92 - t198 * t93 + t239 * t38 - t240 * t39) * mrSges(5,3) - t302 * t45 / 0.2e1 - (t43 / 0.2e1 - t76 / 0.2e1) * t305 + ((Ifges(4,5) * t281 / 0.2e1 + Ifges(4,6) * t286 / 0.2e1 - Ifges(3,6) + Ifges(6,5) * t191 / 0.2e1 + Ifges(6,6) * t304 / 0.2e1 + Ifges(5,5) * t389 + Ifges(5,6) * t390) * t340 - (-Ifges(4,6) * t339 + t272 + t342 + t343) * t287 / 0.2e1) * t276 - (t12 / 0.2e1 - t24 / 0.2e1 - t8 * mrSges(6,3)) * t304 + (t87 / 0.2e1 + t63 * t381 + t62 * t383 - t48 * mrSges(6,3)) * t120 + (t25 / 0.2e1 + t14 * t381 + t13 * t383 - t9 * mrSges(6,3) + (-t283 * t62 / 0.2e1 + t63 * t383) * qJD(6)) * t191 + t362 * t67 + t73 * t389 + t72 * t390 + t44 * t391 + t110 * t393 + (t222 * mrSges(4,2) + t149 / 0.2e1 - pkin(9) * t178 - t130 * mrSges(4,3)) * t281 + (-t222 * mrSges(4,1) + t148 / 0.2e1 + pkin(9) * t177 + t129 * mrSges(4,3)) * t286; (t286 * t255 + (0.2e1 * pkin(3) * t199 - t253) * t281) * qJD(3) + (-t160 * t409 + t161 * t206 + t270 * t273) * t401 + 0.2e1 * (t160 * t239 - t161 * t240 - t197 * t206 + t198 * t409) * mrSges(5,3) + t286 * t247 + t281 * t249 + 0.2e1 * t270 * t154 + t240 * t156 - 0.2e1 * pkin(2) * t244 + t239 * t155 + 0.2e1 * t215 * t75 - t198 * t200 + t197 * t201 + 0.2e1 * t184 * t150 + 0.2e1 * t20 * t145 + 0.2e1 * t21 * t146 + 0.2e1 * t88 * t69 + 0.2e1 * t89 * t70 - (t397 * t66 + t43 - t76) * t304 + (mrSges(6,3) * t395 - t110 * t278 + t111 * t283 + t152) * t120 + (mrSges(6,3) * t396 - t278 * t44 + t283 * t45 + t77 + (-t110 * t283 - t111 * t278) * qJD(6)) * t191 + (t136 * t397 + t109 - t151) * t121 + t60 * t395 + t142 * t396 + (t20 * t89 + t21 * t88 + t368) * t399 + (t136 * t66 + t184 * t215 + t368) * t400; t289 + t293 * t224 + (t164 * t336 + m(5) * (t280 * t38 + t285 * t39 + t336 * t93 - t337 * t92) + t285 * t99 + t280 * t100 - t165 * t337) * pkin(3) + t223 * t15 + t225 * t51 + t226 * t52 + t362 * t188 + t303 * t187 - t129 * mrSges(4,2) + t130 * mrSges(4,1) + m(7) * (t188 * t46 - t22 * t357 + t223 * t6 + t23 * t356) + m(6) * (t187 * t375 - t188 * t48 + t225 * t9 + t226 * t8) + t403; t288 + (-t120 * t225 - t121 * t226 + t187 * t304 + t188 * t191) * mrSges(6,3) + (m(5) * (t160 * t280 + t161 * t285 + (-t206 * t280 - t285 * t409) * qJD(4)) + (-t285 * t197 - t280 * t198 + (t239 * t285 + t240 * t280) * qJD(4)) * mrSges(5,3)) * pkin(3) + t292 * t224 + t223 * t60 + t188 * t142 + (-Ifges(4,6) * t281 + (-mrSges(4,1) * t286 + mrSges(4,2) * t281) * pkin(9)) * qJD(3) + t306 * t187 + m(7) * (t223 * t67 + t356 * t89 - t357 * t88 + t361) + t272 + m(6) * (t136 * t187 - t225 * t67 + t226 * t66 + t361); -0.2e1 * t367 + 0.2e1 * t172 + 0.2e1 * t179 + 0.2e1 * t180 - 0.2e1 * t185 + 0.2e1 * t196 + 0.2e1 * t298 + (t188 * t223 + t224 * t316) * t399 + (t187 * t226 - t188 * t225) * t400 + t299; t289 + (m(6) * (t279 * t8 + t284 * t9) + t284 * t51 + t279 * t52 + (t362 * t279 + t303 * t284 + m(7) * (-t22 * t348 + t23 * t344 + t279 * t46) + m(6) * (-t279 * t48 + t284 * t375)) * qJD(5)) * pkin(4) + t327 * t268 + t293 * t267; t288 + (m(6) * (t279 * t66 - t284 * t67) + (-t120 * t284 - t121 * t279) * mrSges(6,3) + ((t191 * mrSges(6,3) + t142) * t279 + (mrSges(6,3) * t304 + t306) * t284 + m(7) * (t344 * t89 - t348 * t88 + t360) + m(6) * (t136 * t284 + t360)) * qJD(5)) * pkin(4) + t292 * t267 + t323 * t268; t295 + m(7) * (t188 * t268 + (t358 + t359) * t267) + (m(6) * (t187 * t279 - t188 * t284) + (m(7) * (t223 * t279 + t224 * t406) + m(6) * (-t225 * t279 + t226 * t284) + t312) * qJD(5)) * pkin(4) + t298 + t402; 0.2e1 * t219 + 0.2e1 * t236 + 0.2e1 * t256 + 0.2e1 * t257 + 0.2e1 * (m(7) * (t267 * t406 + t268 * t279) + t312) * t369 + t299; -pkin(5) * t327 + pkin(12) * t293 + t291; -pkin(5) * t323 + pkin(12) * t292 + t290; -t378 + m(7) * (-pkin(5) * t188 + pkin(12) * t316) + t295; -t378 + (m(7) * (-pkin(5) * t279 + pkin(12) * t406) + t312) * t369 + t299 + t402; t299 - 0.2e1 * t378; mrSges(7,1) * t3 - mrSges(7,2) * t2 + t12; mrSges(7,1) * t21 - mrSges(7,2) * t20 + t43; t271 - t311 * t187 + (t224 * t250 - t370) * qJD(6); t271 - t311 * t328 + (t250 * t267 - t370) * qJD(6); t271 + (pkin(12) * t250 - t370) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
