% Calculate time derivative of joint inertia matrix for
% S6RRRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 02:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP10_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP10_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP10_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP10_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP10_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:16:22
% EndTime: 2019-03-10 02:16:42
% DurationCPUTime: 8.74s
% Computational Cost: add. (13811->815), mult. (36338->1146), div. (0->0), fcn. (34573->10), ass. (0->309)
t406 = -mrSges(6,1) - mrSges(7,1);
t294 = sin(qJ(3));
t293 = sin(qJ(4));
t298 = cos(qJ(3));
t356 = qJD(3) * t298;
t329 = t293 * t356;
t297 = cos(qJ(4));
t353 = qJD(4) * t297;
t305 = t294 * t353 + t329;
t354 = qJD(4) * t294;
t306 = -t293 * t354 + t297 * t356;
t254 = -pkin(3) * t298 - pkin(10) * t294 - pkin(2);
t362 = t297 * t298;
t276 = pkin(9) * t362;
t216 = t293 * t254 + t276;
t405 = qJD(4) * t216;
t403 = t293 / 0.2e1;
t383 = t297 / 0.2e1;
t291 = cos(pkin(6));
t290 = sin(pkin(6));
t295 = sin(qJ(2));
t366 = t290 * t295;
t229 = t291 * t294 + t298 * t366;
t299 = cos(qJ(2));
t359 = qJD(2) * t290;
t334 = t299 * t359;
t190 = qJD(3) * t229 + t294 * t334;
t292 = sin(qJ(5));
t296 = cos(qJ(5));
t365 = t290 * t299;
t192 = -t229 * t293 - t297 * t365;
t307 = -t229 * t297 + t293 * t365;
t308 = t296 * t192 + t292 * t307;
t228 = -t291 * t298 + t294 * t366;
t191 = -qJD(3) * t228 + t298 * t334;
t358 = qJD(2) * t295;
t335 = t290 * t358;
t96 = qJD(4) * t307 - t191 * t293 + t297 * t335;
t97 = qJD(4) * t192 + t191 * t297 + t293 * t335;
t35 = qJD(5) * t308 + t292 * t96 + t296 * t97;
t112 = t192 * t292 - t296 * t307;
t36 = qJD(5) * t112 + t292 * t97 - t296 * t96;
t10 = Ifges(7,4) * t35 + Ifges(7,2) * t190 + Ifges(7,6) * t36;
t9 = Ifges(6,5) * t35 - Ifges(6,6) * t36 + Ifges(6,3) * t190;
t402 = t10 + t9;
t240 = t292 * t297 + t293 * t296;
t226 = t240 * t294;
t239 = t292 * t293 - t296 * t297;
t398 = qJD(4) + qJD(5);
t129 = -t226 * t398 - t239 * t356;
t352 = qJD(5) * t292;
t363 = t294 * t297;
t364 = t293 * t294;
t130 = -t352 * t364 + (t363 * t398 + t329) * t296 + t306 * t292;
t357 = qJD(3) * t294;
t64 = Ifges(6,5) * t129 - Ifges(6,6) * t130 + Ifges(6,3) * t357;
t65 = Ifges(7,4) * t129 + Ifges(7,2) * t357 + Ifges(7,6) * t130;
t401 = -t64 - t65;
t271 = pkin(8) * t366;
t380 = pkin(1) * t299;
t217 = t271 + (-pkin(2) - t380) * t291;
t139 = t228 * pkin(3) - t229 * pkin(10) + t217;
t234 = t291 * t295 * pkin(1) + pkin(8) * t365;
t218 = pkin(9) * t291 + t234;
t219 = (-pkin(2) * t299 - pkin(9) * t295 - pkin(1)) * t290;
t148 = t298 * t218 + t294 * t219;
t141 = -pkin(10) * t365 + t148;
t72 = t293 * t139 + t297 * t141;
t238 = t297 * t254;
t379 = pkin(9) * t293;
t169 = -pkin(11) * t363 + t238 + (-pkin(4) - t379) * t298;
t194 = -pkin(11) * t364 + t216;
t400 = t292 * t169 + t296 * t194;
t257 = -mrSges(5,1) * t297 + mrSges(5,2) * t293;
t399 = -m(5) * pkin(3) + t257;
t233 = t291 * t380 - t271;
t252 = (pkin(3) * t294 - pkin(10) * t298) * qJD(3);
t360 = t297 * t252 + t357 * t379;
t103 = (pkin(4) * t294 - pkin(11) * t362) * qJD(3) + (-t276 + (pkin(11) * t294 - t254) * t293) * qJD(4) + t360;
t355 = qJD(4) * t293;
t145 = t293 * t252 + t254 * t353 + (-t297 * t357 - t298 * t355) * pkin(9);
t121 = -pkin(11) * t305 + t145;
t43 = -qJD(5) * t400 + t103 * t296 - t121 * t292;
t223 = (pkin(2) * t295 - pkin(9) * t299) * t359;
t224 = t233 * qJD(2);
t78 = -t218 * t357 + t219 * t356 + t294 * t223 + t298 * t224;
t76 = pkin(10) * t335 + t78;
t225 = t234 * qJD(2);
t92 = t190 * pkin(3) - t191 * pkin(10) + t225;
t25 = -t72 * qJD(4) - t293 * t76 + t297 * t92;
t17 = pkin(4) * t190 - pkin(11) * t97 + t25;
t24 = t139 * t353 - t141 * t355 + t293 * t92 + t297 * t76;
t19 = pkin(11) * t96 + t24;
t71 = t297 * t139 - t141 * t293;
t51 = pkin(4) * t228 + pkin(11) * t307 + t71;
t60 = pkin(11) * t192 + t72;
t377 = t292 * t51 + t296 * t60;
t6 = -qJD(5) * t377 + t17 * t296 - t19 * t292;
t397 = 0.2e1 * m(5);
t396 = 2 * m(6);
t395 = 2 * m(7);
t394 = 0.2e1 * pkin(9);
t393 = -2 * mrSges(3,3);
t391 = t96 / 0.2e1;
t390 = -pkin(11) - pkin(10);
t389 = t192 / 0.2e1;
t388 = -t307 / 0.2e1;
t373 = Ifges(5,4) * t293;
t315 = Ifges(5,1) * t297 - t373;
t222 = -Ifges(5,5) * t298 + t294 * t315;
t387 = t222 / 0.2e1;
t372 = Ifges(5,4) * t297;
t261 = Ifges(5,1) * t293 + t372;
t386 = t261 / 0.2e1;
t385 = -t293 / 0.2e1;
t384 = -t297 / 0.2e1;
t382 = m(5) * t298;
t381 = m(7) * t292;
t378 = pkin(9) * t298;
t288 = t294 * pkin(9);
t82 = mrSges(6,1) * t228 - mrSges(6,3) * t112;
t83 = -mrSges(7,1) * t228 + mrSges(7,2) * t112;
t376 = -t82 + t83;
t375 = Ifges(4,4) * t294;
t374 = Ifges(4,4) * t298;
t371 = Ifges(5,6) * t293;
t370 = t224 * mrSges(3,2);
t369 = t225 * mrSges(3,1);
t368 = t225 * mrSges(4,1);
t367 = t225 * mrSges(4,2);
t180 = t398 * t239;
t181 = t398 * t240;
t116 = -Ifges(6,5) * t180 - Ifges(6,6) * t181;
t117 = -Ifges(7,4) * t180 + Ifges(7,6) * t181;
t227 = t239 * t294;
t201 = -mrSges(6,1) * t298 + mrSges(6,3) * t227;
t202 = mrSges(7,1) * t298 - mrSges(7,2) * t227;
t361 = -t201 + t202;
t253 = pkin(4) * t364 + t288;
t351 = qJD(5) * t296;
t11 = Ifges(6,4) * t35 - Ifges(6,2) * t36 + Ifges(6,6) * t190;
t8 = Ifges(7,5) * t35 + Ifges(7,6) * t190 + Ifges(7,3) * t36;
t349 = t8 / 0.2e1 - t11 / 0.2e1;
t44 = Ifges(5,5) * t97 + Ifges(5,6) * t96 + Ifges(5,3) * t190;
t348 = pkin(4) * t355;
t347 = pkin(4) * t352;
t346 = pkin(4) * t351;
t345 = Ifges(4,6) * t365;
t12 = Ifges(7,1) * t35 + Ifges(7,4) * t190 + Ifges(7,5) * t36;
t13 = Ifges(6,1) * t35 - Ifges(6,4) * t36 + Ifges(6,5) * t190;
t344 = t12 / 0.2e1 + t13 / 0.2e1;
t53 = Ifges(7,5) * t112 + Ifges(7,6) * t228 - Ifges(7,3) * t308;
t56 = Ifges(6,4) * t112 + Ifges(6,2) * t308 + Ifges(6,6) * t228;
t343 = t53 / 0.2e1 - t56 / 0.2e1;
t57 = Ifges(7,1) * t112 + Ifges(7,4) * t228 - Ifges(7,5) * t308;
t58 = Ifges(6,1) * t112 + Ifges(6,4) * t308 + Ifges(6,5) * t228;
t342 = t57 / 0.2e1 + t58 / 0.2e1;
t63 = Ifges(7,5) * t129 + Ifges(7,6) * t357 + Ifges(7,3) * t130;
t66 = Ifges(6,4) * t129 - Ifges(6,2) * t130 + Ifges(6,6) * t357;
t341 = t63 / 0.2e1 - t66 / 0.2e1;
t67 = Ifges(7,1) * t129 + Ifges(7,4) * t357 + Ifges(7,5) * t130;
t68 = Ifges(6,1) * t129 - Ifges(6,4) * t130 + Ifges(6,5) * t357;
t340 = t67 / 0.2e1 + t68 / 0.2e1;
t339 = t390 * t293;
t91 = -Ifges(5,1) * t307 + Ifges(5,4) * t192 + Ifges(5,5) * t228;
t338 = t91 * t383;
t337 = Ifges(4,5) * t191 - Ifges(4,6) * t190 + Ifges(4,3) * t335;
t211 = pkin(4) * t305 + pkin(9) * t356;
t281 = -pkin(4) * t297 - pkin(3);
t336 = qJD(4) * t390;
t263 = t390 * t297;
t197 = -t292 * t263 - t296 * t339;
t330 = t197 * t352;
t115 = -Ifges(7,5) * t180 + Ifges(7,3) * t181;
t118 = -Ifges(6,4) * t180 - Ifges(6,2) * t181;
t328 = t115 / 0.2e1 - t118 / 0.2e1;
t119 = -Ifges(7,1) * t180 + Ifges(7,5) * t181;
t120 = -Ifges(6,1) * t180 - Ifges(6,4) * t181;
t327 = t119 / 0.2e1 + t120 / 0.2e1;
t149 = -Ifges(7,5) * t227 - Ifges(7,6) * t298 + Ifges(7,3) * t226;
t152 = -Ifges(6,4) * t227 - Ifges(6,2) * t226 - Ifges(6,6) * t298;
t326 = t149 / 0.2e1 - t152 / 0.2e1;
t153 = -Ifges(7,1) * t227 - Ifges(7,4) * t298 + Ifges(7,5) * t226;
t154 = -Ifges(6,1) * t227 - Ifges(6,4) * t226 - Ifges(6,5) * t298;
t325 = t153 / 0.2e1 + t154 / 0.2e1;
t184 = Ifges(7,5) * t240 + Ifges(7,3) * t239;
t187 = Ifges(6,4) * t240 - Ifges(6,2) * t239;
t324 = t184 / 0.2e1 - t187 / 0.2e1;
t188 = Ifges(7,1) * t240 + Ifges(7,5) * t239;
t189 = Ifges(6,1) * t240 - Ifges(6,4) * t239;
t323 = t189 / 0.2e1 + t188 / 0.2e1;
t27 = -t190 * mrSges(7,1) + t35 * mrSges(7,2);
t251 = t293 * t336;
t320 = t297 * t336;
t134 = -qJD(5) * t197 + t296 * t251 + t292 * t320;
t198 = -t296 * t263 + t292 * t339;
t135 = qJD(5) * t198 + t292 * t251 - t296 * t320;
t322 = t198 * t134 + t135 * t197;
t147 = -t294 * t218 + t219 * t298;
t215 = -t293 * t378 + t238;
t321 = -qJD(4) * t215 + t145;
t319 = t335 / 0.2e1;
t106 = -mrSges(7,1) * t357 + t129 * mrSges(7,2);
t140 = pkin(3) * t365 - t147;
t285 = Ifges(5,5) * t353;
t318 = -Ifges(5,6) * t355 / 0.2e1 + t285 / 0.2e1 + t116 / 0.2e1 + t117 / 0.2e1;
t317 = Ifges(5,5) * t403 + Ifges(5,6) * t383 + (Ifges(6,5) + Ifges(7,4)) * t240 / 0.2e1 + (-Ifges(6,6) / 0.2e1 + Ifges(7,6) / 0.2e1) * t239;
t316 = mrSges(5,1) * t293 + mrSges(5,2) * t297;
t314 = -Ifges(5,2) * t293 + t372;
t259 = Ifges(5,2) * t297 + t373;
t312 = t24 * t297 - t25 * t293;
t22 = -t292 * t60 + t296 * t51;
t108 = t169 * t296 - t194 * t292;
t79 = -t218 * t356 - t219 * t357 + t223 * t298 - t294 * t224;
t5 = t292 * t17 + t296 * t19 + t51 * t351 - t352 * t60;
t42 = t292 * t103 + t296 * t121 + t169 * t351 - t194 * t352;
t86 = -pkin(4) * t192 + t140;
t77 = -pkin(3) * t335 - t79;
t304 = t116 + t117 + t406 * t135 + (-mrSges(6,2) + mrSges(7,3)) * t134;
t2 = qJ(6) * t190 + qJD(6) * t228 + t5;
t3 = -pkin(5) * t190 - t6;
t303 = t6 * mrSges(6,1) - t3 * mrSges(7,1) - t5 * mrSges(6,2) + t2 * mrSges(7,3) + t402;
t47 = -pkin(4) * t96 + t77;
t38 = qJ(6) * t357 - qJD(6) * t298 + t42;
t40 = -pkin(5) * t357 - t43;
t302 = t43 * mrSges(6,1) - t40 * mrSges(7,1) - t42 * mrSges(6,2) + t38 * mrSges(7,3) - t401;
t270 = qJD(6) + t346;
t300 = -mrSges(6,2) * t346 + t270 * mrSges(7,3) + t347 * t406;
t159 = Ifges(5,5) * t306 - Ifges(5,6) * t305 + Ifges(5,3) * t357;
t289 = qJD(6) * mrSges(7,3);
t286 = Ifges(4,5) * t356;
t280 = -pkin(4) * t296 - pkin(5);
t278 = pkin(4) * t292 + qJ(6);
t265 = Ifges(3,5) * t334;
t262 = Ifges(4,1) * t294 + t374;
t260 = Ifges(4,2) * t298 + t375;
t250 = -mrSges(5,1) * t298 - mrSges(5,3) * t363;
t249 = mrSges(5,2) * t298 - mrSges(5,3) * t364;
t248 = (Ifges(4,1) * t298 - t375) * qJD(3);
t247 = t315 * qJD(4);
t246 = (-Ifges(4,2) * t294 + t374) * qJD(3);
t245 = t314 * qJD(4);
t243 = (mrSges(4,1) * t294 + mrSges(4,2) * t298) * qJD(3);
t242 = t316 * qJD(4);
t235 = t316 * t294;
t221 = -Ifges(5,6) * t298 + t294 * t314;
t220 = -Ifges(5,3) * t298 + (Ifges(5,5) * t297 - t371) * t294;
t207 = -mrSges(5,2) * t357 - mrSges(5,3) * t305;
t206 = mrSges(5,1) * t357 - mrSges(5,3) * t306;
t200 = mrSges(6,2) * t298 - mrSges(6,3) * t226;
t199 = -mrSges(7,2) * t226 - mrSges(7,3) * t298;
t196 = -mrSges(4,1) * t365 - t229 * mrSges(4,3);
t195 = mrSges(4,2) * t365 - t228 * mrSges(4,3);
t183 = mrSges(6,1) * t239 + mrSges(6,2) * t240;
t182 = mrSges(7,1) * t239 - mrSges(7,3) * t240;
t167 = mrSges(5,1) * t305 + mrSges(5,2) * t306;
t166 = pkin(5) * t239 - qJ(6) * t240 + t281;
t163 = mrSges(6,1) * t226 - mrSges(6,2) * t227;
t162 = mrSges(7,1) * t226 + mrSges(7,3) * t227;
t161 = -t261 * t354 + (Ifges(5,5) * t294 + t298 * t315) * qJD(3);
t160 = -t259 * t354 + (Ifges(5,6) * t294 + t298 * t314) * qJD(3);
t158 = mrSges(4,1) * t335 - mrSges(4,3) * t191;
t157 = -mrSges(4,2) * t335 - mrSges(4,3) * t190;
t156 = Ifges(4,1) * t229 - Ifges(4,4) * t228 - Ifges(4,5) * t365;
t155 = Ifges(4,4) * t229 - Ifges(4,2) * t228 - t345;
t151 = -Ifges(7,4) * t227 - Ifges(7,2) * t298 + Ifges(7,6) * t226;
t150 = -Ifges(6,5) * t227 - Ifges(6,6) * t226 - Ifges(6,3) * t298;
t146 = t360 - t405;
t144 = mrSges(5,1) * t228 + mrSges(5,3) * t307;
t143 = -mrSges(5,2) * t228 + mrSges(5,3) * t192;
t142 = pkin(5) * t226 + qJ(6) * t227 + t253;
t128 = -mrSges(5,1) * t192 - mrSges(5,2) * t307;
t122 = mrSges(4,1) * t190 + mrSges(4,2) * t191;
t114 = mrSges(6,1) * t181 - mrSges(6,2) * t180;
t113 = mrSges(7,1) * t181 + mrSges(7,3) * t180;
t107 = -mrSges(6,2) * t357 - mrSges(6,3) * t130;
t105 = mrSges(6,1) * t357 - mrSges(6,3) * t129;
t104 = -mrSges(7,2) * t130 + mrSges(7,3) * t357;
t102 = Ifges(4,1) * t191 - Ifges(4,4) * t190 + Ifges(4,5) * t335;
t101 = Ifges(4,4) * t191 - Ifges(4,2) * t190 + Ifges(4,6) * t335;
t100 = pkin(5) * t298 - t108;
t99 = -qJ(6) * t298 + t400;
t90 = -Ifges(5,4) * t307 + Ifges(5,2) * t192 + Ifges(5,6) * t228;
t89 = -Ifges(5,5) * t307 + Ifges(5,6) * t192 + Ifges(5,3) * t228;
t85 = pkin(5) * t181 + qJ(6) * t180 - qJD(6) * t240 + t348;
t81 = -mrSges(6,2) * t228 + mrSges(6,3) * t308;
t80 = mrSges(7,2) * t308 + mrSges(7,3) * t228;
t74 = mrSges(5,1) * t190 - mrSges(5,3) * t97;
t73 = -mrSges(5,2) * t190 + mrSges(5,3) * t96;
t70 = mrSges(6,1) * t130 + mrSges(6,2) * t129;
t69 = mrSges(7,1) * t130 - mrSges(7,3) * t129;
t62 = -mrSges(6,1) * t308 + mrSges(6,2) * t112;
t61 = -mrSges(7,1) * t308 - mrSges(7,3) * t112;
t55 = Ifges(7,4) * t112 + Ifges(7,2) * t228 - Ifges(7,6) * t308;
t54 = Ifges(6,5) * t112 + Ifges(6,6) * t308 + Ifges(6,3) * t228;
t52 = -mrSges(5,1) * t96 + mrSges(5,2) * t97;
t49 = pkin(5) * t130 - qJ(6) * t129 + qJD(6) * t227 + t211;
t46 = Ifges(5,1) * t97 + Ifges(5,4) * t96 + Ifges(5,5) * t190;
t45 = Ifges(5,4) * t97 + Ifges(5,2) * t96 + Ifges(5,6) * t190;
t39 = -pkin(5) * t308 - qJ(6) * t112 + t86;
t29 = -mrSges(7,2) * t36 + mrSges(7,3) * t190;
t28 = -mrSges(6,2) * t190 - mrSges(6,3) * t36;
t26 = mrSges(6,1) * t190 - mrSges(6,3) * t35;
t21 = -pkin(5) * t228 - t22;
t20 = qJ(6) * t228 + t377;
t15 = mrSges(6,1) * t36 + mrSges(6,2) * t35;
t14 = mrSges(7,1) * t36 - mrSges(7,3) * t35;
t7 = pkin(5) * t36 - qJ(6) * t35 - qJD(6) * t112 + t47;
t1 = [(-t299 * t337 + 0.2e1 * (t224 * t299 + t225 * t295) * mrSges(3,3) + ((t233 * t393 + Ifges(3,5) * t291 + 0.2e1 * (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t299) * t290) * t299 + (t234 * t393 + Ifges(4,5) * t229 - 0.2e1 * Ifges(3,6) * t291 - Ifges(4,6) * t228 + (-0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t295 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3)) * t299) * t290) * t295) * qJD(2)) * t290 + 0.2e1 * t377 * t28 + (t22 * t6 + t377 * t5 + t47 * t86) * t396 + (t89 + t54 + t55 - t155) * t190 + (t12 + t13) * t112 - (-t11 + t8) * t308 - t307 * t46 + 0.2e1 * t217 * t122 + 0.2e1 * t78 * t195 + 0.2e1 * t79 * t196 + t191 * t156 + t192 * t45 + 0.2e1 * t147 * t158 + 0.2e1 * t148 * t157 + 0.2e1 * t140 * t52 + 0.2e1 * t24 * t143 + 0.2e1 * t25 * t144 + 0.2e1 * t77 * t128 + t96 * t90 + t97 * t91 + 0.2e1 * t5 * t81 + 0.2e1 * t6 * t82 + 0.2e1 * t3 * t83 + 0.2e1 * t86 * t15 + 0.2e1 * t2 * t80 + 0.2e1 * t72 * t73 + 0.2e1 * t71 * t74 + 0.2e1 * t7 * t61 + 0.2e1 * t47 * t62 + 0.2e1 * t39 * t14 + 0.2e1 * t22 * t26 + 0.2e1 * t21 * t27 + 0.2e1 * t20 * t29 + 0.2e1 * m(4) * (t147 * t79 + t148 * t78 + t217 * t225) + 0.2e1 * m(3) * (t224 * t234 - t225 * t233) + (t102 + 0.2e1 * t367) * t229 + (t265 - 0.2e1 * t369 - 0.2e1 * t370) * t291 + (-t101 + t44 + 0.2e1 * t368 + t402) * t228 + (t2 * t20 + t21 * t3 + t39 * t7) * t395 + (t140 * t77 + t24 * t72 + t25 * t71) * t397 + (t57 + t58) * t35 + (t53 - t56) * t36; t377 * t107 + m(7) * (t100 * t3 + t142 * t7 + t2 * t99 + t20 * t38 + t21 * t40 + t39 * t49) + (-t246 / 0.2e1 + t159 / 0.2e1 + t64 / 0.2e1 + t65 / 0.2e1) * t228 + (-t260 / 0.2e1 + t220 / 0.2e1 + t150 / 0.2e1 + t151 / 0.2e1) * t190 + m(6) * (t108 * t6 + t211 * t86 + t22 * t43 + t253 * t47 + t377 * t42 + t400 * t5) + t400 * t28 - t341 * t308 + t229 * t248 / 0.2e1 + t24 * t249 + t25 * t250 + t253 * t15 + t191 * t262 / 0.2e1 + t217 * t243 + t77 * t235 + t3 * t202 + t71 * t206 + t72 * t207 + t211 * t62 + t215 * t74 + t216 * t73 + t2 * t199 + t5 * t200 + t6 * t201 + t7 * t162 + t47 * t163 + t140 * t167 + t145 * t143 + t146 * t144 + t142 * t14 - pkin(2) * t122 + t21 * t106 + t108 * t26 + t99 * t29 + t100 * t27 + t20 * t104 + t22 * t105 + t42 * t81 + t43 * t82 + t40 * t83 + t86 * t70 + t38 * t80 + t39 * t69 + t49 * t61 + ((t156 / 0.2e1 - t147 * mrSges(4,3) + t90 * t385 + t338) * t298 + (t345 / 0.2e1 + t55 / 0.2e1 - t155 / 0.2e1 + t89 / 0.2e1 + t54 / 0.2e1 - t148 * mrSges(4,3)) * t294 + (-t294 * t195 + (t128 - t196) * t298 + m(4) * (-t147 * t298 - t148 * t294) + t140 * t382) * pkin(9)) * qJD(3) + (t46 * t383 + t45 * t385 - t79 * mrSges(4,3) + t367 + t102 / 0.2e1 + Ifges(4,5) * t319 + (t384 * t90 + t385 * t91) * qJD(4) + (-t158 + t52) * pkin(9)) * t294 + m(4) * (-pkin(2) * t225 - t79 * t288 + t78 * t378) + m(5) * (t145 * t72 + t146 * t71 + t215 * t25 + t216 * t24 + t288 * t77) + (pkin(9) * t157 + t78 * mrSges(4,3) - t368 + t101 / 0.2e1 - t44 / 0.2e1 - t9 / 0.2e1 - t10 / 0.2e1 + Ifges(4,6) * t319) * t298 + (-t299 * t286 / 0.2e1 - Ifges(3,6) * t358) * t290 - t344 * t227 + t349 * t226 + t340 * t112 + t342 * t129 + t343 * t130 + t265 + t97 * t387 + t161 * t388 + t160 * t389 + t221 * t391 - t369 - t370 + t325 * t35 + t326 * t36; -(t67 + t68) * t227 + (t63 - t66) * t226 + (t149 - t152) * t130 + (t153 + t154) * t129 + 0.2e1 * t400 * t107 + (t108 * t43 + t211 * t253 + t400 * t42) * t396 + 0.2e1 * t145 * t249 + 0.2e1 * t146 * t250 + 0.2e1 * t253 * t70 + (t145 * t216 + t146 * t215) * t397 - 0.2e1 * pkin(2) * t243 + 0.2e1 * t211 * t163 + 0.2e1 * t215 * t206 + 0.2e1 * t216 * t207 + 0.2e1 * t38 * t199 + 0.2e1 * t42 * t200 + 0.2e1 * t43 * t201 + 0.2e1 * t40 * t202 + 0.2e1 * t49 * t162 + 0.2e1 * t142 * t69 + 0.2e1 * t108 * t105 + 0.2e1 * t99 * t104 + 0.2e1 * t100 * t106 + (t167 * t394 - t293 * t160 + t297 * t161 + t248 + (-t221 * t297 - t222 * t293) * qJD(4) + (0.2e1 * pkin(9) ^ 2 * t382 + t150 + t151 + t220 - t260) * qJD(3)) * t294 + (-t159 + t246 + (-t293 * t221 + t297 * t222 + t235 * t394 + t262) * qJD(3) + t401) * t298 + (t100 * t40 + t142 * t49 + t38 * t99) * t395; m(6) * (t134 * t377 - t135 * t22 - t197 * t6 + t198 * t5 + t281 * t47 + t348 * t86) + t337 + m(7) * (t134 * t20 + t135 * t21 + t166 * t7 + t197 * t3 + t198 * t2 + t39 * t85) + t46 * t403 + (t28 + t29) * t198 + (t27 - t26) * t197 + (t81 + t80) * t134 + (-t180 * t21 - t181 * t20 - t2 * t239 + t240 * t3) * mrSges(7,2) + (t180 * t22 - t181 * t377 - t239 * t5 - t240 * t6) * mrSges(6,3) - t328 * t308 + t281 * t15 + t140 * t242 + t7 * t182 + t47 * t183 + t166 * t14 + t39 * t113 + t86 * t114 + t85 * t61 - t78 * mrSges(4,2) + t79 * mrSges(4,1) - pkin(3) * t52 + t376 * t135 + (m(5) * (-t353 * t71 - t355 * t72 + t312) + t297 * t73 - t293 * t74 - t143 * t355 - t144 * t353) * pkin(10) + t344 * t240 + t349 * t239 - t342 * t180 + t343 * t181 + t399 * t77 + (t338 + (-t90 / 0.2e1 + pkin(4) * t62) * t293) * qJD(4) + t45 * t383 + t97 * t386 + t247 * t388 + t245 * t389 + t259 * t391 + ((-t293 * t72 - t297 * t71) * qJD(4) + t312) * mrSges(5,3) + t317 * t190 + t318 * t228 + t323 * t35 + t324 * t36 + t327 * t112; (-t146 * mrSges(5,3) + t161 / 0.2e1 - t259 * t356 / 0.2e1 + (-t221 / 0.2e1 - t216 * mrSges(5,3) + (m(6) * t253 + t163) * pkin(4)) * qJD(4) + (m(5) * (-t146 - t405) - t206 - qJD(4) * t249) * pkin(10)) * t293 + (-t100 * t180 - t181 * t99 - t239 * t38 + t240 * t40) * mrSges(7,2) + (t108 * t180 - t181 * t400 - t239 * t42 - t240 * t43) * mrSges(6,3) + m(6) * (-t108 * t135 + t134 * t400 - t197 * t43 + t198 * t42 + t211 * t281) + t281 * t70 + t253 * t114 + t211 * t183 + t49 * t182 + t85 * t162 + t166 * t69 - pkin(3) * t167 + t142 * t113 + (t107 + t104) * t198 + (t106 - t105) * t197 + (t199 + t200) * t134 + m(7) * (t100 * t135 + t134 * t99 + t142 * t85 + t166 * t49 + t197 * t40 + t198 * t38) + (qJD(4) * t387 + t160 / 0.2e1 + t356 * t386 + t321 * mrSges(5,3) + (m(5) * t321 - qJD(4) * t250 + t207) * pkin(10)) * t297 + (t247 * t383 + t245 * t385 + pkin(9) * t242 + (t259 * t384 + t261 * t385) * qJD(4) + (pkin(9) * mrSges(4,2) - Ifges(4,6) + t317) * qJD(3)) * t294 + t361 * t135 + t340 * t240 + t341 * t239 + ((-mrSges(4,1) + t399) * qJD(3) * pkin(9) - t318) * t298 + t286 + t323 * t129 + t324 * t130 - t325 * t180 + t326 * t181 - t327 * t227 + t328 * t226; -0.2e1 * pkin(3) * t242 + 0.2e1 * t166 * t113 + 0.2e1 * t281 * t114 + 0.2e1 * t85 * t182 + t297 * t245 + t293 * t247 + (t119 + t120) * t240 + (t115 - t118) * t239 + (t184 - t187) * t181 - (t188 + t189) * t180 + (t297 * t261 + (0.2e1 * pkin(4) * t183 - t259) * t293) * qJD(4) + (t281 * t348 + t322) * t396 + (t166 * t85 + t322) * t395 + 0.2e1 * (mrSges(6,3) + mrSges(7,2)) * (-t134 * t239 + t135 * t240 - t180 * t197 - t181 * t198); m(7) * (t2 * t278 + t20 * t270 + t280 * t3) + (m(6) * (t292 * t5 + t296 * t6) + t296 * t26 + t292 * t28 + (t296 * t81 + t376 * t292 + m(6) * (-t22 * t292 + t296 * t377) + t21 * t381) * qJD(5)) * pkin(4) + t278 * t29 + t280 * t27 + t270 * t80 + t303 - t24 * mrSges(5,2) + t25 * mrSges(5,1) + t44; t302 + t159 + m(7) * (t270 * t99 + t278 * t38 + t280 * t40) + (m(6) * (t292 * t42 + t296 * t43) + t296 * t105 + t292 * t107 + (t296 * t200 + t361 * t292 + m(6) * (-t108 * t292 + t296 * t400) + t100 * t381) * qJD(5)) * pkin(4) + t278 * t104 + t280 * t106 + t270 * t199 - t145 * mrSges(5,2) + t146 * mrSges(5,1); m(7) * (t134 * t278 + t135 * t280 + t198 * t270) + t285 + (pkin(10) * t257 - t371) * qJD(4) + (-t180 * t280 - t181 * t278 - t239 * t270) * mrSges(7,2) + (t240 * mrSges(7,2) * t352 + m(7) * t330 + m(6) * (t134 * t292 - t135 * t296 + t198 * t351 + t330) + (t296 * t180 - t292 * t181 + (-t239 * t296 + t240 * t292) * qJD(5)) * mrSges(6,3)) * pkin(4) + t304; 0.2e1 * m(7) * (t270 * t278 + t280 * t347) + 0.2e1 * t300; m(7) * (-pkin(5) * t3 + qJ(6) * t2 + qJD(6) * t20) + t303 + qJD(6) * t80 - pkin(5) * t27 + qJ(6) * t29; t302 + m(7) * (-pkin(5) * t40 + qJ(6) * t38 + qJD(6) * t99) + qJD(6) * t199 + qJ(6) * t104 - pkin(5) * t106; m(7) * (-pkin(5) * t135 + qJ(6) * t134 + qJD(6) * t198) + (pkin(5) * t180 - qJ(6) * t181 - qJD(6) * t239) * mrSges(7,2) + t304; t289 + m(7) * (-pkin(5) * t347 + qJ(6) * t270 + qJD(6) * t278) + t300; 0.2e1 * m(7) * qJ(6) * qJD(6) + 0.2e1 * t289; m(7) * t3 + t27; m(7) * t40 + t106; m(7) * t135 - t180 * mrSges(7,2); m(7) * t347; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
