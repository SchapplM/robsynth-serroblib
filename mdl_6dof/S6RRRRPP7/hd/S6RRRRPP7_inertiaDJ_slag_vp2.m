% Calculate time derivative of joint inertia matrix for
% S6RRRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-03-09 21:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP7_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP7_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP7_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP7_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:17:12
% EndTime: 2019-03-09 21:17:30
% DurationCPUTime: 8.63s
% Computational Cost: add. (10296->773), mult. (28025->1083), div. (0->0), fcn. (26180->10), ass. (0->300)
t290 = sin(qJ(4));
t291 = sin(qJ(3));
t293 = cos(qJ(4));
t340 = qJD(4) * t293;
t294 = cos(qJ(3));
t343 = qJD(3) * t294;
t298 = t290 * t343 + t291 * t340;
t341 = qJD(4) * t291;
t299 = -t290 * t341 + t293 * t343;
t344 = qJD(3) * t291;
t150 = Ifges(5,5) * t299 - Ifges(5,6) * t298 + Ifges(5,3) * t344;
t287 = sin(pkin(11));
t355 = cos(pkin(11));
t312 = t355 * t293;
t313 = t355 * t290;
t142 = -t299 * t287 - t312 * t341 - t313 * t343;
t237 = t287 * t293 + t313;
t354 = t287 * t290;
t301 = t312 - t354;
t143 = -t237 * t341 + t301 * t343;
t68 = Ifges(6,5) * t143 + Ifges(6,6) * t142 + Ifges(6,3) * t344;
t69 = Ifges(7,4) * t143 + Ifges(7,2) * t344 - Ifges(7,6) * t142;
t395 = -t68 - t69 - t150;
t255 = -pkin(3) * t294 - pkin(10) * t291 - pkin(2);
t349 = t293 * t294;
t274 = pkin(9) * t349;
t206 = t290 * t255 + t274;
t394 = qJD(4) * t206;
t392 = t290 / 0.2e1;
t372 = t293 / 0.2e1;
t288 = sin(pkin(6));
t292 = sin(qJ(2));
t353 = t288 * t292;
t269 = pkin(8) * t353;
t289 = cos(pkin(6));
t295 = cos(qJ(2));
t370 = pkin(1) * t295;
t208 = t269 + (-pkin(2) - t370) * t289;
t228 = -t289 * t294 + t291 * t353;
t229 = t289 * t291 + t294 * t353;
t113 = t228 * pkin(3) - t229 * pkin(10) + t208;
t352 = t288 * t295;
t233 = t289 * t292 * pkin(1) + pkin(8) * t352;
t209 = pkin(9) * t289 + t233;
t210 = (-pkin(2) * t295 - pkin(9) * t292 - pkin(1)) * t288;
t126 = t294 * t209 + t291 * t210;
t115 = -pkin(10) * t352 + t126;
t61 = t290 * t113 + t293 * t115;
t389 = -m(5) * pkin(3) - mrSges(5,1) * t293 + mrSges(5,2) * t290;
t232 = t289 * t370 - t269;
t346 = qJD(2) * t288;
t325 = t295 * t346;
t180 = qJD(3) * t229 + t291 * t325;
t181 = -qJD(3) * t228 + t294 * t325;
t302 = -t229 * t293 + t290 * t352;
t345 = qJD(2) * t292;
t326 = t288 * t345;
t94 = qJD(4) * t302 - t181 * t290 + t293 * t326;
t182 = -t229 * t290 - t293 * t352;
t95 = qJD(4) * t182 + t181 * t293 + t290 * t326;
t45 = t287 * t95 - t355 * t94;
t47 = t287 * t94 + t355 * t95;
t11 = Ifges(6,5) * t47 - Ifges(6,6) * t45 + Ifges(6,3) * t180;
t12 = Ifges(7,4) * t47 + Ifges(7,2) * t180 + Ifges(7,6) * t45;
t29 = Ifges(5,5) * t95 + Ifges(5,6) * t94 + Ifges(5,3) * t180;
t388 = t11 + t12 + t29;
t369 = pkin(4) * t287;
t275 = qJ(6) + t369;
t387 = m(7) * t275 + mrSges(7,3);
t386 = 0.2e1 * m(5);
t385 = 2 * m(6);
t384 = 0.2e1 * m(7);
t383 = 0.2e1 * pkin(9);
t382 = -2 * mrSges(3,3);
t380 = m(6) * pkin(4);
t379 = t94 / 0.2e1;
t216 = (pkin(2) * t292 - pkin(9) * t295) * t346;
t217 = t232 * qJD(2);
t73 = -t209 * t344 + t210 * t343 + t291 * t216 + t294 * t217;
t65 = pkin(10) * t326 + t73;
t218 = t233 * qJD(2);
t88 = t180 * pkin(3) - t181 * pkin(10) + t218;
t23 = -t61 * qJD(4) - t290 * t65 + t293 * t88;
t7 = pkin(4) * t180 - qJ(5) * t95 + qJD(5) * t302 + t23;
t342 = qJD(4) * t290;
t22 = t113 * t340 - t115 * t342 + t290 * t88 + t293 * t65;
t9 = qJ(5) * t94 + qJD(5) * t182 + t22;
t4 = t287 * t7 + t355 * t9;
t378 = t182 / 0.2e1;
t377 = -t302 / 0.2e1;
t363 = Ifges(5,4) * t290;
t305 = Ifges(5,1) * t293 - t363;
t213 = -Ifges(5,5) * t294 + t291 * t305;
t376 = t213 / 0.2e1;
t362 = Ifges(5,4) * t293;
t261 = Ifges(5,1) * t290 + t362;
t375 = t261 / 0.2e1;
t374 = -t290 / 0.2e1;
t373 = -t293 / 0.2e1;
t371 = m(5) * t294;
t368 = pkin(9) * t290;
t367 = pkin(9) * t294;
t286 = t291 * pkin(9);
t366 = -qJ(5) - pkin(10);
t60 = t293 * t113 - t115 * t290;
t35 = pkin(4) * t228 + qJ(5) * t302 + t60;
t56 = qJ(5) * t182 + t61;
t19 = t287 * t35 + t355 * t56;
t339 = qJD(5) * t293;
t250 = (pkin(3) * t291 - pkin(10) * t294) * qJD(3);
t347 = t293 * t250 + t344 * t368;
t90 = -t291 * t339 + (pkin(4) * t291 - qJ(5) * t349) * qJD(3) + (-t274 + (qJ(5) * t291 - t255) * t290) * qJD(4) + t347;
t348 = t290 * t250 + t255 * t340;
t350 = t291 * t293;
t361 = pkin(9) * qJD(3);
t98 = (-qJ(5) * qJD(4) - t361) * t350 + (-qJD(5) * t291 + (-pkin(9) * qJD(4) - qJ(5) * qJD(3)) * t294) * t290 + t348;
t46 = t287 * t90 + t355 * t98;
t365 = Ifges(4,4) * t291;
t364 = Ifges(4,4) * t294;
t360 = t217 * mrSges(3,2);
t359 = t218 * mrSges(3,1);
t358 = t218 * mrSges(4,1);
t357 = t218 * mrSges(4,2);
t227 = t301 * qJD(4);
t356 = t227 * mrSges(7,2);
t351 = t290 * t291;
t239 = t293 * t255;
t165 = -qJ(5) * t350 + t239 + (-pkin(4) - t368) * t294;
t184 = -qJ(5) * t351 + t206;
t103 = t287 * t165 + t355 * t184;
t226 = t237 * qJD(4);
t156 = Ifges(6,5) * t227 - Ifges(6,6) * t226;
t157 = Ifges(7,4) * t227 + Ifges(7,6) * t226;
t251 = pkin(4) * t351 + t286;
t337 = pkin(4) * t342;
t336 = Ifges(4,6) * t352;
t10 = Ifges(7,5) * t47 + Ifges(7,6) * t180 + Ifges(7,3) * t45;
t13 = Ifges(6,4) * t47 - Ifges(6,2) * t45 + Ifges(6,6) * t180;
t335 = t10 / 0.2e1 - t13 / 0.2e1;
t14 = Ifges(7,1) * t47 + Ifges(7,4) * t180 + Ifges(7,5) * t45;
t15 = Ifges(6,1) * t47 - Ifges(6,4) * t45 + Ifges(6,5) * t180;
t334 = t15 / 0.2e1 + t14 / 0.2e1;
t104 = -t182 * t355 - t287 * t302;
t105 = t287 * t182 - t302 * t355;
t48 = Ifges(7,5) * t105 + Ifges(7,6) * t228 + Ifges(7,3) * t104;
t51 = Ifges(6,4) * t105 - Ifges(6,2) * t104 + Ifges(6,6) * t228;
t333 = t51 / 0.2e1 - t48 / 0.2e1;
t52 = Ifges(7,1) * t105 + Ifges(7,4) * t228 + Ifges(7,5) * t104;
t53 = Ifges(6,1) * t105 - Ifges(6,4) * t104 + Ifges(6,5) * t228;
t332 = t52 / 0.2e1 + t53 / 0.2e1;
t67 = Ifges(7,5) * t143 + Ifges(7,6) * t344 - Ifges(7,3) * t142;
t70 = Ifges(6,4) * t143 + Ifges(6,2) * t142 + Ifges(6,6) * t344;
t331 = t67 / 0.2e1 - t70 / 0.2e1;
t71 = Ifges(7,1) * t143 + Ifges(7,4) * t344 - Ifges(7,5) * t142;
t72 = Ifges(6,1) * t143 + Ifges(6,4) * t142 + Ifges(6,5) * t344;
t330 = t71 / 0.2e1 + t72 / 0.2e1;
t87 = -Ifges(5,1) * t302 + Ifges(5,4) * t182 + Ifges(5,5) * t228;
t329 = t87 * t372;
t328 = Ifges(4,5) * t181 - Ifges(4,6) * t180 + Ifges(4,3) * t326;
t201 = pkin(4) * t298 + pkin(9) * t343;
t279 = -pkin(4) * t293 - pkin(3);
t327 = t355 * pkin(4);
t214 = t237 * t291;
t215 = t301 * t291;
t127 = Ifges(7,5) * t215 - Ifges(7,6) * t294 + Ifges(7,3) * t214;
t130 = Ifges(6,4) * t215 - Ifges(6,2) * t214 - Ifges(6,6) * t294;
t320 = t127 / 0.2e1 - t130 / 0.2e1;
t131 = Ifges(7,1) * t215 - Ifges(7,4) * t294 + Ifges(7,5) * t214;
t132 = Ifges(6,1) * t215 - Ifges(6,4) * t214 - Ifges(6,5) * t294;
t319 = t131 / 0.2e1 + t132 / 0.2e1;
t155 = Ifges(7,5) * t227 + Ifges(7,3) * t226;
t158 = Ifges(6,4) * t227 - Ifges(6,2) * t226;
t318 = -t158 / 0.2e1 + t155 / 0.2e1;
t159 = Ifges(7,1) * t227 + Ifges(7,5) * t226;
t160 = Ifges(6,1) * t227 - Ifges(6,4) * t226;
t317 = t159 / 0.2e1 + t160 / 0.2e1;
t169 = Ifges(7,5) * t237 - Ifges(7,3) * t301;
t172 = Ifges(6,4) * t237 + Ifges(6,2) * t301;
t316 = t169 / 0.2e1 - t172 / 0.2e1;
t173 = Ifges(7,1) * t237 - Ifges(7,5) * t301;
t174 = Ifges(6,1) * t237 + Ifges(6,4) * t301;
t315 = t173 / 0.2e1 + t174 / 0.2e1;
t21 = t45 * mrSges(6,1) + t47 * mrSges(6,2);
t20 = t45 * mrSges(7,1) - t47 * mrSges(7,3);
t28 = -t180 * mrSges(7,1) + t47 * mrSges(7,2);
t314 = qJD(4) * t366;
t76 = -t142 * mrSges(6,1) + t143 * mrSges(6,2);
t154 = t226 * mrSges(6,1) + t227 * mrSges(6,2);
t75 = -t142 * mrSges(7,1) - t143 * mrSges(7,3);
t153 = t226 * mrSges(7,1) - t227 * mrSges(7,3);
t225 = t290 * t314 + t339;
t297 = -qJD(5) * t290 + t293 * t314;
t144 = t225 * t287 - t297 * t355;
t145 = t225 * t355 + t287 * t297;
t257 = t366 * t293;
t185 = -t257 * t287 - t313 * t366;
t186 = -t257 * t355 + t354 * t366;
t311 = t144 * t185 + t186 * t145;
t125 = -t291 * t209 + t210 * t294;
t123 = (-t293 * t344 - t294 * t342) * pkin(9) + t348;
t205 = -t290 * t367 + t239;
t310 = -qJD(4) * t205 + t123;
t309 = t326 / 0.2e1;
t120 = -mrSges(7,1) * t344 + t143 * mrSges(7,2);
t243 = Ifges(5,5) * t340 - Ifges(5,6) * t342;
t114 = pkin(3) * t352 - t125;
t308 = t243 / 0.2e1 + t156 / 0.2e1 + t157 / 0.2e1;
t307 = Ifges(5,5) * t392 + Ifges(5,6) * t372 + (Ifges(6,5) + Ifges(7,4)) * t237 / 0.2e1 + (Ifges(6,6) / 0.2e1 - Ifges(7,6) / 0.2e1) * t301;
t306 = mrSges(5,1) * t290 + mrSges(5,2) * t293;
t304 = -Ifges(5,2) * t290 + t362;
t259 = Ifges(5,2) * t293 + t363;
t303 = t22 * t293 - t23 * t290;
t74 = -t209 * t343 - t210 * t344 + t216 * t294 - t291 * t217;
t3 = -t287 * t9 + t355 * t7;
t18 = -t287 * t56 + t35 * t355;
t44 = -t287 * t98 + t355 * t90;
t102 = t165 * t355 - t287 * t184;
t81 = -pkin(4) * t182 + t114;
t66 = -pkin(3) * t326 - t74;
t32 = -pkin(4) * t94 + t66;
t284 = Ifges(4,5) * t343;
t278 = -t327 - pkin(5);
t264 = Ifges(3,5) * t325;
t262 = Ifges(4,1) * t291 + t364;
t260 = Ifges(4,2) * t294 + t365;
t249 = -mrSges(5,1) * t294 - mrSges(5,3) * t350;
t248 = mrSges(5,2) * t294 - mrSges(5,3) * t351;
t247 = (Ifges(4,1) * t294 - t365) * qJD(3);
t246 = t305 * qJD(4);
t245 = (-Ifges(4,2) * t291 + t364) * qJD(3);
t244 = t304 * qJD(4);
t242 = (mrSges(4,1) * t291 + mrSges(4,2) * t294) * qJD(3);
t241 = t306 * qJD(4);
t234 = t306 * t291;
t212 = -Ifges(5,6) * t294 + t291 * t304;
t211 = -Ifges(5,3) * t294 + (Ifges(5,5) * t293 - Ifges(5,6) * t290) * t291;
t197 = -mrSges(5,2) * t344 - mrSges(5,3) * t298;
t196 = mrSges(5,1) * t344 - mrSges(5,3) * t299;
t192 = mrSges(7,1) * t294 + mrSges(7,2) * t215;
t191 = -mrSges(6,1) * t294 - mrSges(6,3) * t215;
t190 = mrSges(6,2) * t294 - mrSges(6,3) * t214;
t189 = -mrSges(7,2) * t214 - mrSges(7,3) * t294;
t188 = -mrSges(4,1) * t352 - t229 * mrSges(4,3);
t187 = mrSges(4,2) * t352 - t228 * mrSges(4,3);
t168 = -mrSges(6,1) * t301 + mrSges(6,2) * t237;
t167 = -mrSges(7,1) * t301 - mrSges(7,3) * t237;
t163 = mrSges(5,1) * t298 + mrSges(5,2) * t299;
t161 = -pkin(5) * t301 - qJ(6) * t237 + t279;
t152 = -t261 * t341 + (Ifges(5,5) * t291 + t294 * t305) * qJD(3);
t151 = -t259 * t341 + (Ifges(5,6) * t291 + t294 * t304) * qJD(3);
t149 = mrSges(4,1) * t326 - mrSges(4,3) * t181;
t148 = -mrSges(4,2) * t326 - mrSges(4,3) * t180;
t147 = mrSges(6,1) * t214 + mrSges(6,2) * t215;
t146 = mrSges(7,1) * t214 - mrSges(7,3) * t215;
t141 = Ifges(4,1) * t229 - Ifges(4,4) * t228 - Ifges(4,5) * t352;
t140 = Ifges(4,4) * t229 - Ifges(4,2) * t228 - t336;
t129 = Ifges(7,4) * t215 - Ifges(7,2) * t294 + Ifges(7,6) * t214;
t128 = Ifges(6,5) * t215 - Ifges(6,6) * t214 - Ifges(6,3) * t294;
t124 = t347 - t394;
t122 = mrSges(5,1) * t228 + mrSges(5,3) * t302;
t121 = -mrSges(5,2) * t228 + mrSges(5,3) * t182;
t119 = mrSges(6,1) * t344 - mrSges(6,3) * t143;
t118 = -mrSges(6,2) * t344 + mrSges(6,3) * t142;
t117 = mrSges(7,2) * t142 + mrSges(7,3) * t344;
t116 = pkin(5) * t214 - qJ(6) * t215 + t251;
t111 = pkin(5) * t226 - qJ(6) * t227 - qJD(6) * t237 + t337;
t107 = -mrSges(5,1) * t182 - mrSges(5,2) * t302;
t106 = mrSges(4,1) * t180 + mrSges(4,2) * t181;
t100 = Ifges(4,1) * t181 - Ifges(4,4) * t180 + Ifges(4,5) * t326;
t99 = Ifges(4,4) * t181 - Ifges(4,2) * t180 + Ifges(4,6) * t326;
t97 = t294 * pkin(5) - t102;
t96 = -qJ(6) * t294 + t103;
t86 = -Ifges(5,4) * t302 + Ifges(5,2) * t182 + Ifges(5,6) * t228;
t85 = -Ifges(5,5) * t302 + Ifges(5,6) * t182 + Ifges(5,3) * t228;
t80 = -mrSges(7,1) * t228 + mrSges(7,2) * t105;
t79 = mrSges(6,1) * t228 - mrSges(6,3) * t105;
t78 = -mrSges(6,2) * t228 - mrSges(6,3) * t104;
t77 = -mrSges(7,2) * t104 + mrSges(7,3) * t228;
t63 = mrSges(5,1) * t180 - mrSges(5,3) * t95;
t62 = -mrSges(5,2) * t180 + mrSges(5,3) * t94;
t59 = -pkin(5) * t142 - qJ(6) * t143 - qJD(6) * t215 + t201;
t58 = mrSges(6,1) * t104 + mrSges(6,2) * t105;
t57 = mrSges(7,1) * t104 - mrSges(7,3) * t105;
t54 = -mrSges(5,1) * t94 + mrSges(5,2) * t95;
t50 = Ifges(7,4) * t105 + Ifges(7,2) * t228 + Ifges(7,6) * t104;
t49 = Ifges(6,5) * t105 - Ifges(6,6) * t104 + Ifges(6,3) * t228;
t36 = -pkin(5) * t344 - t44;
t34 = qJ(6) * t344 - qJD(6) * t294 + t46;
t31 = Ifges(5,1) * t95 + Ifges(5,4) * t94 + Ifges(5,5) * t180;
t30 = Ifges(5,4) * t95 + Ifges(5,2) * t94 + Ifges(5,6) * t180;
t27 = mrSges(6,1) * t180 - mrSges(6,3) * t47;
t26 = -mrSges(6,2) * t180 - mrSges(6,3) * t45;
t25 = -mrSges(7,2) * t45 + mrSges(7,3) * t180;
t24 = pkin(5) * t104 - qJ(6) * t105 + t81;
t17 = -t228 * pkin(5) - t18;
t16 = qJ(6) * t228 + t19;
t5 = pkin(5) * t45 - qJ(6) * t47 - qJD(6) * t105 + t32;
t2 = -t180 * pkin(5) - t3;
t1 = qJ(6) * t180 + qJD(6) * t228 + t4;
t6 = [-t302 * t31 + (t48 - t51) * t45 + (-t99 + 0.2e1 * t358 + t388) * t228 + (-t295 * t328 + 0.2e1 * (t217 * t295 + t218 * t292) * mrSges(3,3) + ((t232 * t382 + Ifges(3,5) * t289 + 0.2e1 * (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t295) * t288) * t295 + (t233 * t382 + Ifges(4,5) * t229 - 0.2e1 * Ifges(3,6) * t289 - Ifges(4,6) * t228 + (-0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t292 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3)) * t295) * t288) * t292) * qJD(2)) * t288 + (t1 * t16 + t17 * t2 + t24 * t5) * t384 + (t18 * t3 + t19 * t4 + t32 * t81) * t385 + (t114 * t66 + t22 * t61 + t23 * t60) * t386 + (t85 + t49 + t50 - t140) * t180 + (t15 + t14) * t105 + (t10 - t13) * t104 + 0.2e1 * t208 * t106 + 0.2e1 * t73 * t187 + 0.2e1 * t74 * t188 + t181 * t141 + t182 * t30 + 0.2e1 * t126 * t148 + 0.2e1 * t125 * t149 + 0.2e1 * t22 * t121 + 0.2e1 * t23 * t122 + 0.2e1 * t66 * t107 + 0.2e1 * t114 * t54 + t95 * t87 + t94 * t86 + 0.2e1 * t1 * t77 + 0.2e1 * t4 * t78 + 0.2e1 * t3 * t79 + 0.2e1 * t2 * t80 + 0.2e1 * t81 * t21 + 0.2e1 * t61 * t62 + 0.2e1 * t60 * t63 + 0.2e1 * t5 * t57 + 0.2e1 * t32 * t58 + 0.2e1 * t16 * t25 + 0.2e1 * t19 * t26 + 0.2e1 * t18 * t27 + 0.2e1 * t17 * t28 + 0.2e1 * t24 * t20 + 0.2e1 * m(3) * (t217 * t233 - t218 * t232) + 0.2e1 * m(4) * (t125 * t74 + t126 * t73 + t208 * t218) + (t52 + t53) * t47 + (t100 + 0.2e1 * t357) * t229 + (t264 - 0.2e1 * t359 - 0.2e1 * t360) * t289; -t359 - t360 + m(6) * (t102 * t3 + t103 * t4 + t18 * t44 + t19 * t46 + t201 * t81 + t251 * t32) + m(7) * (t1 * t96 + t116 * t5 + t16 * t34 + t17 * t36 + t2 * t97 + t24 * t59) + t95 * t376 + t152 * t377 + t151 * t378 + t212 * t379 + t181 * t262 / 0.2e1 + t208 * t242 + t229 * t247 / 0.2e1 + t22 * t248 + t23 * t249 + t251 * t21 + t66 * t234 + t205 * t63 + t206 * t62 + t1 * t189 + t4 * t190 + t3 * t191 + t2 * t192 + t60 * t196 + t61 * t197 + t201 * t58 + t114 * t163 + t5 * t146 + t32 * t147 + t116 * t20 + t16 * t117 + t19 * t118 + t18 * t119 + t17 * t120 + t123 * t121 + t124 * t122 - pkin(2) * t106 + t96 * t25 + t97 * t28 + t102 * t27 + t103 * t26 + t34 * t77 + t46 * t78 + t44 * t79 + t36 * t80 + t81 * t76 + t24 * t75 + t59 * t57 + t264 + (-t245 / 0.2e1 + t150 / 0.2e1 + t68 / 0.2e1 + t69 / 0.2e1) * t228 + (-t260 / 0.2e1 + t211 / 0.2e1 + t128 / 0.2e1 + t129 / 0.2e1) * t180 + t319 * t47 + t320 * t45 + t330 * t105 + t331 * t104 + t332 * t143 + t333 * t142 + t334 * t215 + t335 * t214 + (-t295 * t284 / 0.2e1 - Ifges(3,6) * t345) * t288 + (Ifges(4,6) * t309 - t358 + t99 / 0.2e1 - t29 / 0.2e1 - t11 / 0.2e1 - t12 / 0.2e1 + pkin(9) * t148 + t73 * mrSges(4,3)) * t294 + m(5) * (t123 * t61 + t124 * t60 + t205 * t23 + t206 * t22 + t286 * t66) + m(4) * (-pkin(2) * t218 - t286 * t74 + t367 * t73) + ((t329 + t86 * t374 - t125 * mrSges(4,3) + t141 / 0.2e1) * t294 + (-t126 * mrSges(4,3) - t140 / 0.2e1 + t85 / 0.2e1 + t49 / 0.2e1 + t50 / 0.2e1 + t336 / 0.2e1) * t291 + (-t291 * t187 + (t107 - t188) * t294 + m(4) * (-t125 * t294 - t126 * t291) + t114 * t371) * pkin(9)) * qJD(3) + (Ifges(4,5) * t309 + t357 + t100 / 0.2e1 - t74 * mrSges(4,3) + t30 * t374 + t31 * t372 + (t373 * t86 + t374 * t87) * qJD(4) + (-t149 + t54) * pkin(9)) * t291; (t116 * t59 + t34 * t96 + t36 * t97) * t384 + (t102 * t44 + t103 * t46 + t201 * t251) * t385 + (t123 * t206 + t124 * t205) * t386 - 0.2e1 * pkin(2) * t242 + 0.2e1 * t123 * t248 + 0.2e1 * t124 * t249 + 0.2e1 * t251 * t76 + 0.2e1 * t205 * t196 + 0.2e1 * t206 * t197 + 0.2e1 * t34 * t189 + 0.2e1 * t46 * t190 + 0.2e1 * t44 * t191 + 0.2e1 * t36 * t192 + 0.2e1 * t201 * t147 + 0.2e1 * t59 * t146 + 0.2e1 * t116 * t75 + 0.2e1 * t96 * t117 + 0.2e1 * t103 * t118 + 0.2e1 * t102 * t119 + 0.2e1 * t97 * t120 + (t71 + t72) * t215 + (-t70 + t67) * t214 + (t131 + t132) * t143 + (t130 - t127) * t142 + (t163 * t383 - t290 * t151 + t293 * t152 + t247 + (-t212 * t293 - t213 * t290) * qJD(4) + (0.2e1 * pkin(9) ^ 2 * t371 + t128 + t129 + t211 - t260) * qJD(3)) * t291 + (t245 + (-t290 * t212 + t213 * t293 + t234 * t383 + t262) * qJD(3) + t395) * t294; -t335 * t301 + (-t18 * t227 - t19 * t226 - t237 * t3 + t301 * t4) * mrSges(6,3) + (t1 * t301 - t16 * t226 + t17 * t227 + t2 * t237) * mrSges(7,2) + t30 * t372 + t95 * t375 + t246 * t377 + t244 * t378 + t259 * t379 + t279 * t21 + t114 * t241 + t161 * t20 + t5 * t167 + t32 * t168 + t24 * t153 + t81 * t154 + t111 * t57 - t73 * mrSges(4,2) + t74 * mrSges(4,1) + m(7) * (t1 * t186 + t111 * t24 + t144 * t17 + t145 * t16 + t161 * t5 + t185 * t2) - pkin(3) * t54 + t31 * t392 + t389 * t66 + t328 + (t25 + t26) * t186 + (-t27 + t28) * t185 + (t77 + t78) * t145 + (t80 - t79) * t144 + ((-t290 * t61 - t293 * t60) * qJD(4) + t303) * mrSges(5,3) + t307 * t180 + t308 * t228 + t315 * t47 + t316 * t45 + t317 * t105 + t318 * t104 + (t329 + (pkin(4) * t58 - t86 / 0.2e1) * t290) * qJD(4) + t332 * t227 - t333 * t226 + t334 * t237 + m(6) * (-t144 * t18 + t145 * t19 - t185 * t3 + t186 * t4 + t279 * t32 + t81 * t337) + (-t122 * t340 - t121 * t342 - t290 * t63 + t293 * t62 + m(5) * (-t340 * t60 - t342 * t61 + t303)) * pkin(10); (-t259 * t343 / 0.2e1 + t152 / 0.2e1 - t124 * mrSges(5,3) + (-t206 * mrSges(5,3) - t212 / 0.2e1 + (m(6) * t251 + t147) * pkin(4)) * qJD(4) + (-qJD(4) * t248 - t196 + m(5) * (-t124 - t394)) * pkin(10)) * t290 + (-t102 * t227 - t103 * t226 - t237 * t44 + t301 * t46) * mrSges(6,3) - t331 * t301 + (-t226 * t96 + t227 * t97 + t237 * t36 + t301 * t34) * mrSges(7,2) + (-t191 + t192) * t144 + t279 * t76 + t251 * t154 + t201 * t168 + t161 * t75 - pkin(3) * t163 + t59 * t167 + t111 * t146 + t116 * t153 + m(7) * (t111 * t116 + t144 * t97 + t145 * t96 + t161 * t59 + t185 * t36 + t186 * t34) + ((-mrSges(4,1) + t389) * t361 - t308) * t294 + t284 + (t117 + t118) * t186 + (t120 - t119) * t185 + (t189 + t190) * t145 + m(6) * (-t102 * t144 + t103 * t145 - t185 * t44 + t186 * t46 + t201 * t279) + t315 * t143 - t316 * t142 + t317 * t215 + t318 * t214 + t319 * t227 + t320 * t226 + t330 * t237 + (pkin(9) * t241 + t244 * t374 + t246 * t372 + (t259 * t373 + t261 * t374) * qJD(4) + (pkin(9) * mrSges(4,2) - Ifges(4,6) + t307) * qJD(3)) * t291 + (t343 * t375 + t151 / 0.2e1 + qJD(4) * t376 + t310 * mrSges(5,3) + (m(5) * t310 - qJD(4) * t249 + t197) * pkin(10)) * t293; -0.2e1 * pkin(3) * t241 + 0.2e1 * t111 * t167 + 0.2e1 * t161 * t153 + 0.2e1 * t279 * t154 + t293 * t244 + t290 * t246 + (t159 + t160) * t237 - (t155 - t158) * t301 + (t173 + t174) * t227 + (t169 - t172) * t226 + (t293 * t261 + (0.2e1 * pkin(4) * t168 - t259) * t290) * qJD(4) + (t111 * t161 + t311) * t384 + (t279 * t337 + t311) * t385 + 0.2e1 * (mrSges(6,3) + mrSges(7,2)) * (t144 * t237 + t145 * t301 + t185 * t227 - t186 * t226); t26 * t369 + (t287 * t4 + t3 * t355) * t380 + t275 * t25 + t278 * t28 + qJD(6) * t77 + t23 * mrSges(5,1) - t22 * mrSges(5,2) + t27 * t327 - t4 * mrSges(6,2) + t3 * mrSges(6,1) + t1 * mrSges(7,3) - t2 * mrSges(7,1) + m(7) * (qJD(6) * t16 + t1 * t275 + t2 * t278) + t388; t118 * t369 + (t287 * t46 + t355 * t44) * t380 + t275 * t117 + t278 * t120 + qJD(6) * t189 - t123 * mrSges(5,2) + t124 * mrSges(5,1) - t46 * mrSges(6,2) + t44 * mrSges(6,1) - t36 * mrSges(7,1) + t34 * mrSges(7,3) + t119 * t327 + m(7) * (qJD(6) * t96 + t275 * t34 + t278 * t36) - t395; m(7) * qJD(6) * t186 + t278 * t356 + t156 + t157 + t243 + (t287 * t380 - mrSges(6,2) + t387) * t145 + (m(7) * t278 - t355 * t380 - mrSges(6,1) - mrSges(7,1)) * t144 + (-mrSges(5,1) * t340 + mrSges(5,2) * t342) * pkin(10) + (-t226 * t369 - t227 * t327) * mrSges(6,3) + (qJD(6) * t301 - t226 * t275) * mrSges(7,2); 0.2e1 * t387 * qJD(6); m(6) * t32 + m(7) * t5 + t20 + t21; m(6) * t201 + m(7) * t59 + t75 + t76; m(6) * t337 + m(7) * t111 + t153 + t154; 0; 0; m(7) * t2 + t28; m(7) * t36 + t120; m(7) * t144 + t356; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
