% Calculate time derivative of joint inertia matrix for
% S6RRRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR12_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR12_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR12_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR12_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR12_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:38:16
% EndTime: 2019-03-09 19:38:33
% DurationCPUTime: 7.12s
% Computational Cost: add. (17562->773), mult. (45860->1141), div. (0->0), fcn. (45991->12), ass. (0->316)
t308 = sin(pkin(12));
t310 = cos(pkin(12));
t313 = sin(qJ(5));
t317 = cos(qJ(5));
t270 = t308 * t317 + t310 * t313;
t312 = sin(qJ(6));
t316 = cos(qJ(6));
t323 = t308 * t313 - t310 * t317;
t193 = -t270 * t312 - t316 * t323;
t378 = t193 / 0.2e1;
t194 = t270 * t316 - t312 * t323;
t377 = t194 / 0.2e1;
t369 = -t323 / 0.2e1;
t368 = t270 / 0.2e1;
t366 = t310 / 0.2e1;
t314 = sin(qJ(3));
t318 = cos(qJ(3));
t277 = -pkin(3) * t318 - qJ(4) * t314 - pkin(2);
t265 = t310 * t277;
t348 = t310 * t314;
t195 = -pkin(10) * t348 + t265 + (-pkin(9) * t308 - pkin(4)) * t318;
t347 = t310 * t318;
t229 = pkin(9) * t347 + t308 * t277;
t352 = t308 * t314;
t212 = -pkin(10) * t352 + t229;
t137 = t313 * t195 + t317 * t212;
t362 = pkin(10) + qJ(4);
t278 = t362 * t308;
t280 = t362 * t310;
t214 = -t313 * t278 + t317 * t280;
t309 = sin(pkin(6));
t315 = sin(qJ(2));
t350 = t309 * t315;
t294 = pkin(8) * t350;
t311 = cos(pkin(6));
t319 = cos(qJ(2));
t364 = pkin(1) * t319;
t259 = t311 * t364 - t294;
t252 = t323 * qJD(5);
t397 = 2 * m(5);
t396 = 2 * m(6);
t395 = 2 * m(7);
t394 = 0.2e1 * pkin(9);
t393 = -2 * mrSges(3,3);
t255 = t311 * t314 + t318 * t350;
t349 = t309 * t319;
t208 = -t255 * t308 - t310 * t349;
t209 = t255 * t310 - t308 * t349;
t138 = t208 * t317 - t209 * t313;
t139 = t208 * t313 + t209 * t317;
t254 = -t311 * t318 + t314 * t350;
t75 = Ifges(6,4) * t139 + Ifges(6,2) * t138 + Ifges(6,6) * t254;
t392 = t75 / 0.2e1;
t83 = t138 * t316 - t139 * t312;
t391 = t83 / 0.2e1;
t84 = t138 * t312 + t139 * t316;
t390 = t84 / 0.2e1;
t253 = t270 * qJD(5);
t128 = qJD(6) * t193 - t252 * t316 - t253 * t312;
t389 = t128 / 0.2e1;
t129 = -qJD(6) * t194 + t252 * t312 - t253 * t316;
t388 = t129 / 0.2e1;
t132 = Ifges(7,4) * t194 + Ifges(7,2) * t193;
t387 = t132 / 0.2e1;
t133 = Ifges(7,1) * t194 + Ifges(7,4) * t193;
t386 = t133 / 0.2e1;
t385 = t138 / 0.2e1;
t384 = t139 / 0.2e1;
t242 = t270 * t314;
t243 = t323 * t314;
t165 = -Ifges(6,4) * t243 - Ifges(6,2) * t242 - Ifges(6,6) * t318;
t383 = t165 / 0.2e1;
t167 = -t242 * t316 + t243 * t312;
t382 = t167 / 0.2e1;
t168 = -t242 * t312 - t243 * t316;
t381 = t168 / 0.2e1;
t345 = qJD(2) * t309;
t329 = t319 * t345;
t211 = -qJD(3) * t254 + t318 * t329;
t344 = qJD(2) * t315;
t330 = t309 * t344;
t170 = -t211 * t308 + t310 * t330;
t380 = t170 / 0.2e1;
t171 = t211 * t310 + t308 * t330;
t379 = t171 / 0.2e1;
t200 = Ifges(6,4) * t270 - Ifges(6,2) * t323;
t376 = t200 / 0.2e1;
t201 = Ifges(6,1) * t270 - Ifges(6,4) * t323;
t375 = t201 / 0.2e1;
t357 = Ifges(5,4) * t310;
t325 = -Ifges(5,2) * t308 + t357;
t223 = (Ifges(5,6) * t314 + t318 * t325) * qJD(3);
t374 = t223 / 0.2e1;
t358 = Ifges(5,4) * t308;
t326 = Ifges(5,1) * t310 - t358;
t224 = (Ifges(5,5) * t314 + t318 * t326) * qJD(3);
t373 = t224 / 0.2e1;
t372 = -t242 / 0.2e1;
t371 = -t243 / 0.2e1;
t370 = -t252 / 0.2e1;
t367 = -t308 / 0.2e1;
t365 = m(5) * t318;
t363 = pkin(5) * t253;
t305 = t314 * pkin(9);
t239 = t294 + (-pkin(2) - t364) * t311;
t148 = t254 * pkin(3) - t255 * qJ(4) + t239;
t260 = t311 * t315 * pkin(1) + pkin(8) * t349;
t240 = pkin(9) * t311 + t260;
t241 = (-pkin(2) * t319 - pkin(9) * t315 - pkin(1)) * t309;
t160 = t318 * t240 + t314 * t241;
t149 = -qJ(4) * t349 + t160;
t94 = t310 * t148 - t149 * t308;
t63 = pkin(4) * t254 - pkin(10) * t209 + t94;
t95 = t308 * t148 + t310 * t149;
t79 = pkin(10) * t208 + t95;
t31 = t313 * t63 + t317 * t79;
t210 = qJD(3) * t255 + t314 * t329;
t247 = t260 * qJD(2);
t101 = t210 * pkin(3) - t211 * qJ(4) - t255 * qJD(4) + t247;
t245 = (pkin(2) * t315 - pkin(9) * t319) * t345;
t246 = t259 * qJD(2);
t342 = qJD(3) * t318;
t343 = qJD(3) * t314;
t113 = -t240 * t343 + t241 * t342 + t314 * t245 + t318 * t246;
t97 = (qJ(4) * t344 - qJD(4) * t319) * t309 + t113;
t50 = t308 * t101 + t310 * t97;
t361 = mrSges(5,2) * t310;
t360 = Ifges(4,4) * t314;
t359 = Ifges(4,4) * t318;
t356 = t246 * mrSges(3,2);
t355 = t247 * mrSges(3,1);
t354 = t247 * mrSges(4,1);
t353 = t247 * mrSges(4,2);
t351 = t308 * t318;
t71 = Ifges(7,5) * t128 + Ifges(7,6) * t129;
t190 = -Ifges(6,5) * t252 - Ifges(6,6) * t253;
t328 = t308 * t342;
t244 = mrSges(5,1) * t328 + t342 * t361;
t251 = -qJD(4) * t314 + (pkin(3) * t314 - qJ(4) * t318) * qJD(3);
t334 = pkin(9) * t343;
t206 = t310 * t251 + t308 * t334;
t304 = pkin(9) * t342;
t261 = pkin(4) * t328 + t304;
t276 = pkin(4) * t352 + t305;
t341 = qJD(4) * t308;
t340 = qJD(4) * t310;
t339 = qJD(5) * t313;
t338 = qJD(5) * t317;
t337 = qJD(6) * t312;
t336 = qJD(6) * t316;
t64 = qJD(5) * t138 + t170 * t313 + t171 * t317;
t65 = -qJD(5) * t139 + t170 * t317 - t171 * t313;
t22 = qJD(6) * t83 + t312 * t65 + t316 * t64;
t23 = -qJD(6) * t84 - t312 * t64 + t316 * t65;
t6 = Ifges(7,5) * t22 + Ifges(7,6) * t23 + Ifges(7,3) * t210;
t26 = Ifges(6,5) * t64 + Ifges(6,6) * t65 + Ifges(6,3) * t210;
t179 = -t253 * t314 - t323 * t342;
t180 = t252 * t314 - t270 * t342;
t88 = qJD(6) * t167 + t179 * t316 + t180 * t312;
t89 = -qJD(6) * t168 - t179 * t312 + t180 * t316;
t39 = Ifges(7,5) * t88 + Ifges(7,6) * t89 + Ifges(7,3) * t343;
t333 = Ifges(4,6) * t349;
t332 = t190 / 0.2e1 + t71 / 0.2e1;
t106 = Ifges(6,5) * t179 + Ifges(6,6) * t180 + Ifges(6,3) * t343;
t331 = Ifges(4,5) * t211 - Ifges(4,6) * t210 + Ifges(4,3) * t330;
t300 = -pkin(4) * t310 - pkin(3);
t29 = -t65 * mrSges(6,1) + t64 * mrSges(6,2);
t9 = -t23 * mrSges(7,1) + t22 * mrSges(7,2);
t43 = -t89 * mrSges(7,1) + t88 * mrSges(7,2);
t49 = t310 * t101 - t308 * t97;
t30 = -t313 * t79 + t317 * t63;
t112 = -t170 * mrSges(5,1) + t171 * mrSges(5,2);
t116 = -t180 * mrSges(6,1) + t179 * mrSges(6,2);
t70 = -mrSges(7,1) * t129 + t128 * mrSges(7,2);
t136 = t317 * t195 - t212 * t313;
t159 = -t314 * t240 + t241 * t318;
t213 = -t317 * t278 - t280 * t313;
t150 = pkin(3) * t349 - t159;
t327 = Ifges(5,5) * t308 / 0.2e1 + Ifges(5,6) * t366 + Ifges(6,5) * t368 + Ifges(6,6) * t369 + Ifges(7,5) * t377 + Ifges(7,6) * t378;
t324 = Ifges(5,5) * t310 - Ifges(5,6) * t308;
t24 = pkin(5) * t254 - pkin(11) * t139 + t30;
t25 = pkin(11) * t138 + t31;
t12 = t24 * t316 - t25 * t312;
t13 = t24 * t312 + t25 * t316;
t105 = -pkin(5) * t318 + pkin(11) * t243 + t136;
t115 = -pkin(11) * t242 + t137;
t55 = t105 * t316 - t115 * t312;
t56 = t105 * t312 + t115 * t316;
t177 = -pkin(11) * t270 + t213;
t178 = -pkin(11) * t323 + t214;
t110 = t177 * t316 - t178 * t312;
t111 = t177 * t312 + t178 * t316;
t161 = -t278 * t338 + t317 * t340 + (-qJD(5) * t280 - t341) * t313;
t140 = -pkin(11) * t253 + t161;
t162 = -t270 * qJD(4) - qJD(5) * t214;
t141 = pkin(11) * t252 + t162;
t47 = qJD(6) * t110 + t140 * t316 + t141 * t312;
t48 = -qJD(6) * t111 - t140 * t312 + t141 * t316;
t322 = t48 * mrSges(7,1) - t47 * mrSges(7,2) + t71;
t114 = -t240 * t342 - t241 * t343 + t245 * t318 - t314 * t246;
t38 = pkin(4) * t210 - pkin(10) * t171 + t49;
t45 = pkin(10) * t170 + t50;
t11 = -qJD(5) * t31 - t313 * t45 + t317 * t38;
t4 = pkin(5) * t210 - pkin(11) * t64 + t11;
t10 = t313 * t38 + t317 * t45 + t63 * t338 - t339 * t79;
t5 = pkin(11) * t65 + t10;
t2 = qJD(6) * t12 + t312 * t4 + t316 * t5;
t3 = -qJD(6) * t13 - t312 * t5 + t316 * t4;
t321 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + t6;
t169 = (pkin(4) * t314 - pkin(10) * t347) * qJD(3) + t206;
t234 = t308 * t251;
t185 = t234 + (-pkin(9) * t348 - pkin(10) * t351) * qJD(3);
t67 = -qJD(5) * t137 + t317 * t169 - t185 * t313;
t51 = pkin(5) * t343 - pkin(11) * t179 + t67;
t66 = t313 * t169 + t317 * t185 + t195 * t338 - t212 * t339;
t52 = pkin(11) * t180 + t66;
t15 = qJD(6) * t55 + t312 * t51 + t316 * t52;
t16 = -qJD(6) * t56 - t312 * t52 + t316 * t51;
t320 = t16 * mrSges(7,1) - t15 * mrSges(7,2) + t39;
t119 = -pkin(4) * t208 + t150;
t100 = -pkin(3) * t330 - t114;
t82 = -pkin(4) * t170 + t100;
t303 = Ifges(4,5) * t342;
t289 = Ifges(3,5) * t329;
t285 = Ifges(4,1) * t314 + t359;
t284 = Ifges(4,2) * t318 + t360;
t283 = Ifges(5,1) * t308 + t357;
t282 = Ifges(5,2) * t310 + t358;
t279 = -mrSges(5,1) * t310 + mrSges(5,2) * t308;
t275 = (Ifges(4,1) * t318 - t360) * qJD(3);
t274 = (-Ifges(4,2) * t314 + t359) * qJD(3);
t273 = (mrSges(4,1) * t314 + mrSges(4,2) * t318) * qJD(3);
t272 = -mrSges(5,1) * t318 - mrSges(5,3) * t348;
t271 = mrSges(5,2) * t318 - mrSges(5,3) * t352;
t263 = (-mrSges(7,1) * t312 - mrSges(7,2) * t316) * qJD(6) * pkin(5);
t258 = (mrSges(5,1) * t314 - mrSges(5,3) * t347) * qJD(3);
t257 = (-mrSges(5,2) * t314 - mrSges(5,3) * t351) * qJD(3);
t256 = (mrSges(5,1) * t308 + t361) * t314;
t248 = t252 * mrSges(6,2);
t238 = -Ifges(5,5) * t318 + t314 * t326;
t237 = -Ifges(5,6) * t318 + t314 * t325;
t236 = -Ifges(5,3) * t318 + t314 * t324;
t233 = pkin(5) * t323 + t300;
t228 = -pkin(9) * t351 + t265;
t222 = (Ifges(5,3) * t314 + t318 * t324) * qJD(3);
t218 = -mrSges(6,1) * t318 + mrSges(6,3) * t243;
t217 = mrSges(6,2) * t318 - mrSges(6,3) * t242;
t216 = -mrSges(4,1) * t349 - t255 * mrSges(4,3);
t215 = mrSges(4,2) * t349 - t254 * mrSges(4,3);
t207 = -t310 * t334 + t234;
t198 = mrSges(6,1) * t323 + mrSges(6,2) * t270;
t197 = pkin(5) * t242 + t276;
t192 = -Ifges(6,1) * t252 - Ifges(6,4) * t253;
t191 = -Ifges(6,4) * t252 - Ifges(6,2) * t253;
t189 = mrSges(6,1) * t253 - t248;
t184 = mrSges(4,1) * t330 - mrSges(4,3) * t211;
t183 = -mrSges(4,2) * t330 - mrSges(4,3) * t210;
t182 = mrSges(6,1) * t242 - mrSges(6,2) * t243;
t176 = Ifges(4,1) * t255 - Ifges(4,4) * t254 - Ifges(4,5) * t349;
t175 = Ifges(4,4) * t255 - Ifges(4,2) * t254 - t333;
t166 = -Ifges(6,1) * t243 - Ifges(6,4) * t242 - Ifges(6,5) * t318;
t164 = -Ifges(6,5) * t243 - Ifges(6,6) * t242 - Ifges(6,3) * t318;
t156 = mrSges(5,1) * t254 - mrSges(5,3) * t209;
t155 = -mrSges(5,2) * t254 + mrSges(5,3) * t208;
t154 = -mrSges(6,2) * t343 + mrSges(6,3) * t180;
t153 = mrSges(6,1) * t343 - mrSges(6,3) * t179;
t152 = -mrSges(7,1) * t318 - mrSges(7,3) * t168;
t151 = mrSges(7,2) * t318 + mrSges(7,3) * t167;
t146 = -pkin(5) * t180 + t261;
t143 = mrSges(4,1) * t210 + mrSges(4,2) * t211;
t142 = -mrSges(5,1) * t208 + mrSges(5,2) * t209;
t135 = Ifges(4,1) * t211 - Ifges(4,4) * t210 + Ifges(4,5) * t330;
t134 = Ifges(4,4) * t211 - Ifges(4,2) * t210 + Ifges(4,6) * t330;
t130 = -mrSges(7,1) * t193 + mrSges(7,2) * t194;
t124 = mrSges(5,1) * t210 - mrSges(5,3) * t171;
t123 = -mrSges(5,2) * t210 + mrSges(5,3) * t170;
t122 = Ifges(5,1) * t209 + Ifges(5,4) * t208 + Ifges(5,5) * t254;
t121 = Ifges(5,4) * t209 + Ifges(5,2) * t208 + Ifges(5,6) * t254;
t120 = Ifges(5,5) * t209 + Ifges(5,6) * t208 + Ifges(5,3) * t254;
t118 = mrSges(6,1) * t254 - mrSges(6,3) * t139;
t117 = -mrSges(6,2) * t254 + mrSges(6,3) * t138;
t109 = -mrSges(7,1) * t167 + mrSges(7,2) * t168;
t108 = Ifges(6,1) * t179 + Ifges(6,4) * t180 + Ifges(6,5) * t343;
t107 = Ifges(6,4) * t179 + Ifges(6,2) * t180 + Ifges(6,6) * t343;
t104 = Ifges(7,1) * t168 + Ifges(7,4) * t167 - Ifges(7,5) * t318;
t103 = Ifges(7,4) * t168 + Ifges(7,2) * t167 - Ifges(7,6) * t318;
t102 = Ifges(7,5) * t168 + Ifges(7,6) * t167 - Ifges(7,3) * t318;
t93 = Ifges(5,1) * t171 + Ifges(5,4) * t170 + Ifges(5,5) * t210;
t92 = Ifges(5,4) * t171 + Ifges(5,2) * t170 + Ifges(5,6) * t210;
t91 = Ifges(5,5) * t171 + Ifges(5,6) * t170 + Ifges(5,3) * t210;
t90 = -mrSges(6,1) * t138 + mrSges(6,2) * t139;
t81 = -mrSges(7,2) * t343 + mrSges(7,3) * t89;
t80 = mrSges(7,1) * t343 - mrSges(7,3) * t88;
t77 = -pkin(5) * t138 + t119;
t76 = Ifges(6,1) * t139 + Ifges(6,4) * t138 + Ifges(6,5) * t254;
t74 = Ifges(6,5) * t139 + Ifges(6,6) * t138 + Ifges(6,3) * t254;
t73 = Ifges(7,1) * t128 + Ifges(7,4) * t129;
t72 = Ifges(7,4) * t128 + Ifges(7,2) * t129;
t69 = mrSges(7,1) * t254 - mrSges(7,3) * t84;
t68 = -mrSges(7,2) * t254 + mrSges(7,3) * t83;
t54 = -mrSges(6,2) * t210 + mrSges(6,3) * t65;
t53 = mrSges(6,1) * t210 - mrSges(6,3) * t64;
t42 = -mrSges(7,1) * t83 + mrSges(7,2) * t84;
t41 = Ifges(7,1) * t88 + Ifges(7,4) * t89 + Ifges(7,5) * t343;
t40 = Ifges(7,4) * t88 + Ifges(7,2) * t89 + Ifges(7,6) * t343;
t35 = Ifges(7,1) * t84 + Ifges(7,4) * t83 + Ifges(7,5) * t254;
t34 = Ifges(7,4) * t84 + Ifges(7,2) * t83 + Ifges(7,6) * t254;
t33 = Ifges(7,5) * t84 + Ifges(7,6) * t83 + Ifges(7,3) * t254;
t32 = -pkin(5) * t65 + t82;
t28 = Ifges(6,1) * t64 + Ifges(6,4) * t65 + Ifges(6,5) * t210;
t27 = Ifges(6,4) * t64 + Ifges(6,2) * t65 + Ifges(6,6) * t210;
t18 = -mrSges(7,2) * t210 + mrSges(7,3) * t23;
t17 = mrSges(7,1) * t210 - mrSges(7,3) * t22;
t8 = Ifges(7,1) * t22 + Ifges(7,4) * t23 + Ifges(7,5) * t210;
t7 = Ifges(7,4) * t22 + Ifges(7,2) * t23 + Ifges(7,6) * t210;
t1 = [0.2e1 * m(4) * (t113 * t160 + t114 * t159 + t239 * t247) + 0.2e1 * m(3) * (t246 * t260 - t247 * t259) + (t33 + t74 + t120 - t175) * t210 + (-t319 * t331 + 0.2e1 * (t246 * t319 + t247 * t315) * mrSges(3,3) + ((t259 * t393 + Ifges(3,5) * t311 + 0.2e1 * (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t319) * t309) * t319 + (t260 * t393 + Ifges(4,5) * t255 - 0.2e1 * Ifges(3,6) * t311 - Ifges(4,6) * t254 + (-0.2e1 * mrSges(3,1) * pkin(1) - 0.2e1 * Ifges(3,4) * t315 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3)) * t319) * t309) * t315) * qJD(2)) * t309 + (t100 * t150 + t49 * t94 + t50 * t95) * t397 + (t12 * t3 + t13 * t2 + t32 * t77) * t395 + (t10 * t31 + t11 * t30 + t119 * t82) * t396 + 0.2e1 * t239 * t143 + 0.2e1 * t113 * t215 + 0.2e1 * t114 * t216 + t208 * t92 + t209 * t93 + t211 * t176 + 0.2e1 * t160 * t183 + 0.2e1 * t159 * t184 + t170 * t121 + t171 * t122 + 0.2e1 * t150 * t112 + 0.2e1 * t50 * t155 + 0.2e1 * t49 * t156 + t139 * t28 + 0.2e1 * t100 * t142 + t138 * t27 + 0.2e1 * t12 * t17 + 0.2e1 * t13 * t18 + t23 * t34 + t22 * t35 + (t135 + 0.2e1 * t353) * t255 + (-t134 + t26 + t6 + t91 + 0.2e1 * t354) * t254 + 0.2e1 * t32 * t42 + (t289 - 0.2e1 * t355 - 0.2e1 * t356) * t311 + 0.2e1 * t30 * t53 + 0.2e1 * t31 * t54 + 0.2e1 * t2 * t68 + 0.2e1 * t3 * t69 + t65 * t75 + t64 * t76 + 0.2e1 * t77 * t9 + t83 * t7 + t84 * t8 + 0.2e1 * t82 * t90 + 0.2e1 * t10 * t117 + 0.2e1 * t11 * t118 + 0.2e1 * t119 * t29 + 0.2e1 * t95 * t123 + 0.2e1 * t94 * t124; (-t274 / 0.2e1 + t222 / 0.2e1 + t39 / 0.2e1 + t106 / 0.2e1) * t254 + (t318 * t183 + (t112 - t184) * t314) * pkin(9) + t238 * t379 + t237 * t380 + t8 * t381 + t7 * t382 + t65 * t383 + t108 * t384 + t107 * t385 + t41 * t390 + t40 * t391 + t180 * t392 + t28 * t371 + t27 * t372 + t209 * t373 + t208 * t374 + (-t284 / 0.2e1 + t236 / 0.2e1 + t164 / 0.2e1 + t102 / 0.2e1) * t210 - t355 - t356 + t289 + t211 * t285 / 0.2e1 + t50 * t271 + t49 * t272 + t239 * t273 + t255 * t275 / 0.2e1 + t276 * t29 + t100 * t256 + t95 * t257 + t94 * t258 + t261 * t90 + t150 * t244 + t10 * t217 + t11 * t218 + t228 * t124 + t229 * t123 + t206 * t156 + t207 * t155 + t197 * t9 + t179 * t76 / 0.2e1 + t82 * t182 + t64 * t166 / 0.2e1 + t2 * t151 + t3 * t152 + t30 * t153 + t31 * t154 - pkin(2) * t143 + t146 * t42 + t136 * t53 + t137 * t54 + (-t319 * t303 / 0.2e1 + (Ifges(4,5) * t314 / 0.2e1 + Ifges(4,6) * t318 / 0.2e1 - Ifges(3,6)) * t344) * t309 + (-t354 + t134 / 0.2e1 - t6 / 0.2e1 - t26 / 0.2e1 - t91 / 0.2e1 + t113 * mrSges(4,3)) * t318 + t55 * t17 + t56 * t18 + m(7) * (t12 * t16 + t13 * t15 + t146 * t77 + t197 * t32 + t2 * t56 + t3 * t55) + m(6) * (t10 * t137 + t11 * t136 + t119 * t261 + t276 * t82 + t30 * t67 + t31 * t66) + m(4) * (pkin(9) * t113 * t318 - pkin(2) * t247 - t114 * t305) + m(5) * (t100 * t305 + t206 * t94 + t207 * t95 + t228 * t49 + t229 * t50) + t15 * t68 + t16 * t69 + t77 * t43 + ((t122 * t366 + t121 * t367 - t159 * mrSges(4,3) + t176 / 0.2e1) * t318 + (-t160 * mrSges(4,3) + t33 / 0.2e1 + t74 / 0.2e1 + t120 / 0.2e1 - t175 / 0.2e1 + t333 / 0.2e1) * t314 + (-t314 * t215 + (t142 - t216) * t318 + t150 * t365 + m(4) * (-t159 * t318 - t160 * t314)) * pkin(9)) * qJD(3) + (t353 + t135 / 0.2e1 - t114 * mrSges(4,3) + t92 * t367 + t93 * t366) * t314 + t12 * t80 + t13 * t81 + t88 * t35 / 0.2e1 + t89 * t34 / 0.2e1 + t23 * t103 / 0.2e1 + t22 * t104 / 0.2e1 + t32 * t109 + t66 * t117 + t67 * t118 + t119 * t116; (-t106 + t274 - t222 - t39) * t318 + (t206 * t228 + t207 * t229) * t397 + (t146 * t197 + t15 * t56 + t16 * t55) * t395 + (t136 * t67 + t137 * t66 + t261 * t276) * t396 + 0.2e1 * t207 * t271 + 0.2e1 * t206 * t272 - 0.2e1 * pkin(2) * t273 + 0.2e1 * t276 * t116 + 0.2e1 * t229 * t257 + 0.2e1 * t228 * t258 + 0.2e1 * t261 * t182 - t243 * t108 - t242 * t107 + 0.2e1 * t66 * t217 + 0.2e1 * t67 * t218 + 0.2e1 * t197 * t43 + t179 * t166 + t180 * t165 + t167 * t40 + t168 * t41 + 0.2e1 * t15 * t151 + 0.2e1 * t16 * t152 + 0.2e1 * t136 * t153 + 0.2e1 * t137 * t154 + 0.2e1 * t146 * t109 + 0.2e1 * t55 * t80 + 0.2e1 * t56 * t81 + t89 * t103 + t88 * t104 + ((-t237 * t308 + t238 * t310 + t256 * t394 + t285) * t318 + (0.2e1 * pkin(9) ^ 2 * t365 + t102 + t164 + t236 - t284) * t314) * qJD(3) + (-t223 * t308 + t224 * t310 + t244 * t394 + t275) * t314; m(5) * (-pkin(3) * t100 + (-t308 * t94 + t310 * t95) * qJD(4) + (-t308 * t49 + t310 * t50) * qJ(4)) + (t92 / 0.2e1 + t50 * mrSges(5,3) + qJ(4) * t123 + qJD(4) * t155) * t310 + (t93 / 0.2e1 - t49 * mrSges(5,3) - qJ(4) * t124 - qJD(4) * t156) * t308 + (-t10 * t323 - t11 * t270 + t252 * t30 - t253 * t31) * mrSges(6,3) + m(7) * (t110 * t3 + t111 * t2 + t12 * t48 + t13 * t47 + t233 * t32 + t363 * t77) + t7 * t378 + t283 * t379 + t282 * t380 + t192 * t384 + t191 * t385 + t22 * t386 + t23 * t387 + t34 * t388 + t35 * t389 + t73 * t390 + t72 * t391 + t28 * t368 + t27 * t369 + t76 * t370 + t64 * t375 + t65 * t376 + t8 * t377 + m(6) * (t10 * t214 + t11 * t213 + t161 * t31 + t162 * t30 + t300 * t82) + t331 + t300 * t29 + t100 * t279 + t233 * t9 + t213 * t53 + t214 * t54 + t82 * t198 + t119 * t189 + t161 * t117 + t162 * t118 + t32 * t130 + (-t12 * t128 + t129 * t13 + t193 * t2 - t194 * t3) * mrSges(7,3) + t327 * t210 + t332 * t254 + t47 * t68 + t48 * t69 + t77 * t70 - (-pkin(5) * t42 + t392) * t253 + t110 * t17 + t111 * t18 - pkin(3) * t112 - t113 * mrSges(4,2) + t114 * mrSges(4,1); (-t128 * t55 + t129 * t56 + t15 * t193 - t16 * t194) * mrSges(7,3) + (t136 * t252 - t137 * t253 - t270 * t67 - t323 * t66) * mrSges(6,3) + m(7) * (t110 * t16 + t111 * t15 + t146 * t233 + t197 * t363 + t47 * t56 + t48 * t55) + m(5) * (-t228 * t341 + t229 * t340 + (-t206 * t308 + t207 * t310) * qJ(4)) + t73 * t381 + t72 * t382 + t88 * t386 + t89 * t387 + t103 * t388 + t104 * t389 + t108 * t368 + t107 * t369 + t166 * t370 + t192 * t371 + t191 * t372 + t179 * t375 + t180 * t376 + t41 * t377 + t40 * t378 + t300 * t116 + t276 * t189 + t261 * t198 - pkin(3) * t244 + t233 * t43 + t213 * t153 + t214 * t154 + t161 * t217 + t162 * t218 + t197 * t70 + t47 * t151 + t48 * t152 + t146 * t130 + m(6) * (t136 * t162 + t137 * t161 + t213 * t67 + t214 * t66 + t261 * t300) + t303 - t332 * t318 + ((t283 * t366 + t282 * t367 + (-m(5) * pkin(3) - mrSges(4,1) + t279) * pkin(9)) * t318 + (pkin(9) * mrSges(4,2) - Ifges(4,6) + t327) * t314) * qJD(3) + (-t206 * mrSges(5,3) - qJ(4) * t258 - qJD(4) * t272 + t373) * t308 + (t207 * mrSges(5,3) + qJ(4) * t257 + qJD(4) * t271 + t374) * t310 - (-pkin(5) * t109 + t383) * t253 + t110 * t80 + t111 * t81; t128 * t133 + t129 * t132 + 0.2e1 * t300 * t189 - t323 * t191 + t270 * t192 + t193 * t72 + t194 * t73 - t252 * t201 + 0.2e1 * t233 * t70 - (-0.2e1 * pkin(5) * t130 + t200) * t253 + (t110 * t48 + t111 * t47 + t233 * t363) * t395 + (t161 * t214 + t162 * t213) * t396 + 0.2e1 * (-t110 * t128 + t111 * t129 + t193 * t47 - t194 * t48) * mrSges(7,3) + 0.2e1 * (-t161 * t323 - t162 * t270 + t213 * t252 - t214 * t253) * mrSges(6,3) + (qJ(4) * t397 + 0.2e1 * mrSges(5,3)) * qJD(4) * (t308 ^ 2 + t310 ^ 2); m(5) * t100 + m(6) * t82 + m(7) * t32 + t112 + t29 + t9; m(5) * t304 + m(6) * t261 + m(7) * t146 + t116 + t244 + t43; -t248 - (-m(7) * pkin(5) - mrSges(6,1)) * t253 + t70; 0; t11 * mrSges(6,1) - t10 * mrSges(6,2) + (m(7) * (-t12 * t337 + t13 * t336 + t2 * t312 + t3 * t316) + t68 * t336 + t312 * t18 - t69 * t337 + t316 * t17) * pkin(5) + t321 + t26; t67 * mrSges(6,1) - t66 * mrSges(6,2) + (m(7) * (t15 * t312 + t16 * t316 + t336 * t56 - t337 * t55) + t151 * t336 + t312 * t81 - t152 * t337 + t316 * t80) * pkin(5) + t320 + t106; t162 * mrSges(6,1) - t161 * mrSges(6,2) + (m(7) * (t312 * t47 + t316 * t48 + (-t110 * t312 + t111 * t316) * qJD(6)) + (-t316 * t128 + t312 * t129 + (t193 * t316 + t194 * t312) * qJD(6)) * mrSges(7,3)) * pkin(5) + t322 + t190; 0; 0.2e1 * t263; t321; t320; t322; 0; t263; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
