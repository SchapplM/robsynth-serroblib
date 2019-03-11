% Calculate time derivative of joint inertia matrix for
% S6RRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR7_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR7_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR7_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR7_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:26:43
% EndTime: 2019-03-09 22:27:02
% DurationCPUTime: 9.00s
% Computational Cost: add. (19638->675), mult. (50277->992), div. (0->0), fcn. (50815->12), ass. (0->284)
t398 = Ifges(5,3) + Ifges(6,3);
t271 = sin(qJ(6));
t275 = cos(qJ(6));
t244 = -mrSges(7,1) * t275 + mrSges(7,2) * t271;
t396 = -mrSges(6,1) + t244;
t272 = sin(qJ(4));
t273 = sin(qJ(3));
t276 = cos(qJ(4));
t277 = cos(qJ(3));
t233 = -t272 * t273 + t276 * t277;
t389 = qJD(3) + qJD(4);
t191 = t389 * t233;
t234 = t272 * t277 + t273 * t276;
t192 = t389 * t234;
t268 = sin(pkin(12));
t346 = cos(pkin(12));
t153 = t191 * t346 - t268 * t192;
t186 = t268 * t233 + t234 * t346;
t321 = qJD(6) * t275;
t287 = t271 * t153 + t186 * t321;
t322 = qJD(6) * t271;
t332 = t275 * t153;
t286 = t186 * t322 - t332;
t377 = -pkin(10) - pkin(9);
t317 = t377 * t273;
t395 = t377 * t277;
t397 = -t272 * t317 + t276 * t395;
t200 = t272 * t395 + t276 * t317;
t367 = t271 / 0.2e1;
t366 = t275 / 0.2e1;
t302 = -t322 / 0.2e1;
t336 = t268 * t272;
t355 = pkin(3) * qJD(4);
t217 = (t276 * t346 - t336) * t355;
t299 = (t271 ^ 2 + t275 ^ 2) * t217;
t270 = cos(pkin(6));
t274 = sin(qJ(2));
t269 = sin(pkin(6));
t278 = cos(qJ(2));
t334 = t269 * t278;
t230 = t270 * t274 * pkin(1) + pkin(8) * t334;
t211 = pkin(9) * t270 + t230;
t212 = (-pkin(2) * t278 - pkin(9) * t274 - pkin(1)) * t269;
t171 = -t273 * t211 + t277 * t212;
t335 = t269 * t274;
t224 = t270 * t273 + t277 * t335;
t142 = -pkin(3) * t334 - t224 * pkin(10) + t171;
t172 = t277 * t211 + t273 * t212;
t223 = t270 * t277 - t273 * t335;
t158 = pkin(10) * t223 + t172;
t90 = t272 * t142 + t276 * t158;
t255 = pkin(8) * t335;
t364 = pkin(1) * t278;
t229 = t270 * t364 - t255;
t185 = -t233 * t346 + t234 * t268;
t344 = t186 * t271;
t139 = -mrSges(7,2) * t185 - mrSges(7,3) * t344;
t343 = t186 * t275;
t140 = mrSges(7,1) * t185 - mrSges(7,3) * t343;
t394 = -t139 * t322 - t140 * t321;
t262 = -pkin(3) * t277 - pkin(2);
t209 = -pkin(4) * t233 + t262;
t125 = pkin(5) * t185 - pkin(11) * t186 + t209;
t179 = qJ(5) * t233 - t397;
t281 = -t234 * qJ(5) + t200;
t127 = t179 * t346 + t268 * t281;
t79 = t125 * t275 - t127 * t271;
t80 = t125 * t271 + t127 * t275;
t393 = -t79 * t321 - t80 * t322;
t180 = t223 * t276 - t224 * t272;
t181 = t223 * t272 + t224 * t276;
t129 = t268 * t180 + t181 * t346;
t117 = -t271 * t129 - t275 * t334;
t128 = -t180 * t346 + t181 * t268;
t76 = -mrSges(7,2) * t128 + mrSges(7,3) * t117;
t288 = -t275 * t129 + t271 * t334;
t77 = mrSges(7,1) * t128 + mrSges(7,3) * t288;
t392 = -t77 * t321 - t76 * t322;
t89 = t276 * t142 - t272 * t158;
t75 = -pkin(4) * t334 - t181 * qJ(5) + t89;
t83 = qJ(5) * t180 + t90;
t38 = t268 * t75 + t346 * t83;
t36 = -pkin(11) * t334 + t38;
t210 = t255 + (-pkin(2) - t364) * t270;
t183 = -t223 * pkin(3) + t210;
t134 = -t180 * pkin(4) + t183;
t55 = t128 * pkin(5) - t129 * pkin(11) + t134;
t16 = -t271 * t36 + t275 * t55;
t17 = t271 * t55 + t275 * t36;
t391 = -t16 * t321 - t17 * t322;
t152 = t191 * t268 + t192 * t346;
t390 = Ifges(5,5) * t191 + Ifges(6,5) * t153 - Ifges(5,6) * t192 - Ifges(6,6) * t152;
t328 = qJD(2) * t269;
t307 = t278 * t328;
t196 = -qJD(3) * t224 - t273 * t307;
t197 = qJD(3) * t223 + t277 * t307;
t327 = qJD(2) * t274;
t308 = t269 * t327;
t388 = Ifges(4,5) * t197 + Ifges(4,6) * t196 + Ifges(4,3) * t308;
t387 = 2 * m(5);
t386 = 2 * m(6);
t385 = 2 * m(7);
t384 = -2 * mrSges(3,3);
t383 = -2 * mrSges(6,3);
t296 = qJD(4) * t395;
t297 = qJD(4) * t317;
t161 = t200 * qJD(3) + t272 * t296 + t276 * t297;
t114 = -t192 * qJ(5) + t233 * qJD(5) + t161;
t282 = t397 * qJD(3);
t162 = qJD(4) * t397 + t282;
t290 = -t191 * qJ(5) - t234 * qJD(5);
t70 = t114 * t268 - t346 * (t162 + t290);
t382 = 0.2e1 * t70;
t126 = t179 * t268 - t281 * t346;
t381 = 0.2e1 * t126;
t219 = t230 * qJD(2);
t380 = 0.2e1 * t219;
t379 = m(6) * pkin(4);
t112 = qJD(4) * t180 + t196 * t272 + t197 * t276;
t113 = -qJD(4) * t181 + t196 * t276 - t197 * t272;
t68 = t112 * t346 + t268 * t113;
t44 = qJD(6) * t288 - t271 * t68 + t275 * t308;
t378 = t44 / 0.2e1;
t376 = t117 / 0.2e1;
t375 = t233 / 0.2e1;
t374 = t234 / 0.2e1;
t263 = Ifges(7,5) * t321;
t373 = Ifges(7,6) * t302 + t263 / 0.2e1;
t358 = Ifges(7,4) * t271;
t294 = Ifges(7,1) * t275 - t358;
t242 = t294 * qJD(6);
t372 = t242 / 0.2e1;
t371 = Ifges(7,5) * t367 + Ifges(7,6) * t366;
t357 = Ifges(7,4) * t275;
t248 = Ifges(7,1) * t271 + t357;
t369 = t248 / 0.2e1;
t368 = -t271 / 0.2e1;
t363 = pkin(4) * t268;
t67 = t112 * t268 - t113 * t346;
t169 = -t196 * pkin(3) + t219;
t88 = -t113 * pkin(4) + t169;
t18 = t67 * pkin(5) - t68 * pkin(11) + t88;
t216 = (pkin(2) * t274 - pkin(9) * t278) * t328;
t218 = t229 * qJD(2);
t122 = -qJD(3) * t172 + t277 * t216 - t218 * t273;
t96 = pkin(3) * t308 - pkin(10) * t197 + t122;
t325 = qJD(3) * t277;
t326 = qJD(3) * t273;
t121 = -t211 * t326 + t212 * t325 + t273 * t216 + t277 * t218;
t98 = pkin(10) * t196 + t121;
t33 = -qJD(4) * t90 - t272 * t98 + t276 * t96;
t23 = pkin(4) * t308 - qJ(5) * t112 - qJD(5) * t181 + t33;
t323 = qJD(4) * t276;
t324 = qJD(4) * t272;
t32 = t142 * t323 - t158 * t324 + t272 * t96 + t276 * t98;
t27 = qJ(5) * t113 + qJD(5) * t180 + t32;
t9 = t268 * t23 + t346 * t27;
t6 = pkin(11) * t308 + t9;
t2 = qJD(6) * t16 + t18 * t271 + t275 * t6;
t362 = t2 * t275;
t3 = -qJD(6) * t17 + t18 * t275 - t271 * t6;
t361 = t3 * t271;
t360 = Ifges(4,4) * t273;
t359 = Ifges(4,4) * t277;
t356 = Ifges(7,6) * t271;
t354 = t126 * t70;
t71 = t346 * t114 + (-t272 * t297 + t276 * t296 + t282 + t290) * t268;
t265 = pkin(3) * t326;
t182 = pkin(4) * t192 + t265;
t86 = pkin(5) * t152 - pkin(11) * t153 + t182;
t20 = qJD(6) * t79 + t271 * t86 + t275 * t71;
t353 = t20 * t275;
t21 = -qJD(6) * t80 - t271 * t71 + t275 * t86;
t352 = t21 * t271;
t350 = t218 * mrSges(3,2);
t349 = t219 * mrSges(3,1);
t348 = t271 * mrSges(7,3);
t120 = -mrSges(6,1) * t334 - t129 * mrSges(6,3);
t73 = -mrSges(7,1) * t117 - mrSges(7,2) * t288;
t347 = t73 - t120;
t300 = t346 * t272;
t215 = (t268 * t276 + t300) * t355;
t345 = t126 * t215;
t261 = pkin(3) * t276 + pkin(4);
t220 = -pkin(3) * t336 + t261 * t346;
t213 = -pkin(5) - t220;
t295 = mrSges(7,1) * t271 + mrSges(7,2) * t275;
t237 = t295 * qJD(6);
t342 = t213 * t237;
t341 = t217 * t271;
t340 = t217 * t275;
t259 = pkin(11) + t363;
t339 = t259 * t271;
t338 = t259 * t275;
t309 = t346 * pkin(4);
t260 = -t309 - pkin(5);
t337 = t260 * t237;
t221 = pkin(3) * t300 + t268 * t261;
t320 = 0.2e1 * t269;
t43 = qJD(6) * t117 + t271 * t308 + t275 * t68;
t12 = Ifges(7,5) * t43 + Ifges(7,6) * t44 + Ifges(7,3) * t67;
t30 = t67 * mrSges(6,1) + t68 * mrSges(6,2);
t301 = t321 / 0.2e1;
t91 = t152 * mrSges(6,1) + t153 * mrSges(6,2);
t293 = -Ifges(7,2) * t271 + t357;
t292 = -t271 * t77 + t275 * t76;
t291 = t139 * t275 - t140 * t271;
t289 = Ifges(5,5) * t112 + Ifges(6,5) * t68 + Ifges(5,6) * t113 - Ifges(6,6) * t67 + t398 * t308;
t8 = t23 * t346 - t268 * t27;
t37 = -t268 * t83 + t346 * t75;
t240 = t293 * qJD(6);
t246 = Ifges(7,2) * t275 + t358;
t285 = t275 * t240 + t271 * t242 - t246 * t322 + t248 * t321;
t284 = -t361 + (-t16 * t275 - t17 * t271) * qJD(6);
t283 = -t352 + (-t271 * t80 - t275 * t79) * qJD(6);
t49 = -t286 * Ifges(7,5) - Ifges(7,6) * t287 + Ifges(7,3) * t152;
t13 = Ifges(7,4) * t43 + Ifges(7,2) * t44 + Ifges(7,6) * t67;
t14 = Ifges(7,1) * t43 + Ifges(7,4) * t44 + Ifges(7,5) * t67;
t35 = pkin(5) * t334 - t37;
t5 = -pkin(5) * t308 - t8;
t53 = -Ifges(7,4) * t288 + Ifges(7,2) * t117 + Ifges(7,6) * t128;
t54 = -Ifges(7,1) * t288 + Ifges(7,4) * t117 + Ifges(7,5) * t128;
t280 = t33 * mrSges(5,1) + t8 * mrSges(6,1) - t32 * mrSges(5,2) - t9 * mrSges(6,2) + mrSges(7,3) * t362 + t13 * t366 + t14 * t367 + t35 * t237 + t5 * t244 + t246 * t378 + t54 * t301 + t53 * t302 + t43 * t369 + t240 * t376 - t288 * t372 + t128 * t373 + t289 + t67 * t371;
t106 = Ifges(7,6) * t185 + t186 * t293;
t107 = Ifges(7,5) * t185 + t186 * t294;
t50 = -Ifges(7,4) * t286 - Ifges(7,2) * t287 + Ifges(7,6) * t152;
t51 = -Ifges(7,1) * t286 - Ifges(7,4) * t287 + Ifges(7,5) * t152;
t279 = -t161 * mrSges(5,2) - t71 * mrSges(6,2) + mrSges(7,3) * t353 + t107 * t301 + t126 * t237 + t50 * t366 + t51 * t367 + t332 * t369 + t152 * t371 + t162 * mrSges(5,1) - t240 * t344 / 0.2e1 + t343 * t372 + t185 * t373 + t390 + t396 * t70 - t287 * t246 / 0.2e1 + (t186 * t248 + t106) * t302;
t264 = Ifges(4,5) * t325;
t254 = Ifges(3,5) * t307;
t249 = Ifges(4,1) * t273 + t359;
t247 = Ifges(4,2) * t277 + t360;
t243 = (Ifges(4,1) * t277 - t360) * qJD(3);
t241 = (-Ifges(4,2) * t273 + t359) * qJD(3);
t238 = (mrSges(4,1) * t273 + mrSges(4,2) * t277) * qJD(3);
t214 = pkin(11) + t221;
t199 = -mrSges(4,1) * t334 - t224 * mrSges(4,3);
t198 = mrSges(4,2) * t334 + t223 * mrSges(4,3);
t195 = Ifges(5,1) * t234 + Ifges(5,4) * t233;
t194 = Ifges(5,4) * t234 + Ifges(5,2) * t233;
t193 = -mrSges(5,1) * t233 + mrSges(5,2) * t234;
t178 = mrSges(4,1) * t308 - mrSges(4,3) * t197;
t177 = -mrSges(4,2) * t308 + mrSges(4,3) * t196;
t175 = Ifges(4,1) * t224 + Ifges(4,4) * t223 - Ifges(4,5) * t334;
t174 = Ifges(4,4) * t224 + Ifges(4,2) * t223 - Ifges(4,6) * t334;
t166 = -mrSges(5,1) * t334 - t181 * mrSges(5,3);
t165 = mrSges(5,2) * t334 + t180 * mrSges(5,3);
t159 = -mrSges(4,1) * t196 + mrSges(4,2) * t197;
t157 = Ifges(5,1) * t191 - Ifges(5,4) * t192;
t156 = Ifges(5,4) * t191 - Ifges(5,2) * t192;
t155 = mrSges(5,1) * t192 + mrSges(5,2) * t191;
t147 = Ifges(6,1) * t186 - Ifges(6,4) * t185;
t146 = Ifges(6,4) * t186 - Ifges(6,2) * t185;
t145 = mrSges(6,1) * t185 + mrSges(6,2) * t186;
t144 = Ifges(4,1) * t197 + Ifges(4,4) * t196 + Ifges(4,5) * t308;
t143 = Ifges(4,4) * t197 + Ifges(4,2) * t196 + Ifges(4,6) * t308;
t135 = t295 * t186;
t130 = -mrSges(5,1) * t180 + mrSges(5,2) * t181;
t124 = Ifges(5,1) * t181 + Ifges(5,4) * t180 - Ifges(5,5) * t334;
t123 = Ifges(5,4) * t181 + Ifges(5,2) * t180 - Ifges(5,6) * t334;
t119 = mrSges(6,2) * t334 - t128 * mrSges(6,3);
t105 = Ifges(7,3) * t185 + (Ifges(7,5) * t275 - t356) * t186;
t100 = -mrSges(5,2) * t308 + mrSges(5,3) * t113;
t99 = mrSges(5,1) * t308 - mrSges(5,3) * t112;
t93 = Ifges(6,1) * t153 - Ifges(6,4) * t152;
t92 = Ifges(6,4) * t153 - Ifges(6,2) * t152;
t87 = mrSges(6,1) * t128 + mrSges(6,2) * t129;
t85 = -mrSges(7,2) * t152 - mrSges(7,3) * t287;
t84 = mrSges(7,1) * t152 + mrSges(7,3) * t286;
t82 = Ifges(6,1) * t129 - Ifges(6,4) * t128 - Ifges(6,5) * t334;
t81 = Ifges(6,4) * t129 - Ifges(6,2) * t128 - Ifges(6,6) * t334;
t72 = -mrSges(5,1) * t113 + mrSges(5,2) * t112;
t62 = mrSges(7,1) * t287 - mrSges(7,2) * t286;
t61 = Ifges(5,1) * t112 + Ifges(5,4) * t113 + Ifges(5,5) * t308;
t60 = Ifges(5,4) * t112 + Ifges(5,2) * t113 + Ifges(5,6) * t308;
t58 = mrSges(6,1) * t308 - mrSges(6,3) * t68;
t57 = -mrSges(6,2) * t308 - mrSges(6,3) * t67;
t52 = -Ifges(7,5) * t288 + Ifges(7,6) * t117 + Ifges(7,3) * t128;
t29 = Ifges(6,1) * t68 - Ifges(6,4) * t67 + Ifges(6,5) * t308;
t28 = Ifges(6,4) * t68 - Ifges(6,2) * t67 + Ifges(6,6) * t308;
t25 = -mrSges(7,2) * t67 + mrSges(7,3) * t44;
t24 = mrSges(7,1) * t67 - mrSges(7,3) * t43;
t15 = -mrSges(7,1) * t44 + mrSges(7,2) * t43;
t1 = [(t254 - 0.2e1 * t349 - 0.2e1 * t350) * t270 + (t12 - t28) * t128 + 0.2e1 * m(4) * (t121 * t172 + t122 * t171 + t210 * t219) + 0.2e1 * m(3) * (t218 * t230 - t219 * t229) + (-t81 + t52) * t67 + t224 * t144 + t223 * t143 + 0.2e1 * t210 * t159 + t196 * t174 + t197 * t175 + 0.2e1 * t121 * t198 + 0.2e1 * t122 * t199 + 0.2e1 * t183 * t72 + 0.2e1 * t172 * t177 + 0.2e1 * t171 * t178 + t180 * t60 + t181 * t61 + 0.2e1 * t32 * t165 + 0.2e1 * t33 * t166 + 0.2e1 * t169 * t130 + t129 * t29 + 0.2e1 * t134 * t30 + t117 * t13 + 0.2e1 * t9 * t119 + 0.2e1 * t8 * t120 + t113 * t123 + t112 * t124 + 0.2e1 * t89 * t99 + 0.2e1 * t90 * t100 + 0.2e1 * t88 * t87 + t68 * t82 + 0.2e1 * t2 * t76 + 0.2e1 * t3 * t77 + 0.2e1 * t5 * t73 + 0.2e1 * t38 * t57 + 0.2e1 * t37 * t58 + t44 * t53 + t43 * t54 + 0.2e1 * t35 * t15 + 0.2e1 * t16 * t24 + 0.2e1 * t17 * t25 + (mrSges(3,3) * t274 * t380 + (0.2e1 * mrSges(3,3) * t218 - t289 - t388) * t278 + ((t229 * t384 + Ifges(3,5) * t270 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t278) * t320) * t278 + (t230 * t384 + Ifges(4,5) * t224 + Ifges(5,5) * t181 + Ifges(6,5) * t129 - 0.2e1 * Ifges(3,6) * t270 + Ifges(4,6) * t223 + Ifges(5,6) * t180 - Ifges(6,6) * t128 + (-mrSges(3,1) * pkin(1) - Ifges(3,4) * t274) * t320 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3) - t398) * t334) * t274) * qJD(2)) * t269 + (-mrSges(4,1) * t223 + mrSges(4,2) * t224) * t380 + (t16 * t3 + t17 * t2 + t35 * t5) * t385 + (t134 * t88 + t37 * t8 + t38 * t9) * t386 + (t169 * t183 + t32 * t90 + t33 * t89) * t387 - t288 * t14; ((Ifges(4,5) * t273 / 0.2e1 + Ifges(4,6) * t277 / 0.2e1 - Ifges(3,6) + Ifges(6,5) * t186 / 0.2e1 - Ifges(6,6) * t185 / 0.2e1 + Ifges(5,5) * t374 + Ifges(5,6) * t375) * t327 - (-Ifges(4,6) * t326 + t264 + t390) * t278 / 0.2e1) * t269 - t349 - t350 + (-t38 * mrSges(6,3) + t52 / 0.2e1 - t81 / 0.2e1) * t152 + m(5) * (t161 * t90 + t162 * t89 + t169 * t262 + t183 * t265 + t200 * t33 - t32 * t397) - t397 * t100 + ((t175 / 0.2e1 - t171 * mrSges(4,3) - pkin(9) * t199) * t277 + (-t174 / 0.2e1 + pkin(3) * t130 - t172 * mrSges(4,3) - pkin(9) * t198) * t273) * qJD(3) + (t54 * t366 + t53 * t368 - t37 * mrSges(6,3) + t82 / 0.2e1) * t153 + (t14 * t366 + t13 * t368 - t8 * mrSges(6,3) + t29 / 0.2e1 + (-t275 * t53 / 0.2e1 + t54 * t368) * qJD(6)) * t186 + t347 * t70 + (-t9 * mrSges(6,3) + t12 / 0.2e1 - t28 / 0.2e1) * t185 + t254 + (pkin(9) * t177 + t121 * mrSges(4,3) - t219 * mrSges(4,1) + t143 / 0.2e1) * t277 + (-pkin(9) * t178 - t122 * mrSges(4,3) + t219 * mrSges(4,2) + t144 / 0.2e1) * t273 + ((t121 * t277 - t122 * t273 + (-t171 * t277 - t172 * t273) * qJD(3)) * pkin(9) - pkin(2) * t219) * m(4) + t262 * t72 + t196 * t247 / 0.2e1 + t197 * t249 / 0.2e1 + t210 * t238 + t223 * t241 / 0.2e1 + t224 * t243 / 0.2e1 + t209 * t30 + t169 * t193 + t113 * t194 / 0.2e1 + t112 * t195 / 0.2e1 + t200 * t99 + t191 * t124 / 0.2e1 - t192 * t123 / 0.2e1 + t180 * t156 / 0.2e1 + t181 * t157 / 0.2e1 + t182 * t87 + t183 * t155 + t161 * t165 + t162 * t166 - pkin(2) * t159 + t2 * t139 + t3 * t140 + t88 * t145 + t68 * t147 / 0.2e1 + t129 * t93 / 0.2e1 + t134 * t91 + t5 * t135 + t127 * t57 + t71 * t119 + t43 * t107 / 0.2e1 + t16 * t84 + t17 * t85 + t20 * t76 + t21 * t77 + t79 * t24 + t80 * t25 + t35 * t62 + m(6) * (-t126 * t8 + t127 * t9 + t134 * t182 + t209 * t88 - t37 * t70 + t38 * t71) + m(7) * (t126 * t5 + t16 * t21 + t17 * t20 + t2 * t80 + t3 * t79 + t35 * t70) + (t15 - t58) * t126 + (-t146 / 0.2e1 + t105 / 0.2e1) * t67 + (t49 / 0.2e1 - t92 / 0.2e1) * t128 + t61 * t374 + t60 * t375 + t50 * t376 + t106 * t378 - t288 * t51 / 0.2e1 + (-t191 * t89 - t192 * t90 + t233 * t32 - t234 * t33) * mrSges(5,3); (-t161 * t397 + t162 * t200 + t262 * t265) * t387 + 0.2e1 * (t161 * t233 - t162 * t234 - t191 * t200 + t192 * t397) * mrSges(5,3) + (t127 * t383 + t105 - t146) * t152 + (t383 * t71 + t49 - t92) * t185 + (mrSges(6,3) * t381 - t106 * t271 + t107 * t275 + t147) * t153 + (mrSges(6,3) * t382 - t271 * t50 + t275 * t51 + t93 + (-t106 * t275 - t107 * t271) * qJD(6)) * t186 + t277 * t241 + t273 * t243 + 0.2e1 * t262 * t155 + t234 * t157 - 0.2e1 * pkin(2) * t238 + t233 * t156 + 0.2e1 * t209 * t91 - t192 * t194 + t191 * t195 + 0.2e1 * t182 * t145 + 0.2e1 * t20 * t139 + 0.2e1 * t21 * t140 + 0.2e1 * t79 * t84 + 0.2e1 * t80 * t85 + (t277 * t249 + (0.2e1 * pkin(3) * t193 - t247) * t273) * qJD(3) + t62 * t381 + t135 * t382 + (t20 * t80 + t21 * t79 + t354) * t385 + (t127 * t71 + t182 * t209 + t354) * t386; t280 + t220 * t58 + t221 * t57 + t213 * t15 - t121 * mrSges(4,2) + t122 * mrSges(4,1) + m(6) * (-t215 * t37 + t217 * t38 + t220 * t8 + t221 * t9) + m(7) * (-t16 * t341 + t17 * t340 + t213 * t5 + t215 * t35) + t347 * t215 + (t119 + t292) * t217 + (m(5) * (t272 * t32 + t276 * t33 + t323 * t90 - t324 * t89) + t276 * t99 + t272 * t100 + t165 * t323 - t166 * t324) * pkin(3) + (m(7) * (-t361 + t362 + t391) + t275 * t25 - t271 * t24 + t392) * t214 + t284 * mrSges(7,3) + t388; t279 + t283 * mrSges(7,3) + m(7) * (t213 * t70 + t80 * t340 - t79 * t341 + t345) + t264 + (-Ifges(4,6) * t273 + (-mrSges(4,1) * t277 + mrSges(4,2) * t273) * pkin(9)) * qJD(3) + m(6) * (t127 * t217 - t220 * t70 + t221 * t71 + t345) + t213 * t62 + t215 * t135 + (-t152 * t221 - t153 * t220 - t185 * t217 + t186 * t215) * mrSges(6,3) + t291 * t217 + (m(5) * (t161 * t272 + t162 * t276 + (-t200 * t272 - t276 * t397) * qJD(4)) + (-t276 * t191 - t272 * t192 + (t233 * t276 + t234 * t272) * qJD(4)) * mrSges(5,3)) * pkin(3) + (m(7) * (-t352 + t353 + t393) + t275 * t85 - t271 * t84 + t394) * t214; 0.2e1 * t342 - 0.2e1 * t217 * mrSges(6,2) + (t213 * t215 + t214 * t299) * t385 + (-t215 * t220 + t217 * t221) * t386 + t285 + 0.2e1 * t396 * t215 + 0.2e1 * (-mrSges(5,1) * t272 - mrSges(5,2) * t276) * t355 + 0.2e1 * mrSges(7,3) * t299; t280 + (t268 * t9 + t346 * t8) * t379 + t25 * t338 - t24 * t339 - t3 * t348 + t57 * t363 + t58 * t309 + (m(7) * t5 + t15) * t260 + (m(7) * (t284 + t362) + t392) * t259 + t391 * mrSges(7,3); t279 + (t268 * t71 - t346 * t70) * t379 + t85 * t338 - t84 * t339 - t21 * t348 + (m(7) * t70 + t62) * t260 + (m(7) * (t283 + t353) + t394) * t259 + t393 * mrSges(7,3) + (-t152 * t363 - t153 * t309) * mrSges(6,3); t285 + t337 + t342 + (t268 * t379 - mrSges(6,2)) * t217 + (m(7) * t260 - t346 * t379 + t396) * t215 + (-mrSges(5,1) * t324 - mrSges(5,2) * t323) * pkin(3) + (m(7) * t259 + mrSges(7,3)) * t299; t285 + 0.2e1 * t337; t275 * t24 + t271 * t25 + t292 * qJD(6) + m(7) * (t2 * t271 + t275 * t3 + (-t16 * t271 + t17 * t275) * qJD(6)) + m(6) * t88 + t30; t271 * t85 + t275 * t84 + t291 * qJD(6) + m(7) * (t20 * t271 + t21 * t275 + (-t271 * t79 + t275 * t80) * qJD(6)) + m(6) * t182 + t91; 0; 0; 0; mrSges(7,1) * t3 - mrSges(7,2) * t2 + t12; mrSges(7,1) * t21 - mrSges(7,2) * t20 + t49; t263 - t295 * t217 + (t214 * t244 - t356) * qJD(6); t263 + (t244 * t259 - t356) * qJD(6); -t237; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
