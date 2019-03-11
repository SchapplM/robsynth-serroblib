% Calculate time derivative of joint inertia matrix for
% S6RRRPRR7
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
% Datum: 2019-03-09 18:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR7_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR7_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR7_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR7_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:32:53
% EndTime: 2019-03-09 18:33:11
% DurationCPUTime: 8.18s
% Computational Cost: add. (18684->640), mult. (47456->945), div. (0->0), fcn. (48825->12), ass. (0->267)
t375 = Ifges(4,3) + Ifges(5,3);
t271 = sin(qJ(6));
t275 = cos(qJ(6));
t243 = -mrSges(7,1) * t275 + mrSges(7,2) * t271;
t371 = -mrSges(6,1) + t243;
t267 = sin(pkin(12));
t269 = cos(pkin(12));
t273 = sin(qJ(3));
t276 = cos(qJ(3));
t233 = t267 * t276 + t269 * t273;
t223 = t233 * qJD(3);
t321 = t269 * t276;
t232 = -t267 * t273 + t321;
t224 = t232 * qJD(3);
t272 = sin(qJ(5));
t347 = cos(qJ(5));
t288 = t232 * t347 - t272 * t233;
t134 = qJD(5) * t288 - t272 * t223 + t224 * t347;
t183 = t272 * t232 + t233 * t347;
t310 = qJD(6) * t275;
t287 = t271 * t134 + t183 * t310;
t311 = qJD(6) * t271;
t328 = t134 * t275;
t286 = t183 * t311 - t328;
t342 = -qJ(4) - pkin(9);
t299 = qJD(3) * t342;
t172 = t269 * (qJD(4) * t276 + t273 * t299) + t267 * (-qJD(4) * t273 + t276 * t299);
t156 = -pkin(10) * t223 + t172;
t244 = t342 * t273;
t234 = t267 * t244;
t245 = t342 * t276;
t192 = -t269 * t245 + t234;
t173 = pkin(10) * t232 + t192;
t171 = -t233 * qJD(4) + (t321 * t342 - t234) * qJD(3);
t282 = -t224 * pkin(10) + t171;
t191 = t244 * t269 + t245 * t267;
t284 = -pkin(10) * t233 + t191;
t283 = t347 * t284;
t312 = qJD(5) * t272;
t67 = qJD(5) * t283 + t156 * t347 - t173 * t312 + t272 * t282;
t261 = -pkin(3) * t276 - pkin(2);
t203 = -pkin(4) * t232 + t261;
t112 = -pkin(5) * t288 - pkin(11) * t183 + t203;
t118 = t173 * t347 + t272 * t284;
t77 = t112 * t275 - t118 * t271;
t135 = qJD(5) * t183 + t223 * t347 + t272 * t224;
t314 = qJD(3) * t273;
t264 = pkin(3) * t314;
t196 = pkin(4) * t223 + t264;
t84 = pkin(5) * t135 - pkin(11) * t134 + t196;
t19 = qJD(6) * t77 + t271 * t84 + t275 * t67;
t78 = t112 * t271 + t118 * t275;
t20 = -qJD(6) * t78 - t271 * t67 + t275 * t84;
t374 = t19 * t275 - t20 * t271;
t268 = sin(pkin(6));
t277 = cos(qJ(2));
t322 = t268 * t277;
t270 = cos(pkin(6));
t274 = sin(qJ(2));
t323 = t268 * t274;
t225 = t270 * t276 - t273 * t323;
t226 = t270 * t273 + t276 * t323;
t177 = t225 * t267 + t226 * t269;
t230 = t270 * t274 * pkin(1) + pkin(8) * t322;
t209 = pkin(9) * t270 + t230;
t210 = (-pkin(2) * t277 - pkin(9) * t274 - pkin(1)) * t268;
t165 = -t273 * t209 + t276 * t210;
t143 = -pkin(3) * t322 - t226 * qJ(4) + t165;
t166 = t276 * t209 + t273 * t210;
t153 = qJ(4) * t225 + t166;
t94 = t269 * t143 - t267 * t153;
t74 = -pkin(4) * t322 - t177 * pkin(10) + t94;
t176 = t225 * t269 - t226 * t267;
t95 = t267 * t143 + t269 * t153;
t82 = pkin(10) * t176 + t95;
t341 = t272 * t74 + t347 * t82;
t29 = -pkin(11) * t322 + t341;
t122 = t272 * t176 + t177 * t347;
t256 = pkin(8) * t323;
t346 = pkin(1) * t277;
t208 = t256 + (-pkin(2) - t346) * t270;
t178 = -t225 * pkin(3) + t208;
t126 = -t176 * pkin(4) + t178;
t289 = t176 * t347 - t272 * t177;
t54 = -pkin(5) * t289 - t122 * pkin(11) + t126;
t16 = -t271 * t29 + t275 * t54;
t316 = qJD(2) * t268;
t304 = t277 * t316;
t189 = -qJD(3) * t226 - t273 * t304;
t190 = qJD(3) * t225 + t276 * t304;
t151 = t189 * t269 - t190 * t267;
t216 = t230 * qJD(2);
t163 = -t189 * pkin(3) + t216;
t100 = -t151 * pkin(4) + t163;
t152 = t189 * t267 + t190 * t269;
t64 = qJD(5) * t289 + t272 * t151 + t152 * t347;
t65 = qJD(5) * t122 - t151 * t347 + t272 * t152;
t23 = t65 * pkin(5) - t64 * pkin(11) + t100;
t315 = qJD(2) * t274;
t305 = t268 * t315;
t306 = t347 * t74;
t214 = (pkin(2) * t274 - pkin(9) * t277) * t316;
t229 = t270 * t346 - t256;
t215 = t229 * qJD(2);
t114 = -qJD(3) * t166 + t276 * t214 - t215 * t273;
t90 = pkin(3) * t305 - qJ(4) * t190 - qJD(4) * t226 + t114;
t313 = qJD(3) * t276;
t113 = -t209 * t314 + t210 * t313 + t273 * t214 + t276 * t215;
t96 = qJ(4) * t189 + qJD(4) * t225 + t113;
t49 = -t267 * t96 + t269 * t90;
t33 = pkin(4) * t305 - pkin(10) * t152 + t49;
t50 = t267 * t90 + t269 * t96;
t35 = pkin(10) * t151 + t50;
t8 = qJD(5) * t306 + t272 * t33 - t312 * t82 + t347 * t35;
t5 = pkin(11) * t305 + t8;
t2 = qJD(6) * t16 + t23 * t271 + t275 * t5;
t17 = t271 * t54 + t275 * t29;
t3 = -qJD(6) * t17 + t23 * t275 - t271 * t5;
t373 = t2 * t275 - t3 * t271;
t350 = t271 / 0.2e1;
t349 = t275 / 0.2e1;
t301 = -t311 / 0.2e1;
t260 = pkin(3) * t269 + pkin(4);
t345 = pkin(3) * t267;
t220 = t260 * t347 - t272 * t345;
t206 = t220 * qJD(5);
t372 = t206 * mrSges(6,2);
t298 = (t271 ^ 2 + t275 ^ 2) * t206;
t221 = t272 * t260 + t347 * t345;
t370 = Ifges(4,5) * t313 + Ifges(5,5) * t224 - Ifges(5,6) * t223;
t9 = -qJD(5) * t341 - t272 * t35 + t33 * t347;
t369 = Ifges(4,5) * t190 + Ifges(5,5) * t152 + Ifges(4,6) * t189 + Ifges(5,6) * t151 + t305 * t375;
t368 = 2 * m(5);
t367 = 2 * m(6);
t366 = 2 * m(7);
t365 = -2 * mrSges(3,3);
t364 = -2 * mrSges(6,3);
t68 = qJD(5) * t118 + t272 * t156 - t282 * t347;
t363 = 0.2e1 * t68;
t117 = t173 * t272 - t283;
t362 = 0.2e1 * t117;
t361 = 0.2e1 * t216;
t290 = -t275 * t122 + t271 * t322;
t41 = qJD(6) * t290 - t271 * t64 + t275 * t305;
t360 = t41 / 0.2e1;
t108 = -t271 * t122 - t275 * t322;
t359 = t108 / 0.2e1;
t358 = t232 / 0.2e1;
t357 = t233 / 0.2e1;
t262 = Ifges(7,5) * t310;
t356 = Ifges(7,6) * t301 + t262 / 0.2e1;
t338 = Ifges(7,4) * t271;
t296 = Ifges(7,1) * t275 - t338;
t241 = t296 * qJD(6);
t355 = t241 / 0.2e1;
t354 = Ifges(7,5) * t350 + Ifges(7,6) * t349;
t337 = Ifges(7,4) * t275;
t249 = Ifges(7,1) * t271 + t337;
t352 = t249 / 0.2e1;
t351 = -t271 / 0.2e1;
t340 = Ifges(4,4) * t273;
t339 = Ifges(4,4) * t276;
t336 = Ifges(7,6) * t271;
t335 = t117 * t68;
t332 = t215 * mrSges(3,2);
t331 = t216 * mrSges(3,1);
t111 = -mrSges(6,1) * t322 - t122 * mrSges(6,3);
t69 = -mrSges(7,1) * t108 - mrSges(7,2) * t290;
t330 = t69 - t111;
t207 = t221 * qJD(5);
t329 = t117 * t207;
t327 = t183 * t271;
t326 = t183 * t275;
t325 = t206 * t271;
t324 = t206 * t275;
t319 = Ifges(6,5) * t134 - Ifges(6,6) * t135;
t309 = 0.2e1 * t268;
t40 = qJD(6) * t108 + t271 * t305 + t275 * t64;
t12 = Ifges(7,5) * t40 + Ifges(7,6) * t41 + Ifges(7,3) * t65;
t307 = Ifges(6,5) * t64 - Ifges(6,6) * t65 + Ifges(6,3) * t305;
t26 = t65 * mrSges(6,1) + t64 * mrSges(6,2);
t300 = t310 / 0.2e1;
t97 = -t151 * mrSges(5,1) + t152 * mrSges(5,2);
t179 = t223 * mrSges(5,1) + t224 * mrSges(5,2);
t85 = t135 * mrSges(6,1) + t134 * mrSges(6,2);
t297 = mrSges(7,1) * t271 + mrSges(7,2) * t275;
t295 = -Ifges(7,2) * t271 + t337;
t75 = mrSges(7,2) * t289 + mrSges(7,3) * t108;
t76 = -mrSges(7,1) * t289 + mrSges(7,3) * t290;
t294 = -t271 * t76 + t275 * t75;
t138 = mrSges(7,2) * t288 - mrSges(7,3) * t327;
t139 = -mrSges(7,1) * t288 - mrSges(7,3) * t326;
t293 = t138 * t275 - t139 * t271;
t30 = -t272 * t82 + t306;
t239 = t295 * qJD(6);
t247 = Ifges(7,2) * t275 + t338;
t285 = t275 * t239 + t271 * t241 - t247 * t311 + t249 * t310;
t44 = -Ifges(7,5) * t286 - Ifges(7,6) * t287 + Ifges(7,3) * t135;
t21 = mrSges(7,1) * t65 - mrSges(7,3) * t40;
t22 = -mrSges(7,2) * t65 + mrSges(7,3) * t41;
t281 = m(7) * (-t16 * t310 - t17 * t311 + t373) + t275 * t22 - t271 * t21 - t75 * t311 - t76 * t310;
t70 = mrSges(7,1) * t135 + mrSges(7,3) * t286;
t71 = -mrSges(7,2) * t135 - mrSges(7,3) * t287;
t280 = m(7) * (-t310 * t77 - t311 * t78 + t374) + t275 * t71 - t271 * t70 - t138 * t311 - t139 * t310;
t13 = Ifges(7,4) * t40 + Ifges(7,2) * t41 + Ifges(7,6) * t65;
t14 = Ifges(7,1) * t40 + Ifges(7,4) * t41 + Ifges(7,5) * t65;
t236 = t297 * qJD(6);
t28 = pkin(5) * t322 - t30;
t52 = -Ifges(7,4) * t290 + Ifges(7,2) * t108 - Ifges(7,6) * t289;
t53 = -Ifges(7,1) * t290 + Ifges(7,4) * t108 - Ifges(7,5) * t289;
t6 = -pkin(5) * t305 - t9;
t279 = t9 * mrSges(6,1) - t8 * mrSges(6,2) - t290 * t355 - t289 * t356 + t13 * t349 + t14 * t350 + t28 * t236 + t239 * t359 + t6 * t243 + t247 * t360 + t53 * t300 + t52 * t301 + t40 * t352 + t65 * t354 + t307 + ((-t16 * t275 - t17 * t271) * qJD(6) + t373) * mrSges(7,3);
t104 = -Ifges(7,6) * t288 + t183 * t295;
t105 = -Ifges(7,5) * t288 + t183 * t296;
t45 = -Ifges(7,4) * t286 - Ifges(7,2) * t287 + Ifges(7,6) * t135;
t46 = -Ifges(7,1) * t286 - Ifges(7,4) * t287 + Ifges(7,5) * t135;
t278 = t117 * t236 + t328 * t352 + t135 * t354 - t239 * t327 / 0.2e1 + t326 * t355 - t288 * t356 + t46 * t350 + t45 * t349 + t105 * t300 - t67 * mrSges(6,2) + t319 + t371 * t68 - t287 * t247 / 0.2e1 + (t183 * t249 + t104) * t301 + ((-t271 * t78 - t275 * t77) * qJD(6) + t374) * mrSges(7,3);
t255 = Ifges(3,5) * t304;
t250 = Ifges(4,1) * t273 + t339;
t248 = Ifges(4,2) * t276 + t340;
t242 = (Ifges(4,1) * t276 - t340) * qJD(3);
t240 = (-Ifges(4,2) * t273 + t339) * qJD(3);
t237 = (mrSges(4,1) * t273 + mrSges(4,2) * t276) * qJD(3);
t213 = pkin(11) + t221;
t212 = -pkin(5) - t220;
t194 = -mrSges(4,1) * t322 - t226 * mrSges(4,3);
t193 = mrSges(4,2) * t322 + t225 * mrSges(4,3);
t186 = Ifges(5,1) * t233 + Ifges(5,4) * t232;
t185 = Ifges(5,4) * t233 + Ifges(5,2) * t232;
t184 = -mrSges(5,1) * t232 + mrSges(5,2) * t233;
t181 = Ifges(5,1) * t224 - Ifges(5,4) * t223;
t180 = Ifges(5,4) * t224 - Ifges(5,2) * t223;
t175 = mrSges(4,1) * t305 - mrSges(4,3) * t190;
t174 = -mrSges(4,2) * t305 + mrSges(4,3) * t189;
t169 = Ifges(4,1) * t226 + Ifges(4,4) * t225 - Ifges(4,5) * t322;
t168 = Ifges(4,4) * t226 + Ifges(4,2) * t225 - Ifges(4,6) * t322;
t160 = -mrSges(5,1) * t322 - t177 * mrSges(5,3);
t159 = mrSges(5,2) * t322 + t176 * mrSges(5,3);
t155 = -mrSges(4,1) * t189 + mrSges(4,2) * t190;
t145 = Ifges(4,1) * t190 + Ifges(4,4) * t189 + Ifges(4,5) * t305;
t144 = Ifges(4,4) * t190 + Ifges(4,2) * t189 + Ifges(4,6) * t305;
t142 = Ifges(6,1) * t183 + Ifges(6,4) * t288;
t141 = Ifges(6,4) * t183 + Ifges(6,2) * t288;
t140 = -mrSges(6,1) * t288 + mrSges(6,2) * t183;
t129 = mrSges(5,1) * t305 - mrSges(5,3) * t152;
t128 = -mrSges(5,2) * t305 + mrSges(5,3) * t151;
t127 = t297 * t183;
t124 = -mrSges(5,1) * t176 + mrSges(5,2) * t177;
t116 = Ifges(5,1) * t177 + Ifges(5,4) * t176 - Ifges(5,5) * t322;
t115 = Ifges(5,4) * t177 + Ifges(5,2) * t176 - Ifges(5,6) * t322;
t110 = mrSges(6,2) * t322 + mrSges(6,3) * t289;
t103 = -Ifges(7,3) * t288 + (Ifges(7,5) * t275 - t336) * t183;
t92 = Ifges(5,1) * t152 + Ifges(5,4) * t151 + Ifges(5,5) * t305;
t91 = Ifges(5,4) * t152 + Ifges(5,2) * t151 + Ifges(5,6) * t305;
t87 = Ifges(6,1) * t134 - Ifges(6,4) * t135;
t86 = Ifges(6,4) * t134 - Ifges(6,2) * t135;
t83 = -mrSges(6,1) * t289 + mrSges(6,2) * t122;
t81 = Ifges(6,1) * t122 + Ifges(6,4) * t289 - Ifges(6,5) * t322;
t80 = Ifges(6,4) * t122 + Ifges(6,2) * t289 - Ifges(6,6) * t322;
t58 = -mrSges(6,2) * t305 - mrSges(6,3) * t65;
t57 = mrSges(6,1) * t305 - mrSges(6,3) * t64;
t56 = mrSges(7,1) * t287 - mrSges(7,2) * t286;
t51 = -Ifges(7,5) * t290 + Ifges(7,6) * t108 - Ifges(7,3) * t289;
t25 = Ifges(6,1) * t64 - Ifges(6,4) * t65 + Ifges(6,5) * t305;
t24 = Ifges(6,4) * t64 - Ifges(6,2) * t65 + Ifges(6,6) * t305;
t15 = -mrSges(7,1) * t41 + mrSges(7,2) * t40;
t1 = [(-mrSges(4,1) * t225 + mrSges(4,2) * t226) * t361 + (t16 * t3 + t17 * t2 + t28 * t6) * t366 + (t163 * t178 + t49 * t94 + t50 * t95) * t368 - t290 * t14 - (t12 - t24) * t289 + (mrSges(3,3) * t274 * t361 + (0.2e1 * mrSges(3,3) * t215 - t307 - t369) * t277 + ((t229 * t365 + Ifges(3,5) * t270 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t277) * t309) * t277 + (t230 * t365 + Ifges(4,5) * t226 + Ifges(5,5) * t177 + Ifges(6,5) * t122 - 0.2e1 * Ifges(3,6) * t270 + Ifges(4,6) * t225 + Ifges(5,6) * t176 + Ifges(6,6) * t289 + (-mrSges(3,1) * pkin(1) - Ifges(3,4) * t274) * t309 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(6,3) - t375) * t322) * t274) * qJD(2)) * t268 + t226 * t145 + t225 * t144 + 0.2e1 * t208 * t155 + t189 * t168 + t190 * t169 + 0.2e1 * t113 * t193 + 0.2e1 * t114 * t194 + 0.2e1 * t166 * t174 + 0.2e1 * t165 * t175 + t176 * t91 + t177 * t92 + 0.2e1 * t178 * t97 + 0.2e1 * t163 * t124 + t151 * t115 + t152 * t116 + 0.2e1 * t50 * t159 + 0.2e1 * t49 * t160 + t122 * t25 + 0.2e1 * t126 * t26 + 0.2e1 * t95 * t128 + 0.2e1 * t94 * t129 + 0.2e1 * t8 * t110 + 0.2e1 * t9 * t111 + 0.2e1 * t100 * t83 + t108 * t13 + 0.2e1 * t2 * t75 + 0.2e1 * t3 * t76 + t64 * t81 + 0.2e1 * t6 * t69 + 0.2e1 * t30 * t57 + t41 * t52 + t40 * t53 + 0.2e1 * t28 * t15 + 0.2e1 * t16 * t21 + 0.2e1 * t17 * t22 + (t100 * t126 + t30 * t9 + t341 * t8) * t367 + 0.2e1 * t341 * t58 + (-t80 + t51) * t65 + (t255 - 0.2e1 * t331 - 0.2e1 * t332) * t270 + 0.2e1 * m(3) * (t215 * t230 - t216 * t229) + 0.2e1 * m(4) * (t113 * t166 + t114 * t165 + t208 * t216); (-t223 * t95 - t224 * t94 + t232 * t50 - t233 * t49) * mrSges(5,3) + t92 * t357 + t91 * t358 + t45 * t359 + t104 * t360 + t255 + ((t113 * t276 - t114 * t273 + (-t165 * t276 - t166 * t273) * qJD(3)) * pkin(9) - pkin(2) * t216) * m(4) + (t15 - t57) * t117 + ((t169 / 0.2e1 - pkin(9) * t194 - t165 * mrSges(4,3)) * t276 + (-t168 / 0.2e1 - pkin(9) * t193 - t166 * mrSges(4,3) + pkin(3) * t124) * t273) * qJD(3) - t290 * t46 / 0.2e1 - (t44 / 0.2e1 - t86 / 0.2e1) * t289 - (t12 / 0.2e1 - t24 / 0.2e1 - t8 * mrSges(6,3)) * t288 + ((Ifges(5,5) * t357 + Ifges(5,6) * t358 + Ifges(6,5) * t183 / 0.2e1 + Ifges(6,6) * t288 / 0.2e1 - Ifges(3,6) + Ifges(4,5) * t273 / 0.2e1 + Ifges(4,6) * t276 / 0.2e1) * t315 - (-Ifges(4,6) * t314 + t319 + t370) * t277 / 0.2e1) * t268 + t261 * t97 + t208 * t237 + t225 * t240 / 0.2e1 + t226 * t242 / 0.2e1 + t189 * t248 / 0.2e1 + t190 * t250 / 0.2e1 - t223 * t115 / 0.2e1 + t224 * t116 / 0.2e1 + t203 * t26 + t152 * t186 / 0.2e1 + t191 * t129 + t192 * t128 + t196 * t83 + t163 * t184 + t151 * t185 / 0.2e1 + t178 * t179 + t176 * t180 / 0.2e1 + t177 * t181 / 0.2e1 + t171 * t160 + t172 * t159 - pkin(2) * t155 + t3 * t139 + t100 * t140 + t64 * t142 / 0.2e1 + t2 * t138 + t122 * t87 / 0.2e1 + t126 * t85 + t6 * t127 + t118 * t58 + t67 * t110 + t40 * t105 / 0.2e1 + t19 * t75 + t20 * t76 + t77 * t21 + t78 * t22 + t16 * t70 + t17 * t71 + t28 * t56 + (-t216 * mrSges(4,1) + t144 / 0.2e1 + pkin(9) * t174 + t113 * mrSges(4,3)) * t276 + (t216 * mrSges(4,2) + t145 / 0.2e1 - pkin(9) * t175 - t114 * mrSges(4,3)) * t273 + m(7) * (t117 * t6 + t16 * t20 + t17 * t19 + t2 * t78 + t28 * t68 + t3 * t77) + (t51 / 0.2e1 - t80 / 0.2e1 - t341 * mrSges(6,3)) * t135 + m(6) * (t100 * t203 - t117 * t9 + t118 * t8 + t126 * t196 - t30 * t68 + t341 * t67) + (-t141 / 0.2e1 + t103 / 0.2e1) * t65 - t331 - t332 + m(5) * (t163 * t261 + t171 * t94 + t172 * t95 + t178 * t264 + t191 * t49 + t192 * t50) + t330 * t68 + (t81 / 0.2e1 + t53 * t349 + t52 * t351 - t30 * mrSges(6,3)) * t134 + (t25 / 0.2e1 + t14 * t349 + t13 * t351 - t9 * mrSges(6,3) + (t53 * t351 - t275 * t52 / 0.2e1) * qJD(6)) * t183; t56 * t362 + t127 * t363 + (t19 * t78 + t20 * t77 + t335) * t366 + (t118 * t67 + t196 * t203 + t335) * t367 + (t171 * t191 + t172 * t192 + t261 * t264) * t368 - (t364 * t67 + t44 - t86) * t288 + t276 * t240 + t273 * t242 + 0.2e1 * t261 * t179 - 0.2e1 * pkin(2) * t237 + t232 * t180 + t233 * t181 - t223 * t185 + t224 * t186 + 0.2e1 * t203 * t85 + 0.2e1 * t196 * t140 + 0.2e1 * t20 * t139 + 0.2e1 * t19 * t138 + 0.2e1 * t77 * t70 + 0.2e1 * t78 * t71 + 0.2e1 * (-t171 * t233 + t172 * t232 - t191 * t224 - t192 * t223) * mrSges(5,3) + (t276 * t250 + (0.2e1 * pkin(3) * t184 - t248) * t273) * qJD(3) + (mrSges(6,3) * t362 - t104 * t271 + t105 * t275 + t142) * t134 + (mrSges(6,3) * t363 - t271 * t45 + t275 * t46 + t87 + (-t104 * t275 - t105 * t271) * qJD(6)) * t183 + (t118 * t364 + t103 - t141) * t135; t369 + (m(5) * (t267 * t50 + t269 * t49) + t269 * t129 + t267 * t128) * pkin(3) + t281 * t213 + t279 + t220 * t57 + t221 * t58 + t212 * t15 - t113 * mrSges(4,2) + t114 * mrSges(4,1) + t49 * mrSges(5,1) - t50 * mrSges(5,2) + m(6) * (t206 * t341 - t207 * t30 + t220 * t9 + t221 * t8) + (t110 + t294) * t206 + m(7) * (-t16 * t325 + t17 * t324 + t207 * t28 + t212 * t6) + t330 * t207; (m(5) * (t171 * t269 + t172 * t267) + (-t223 * t267 - t224 * t269) * mrSges(5,3)) * pkin(3) + t280 * t213 + (-Ifges(4,6) * t273 + (-mrSges(4,1) * t276 + mrSges(4,2) * t273) * pkin(9)) * qJD(3) + t293 * t206 + t278 + m(7) * (t212 * t68 + t324 * t78 - t325 * t77 + t329) + t212 * t56 + t207 * t127 + t171 * mrSges(5,1) - t172 * mrSges(5,2) + m(6) * (t118 * t206 - t220 * t68 + t221 * t67 + t329) + (-t134 * t220 - t135 * t221 + t183 * t207 + t206 * t288) * mrSges(6,3) + t370; 0.2e1 * t212 * t236 - 0.2e1 * t372 + (t207 * t212 + t213 * t298) * t366 + (t206 * t221 - t207 * t220) * t367 + t285 + 0.2e1 * t371 * t207 + 0.2e1 * mrSges(7,3) * t298; t275 * t21 + t271 * t22 + t294 * qJD(6) + m(7) * (t2 * t271 + t275 * t3 + (-t16 * t271 + t17 * t275) * qJD(6)) + m(6) * t100 + m(5) * t163 + t97 + t26; m(5) * t264 + t271 * t71 + t275 * t70 + t293 * qJD(6) + m(7) * (t19 * t271 + t20 * t275 + (-t271 * t77 + t275 * t78) * qJD(6)) + m(6) * t196 + t85 + t179; 0; 0; (-m(7) * t6 - t15) * pkin(5) + t281 * pkin(11) + t279; (-m(7) * t68 - t56) * pkin(5) + t280 * pkin(11) + t278; -t372 + (t212 - pkin(5)) * t236 + t285 + (m(7) * pkin(11) + mrSges(7,3)) * t298 + (-m(7) * pkin(5) + t371) * t207; 0; -0.2e1 * pkin(5) * t236 + t285; mrSges(7,1) * t3 - mrSges(7,2) * t2 + t12; mrSges(7,1) * t20 - mrSges(7,2) * t19 + t44; t262 - t297 * t206 + (t213 * t243 - t336) * qJD(6); -t236; t262 + (pkin(11) * t243 - t336) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
