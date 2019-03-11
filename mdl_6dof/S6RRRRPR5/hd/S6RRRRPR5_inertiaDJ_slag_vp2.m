% Calculate time derivative of joint inertia matrix for
% S6RRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 22:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:12:04
% EndTime: 2019-03-09 22:12:18
% DurationCPUTime: 6.40s
% Computational Cost: add. (7660->511), mult. (16926->710), div. (0->0), fcn. (15357->8), ass. (0->226)
t375 = Ifges(5,1) + Ifges(6,1);
t232 = sin(qJ(3));
t233 = sin(qJ(2));
t236 = cos(qJ(3));
t237 = cos(qJ(2));
t184 = t232 * t233 - t236 * t237;
t328 = Ifges(6,4) + Ifges(5,5);
t374 = t184 * t328;
t231 = sin(qJ(4));
t235 = cos(qJ(4));
t373 = t231 ^ 2 + t235 ^ 2;
t321 = Ifges(6,5) * t231;
t323 = Ifges(5,4) * t231;
t372 = t375 * t235 + t321 - t323;
t320 = Ifges(6,5) * t235;
t322 = Ifges(5,4) * t235;
t371 = t375 * t231 - t320 + t322;
t370 = Ifges(5,6) * t235 + t231 * t328;
t369 = Ifges(6,2) + Ifges(5,3);
t186 = t232 * t237 + t233 * t236;
t355 = qJD(2) + qJD(3);
t136 = t355 * t186;
t290 = qJD(4) * t231;
t135 = t355 * t184;
t305 = t135 * t235;
t247 = t186 * t290 + t305;
t289 = qJD(4) * t235;
t273 = t186 * t289;
t306 = t135 * t231;
t248 = t273 - t306;
t368 = (-Ifges(5,4) + Ifges(6,5)) * t248 - t375 * t247 + t328 * t136;
t325 = t372 * t186 + t374;
t201 = -t235 * mrSges(5,1) + t231 * mrSges(5,2);
t367 = -mrSges(4,1) + t201;
t366 = t372 * qJD(4);
t365 = Ifges(6,6) * t290 + t328 * t289;
t316 = pkin(2) * qJD(3);
t278 = t236 * t316;
t363 = t373 * t278;
t216 = -pkin(2) * t237 - pkin(1);
t118 = pkin(3) * t184 - pkin(9) * t186 + t216;
t341 = -pkin(8) - pkin(7);
t207 = t341 * t233;
t209 = t341 * t237;
t156 = t207 * t232 - t209 * t236;
t280 = pkin(2) * qJD(2) * t233;
t68 = pkin(3) * t136 + pkin(9) * t135 + t280;
t275 = qJD(2) * t341;
t196 = t233 * t275;
t262 = t237 * t275;
t359 = t236 * t207 + t209 * t232;
t84 = qJD(3) * t359 + t236 * t196 + t232 * t262;
t21 = t118 * t289 - t156 * t290 + t231 * t68 + t235 * t84;
t22 = -t118 * t290 - t156 * t289 - t231 * t84 + t235 * t68;
t362 = t21 * t235 - t22 * t231;
t13 = t136 * qJ(5) + t184 * qJD(5) + t21;
t16 = -pkin(4) * t136 - t22;
t71 = t231 * t118 + t235 * t156;
t59 = t184 * qJ(5) + t71;
t144 = t231 * t156;
t70 = t118 * t235 - t144;
t60 = -pkin(4) * t184 - t70;
t361 = t13 * t235 + t16 * t231 + t60 * t289 - t59 * t290;
t46 = mrSges(5,1) * t248 - mrSges(5,2) * t247;
t85 = t156 * qJD(3) + t196 * t232 - t236 * t262;
t360 = m(5) * t85 + t46;
t230 = sin(qJ(6));
t234 = cos(qJ(6));
t252 = t230 * t235 - t231 * t234;
t106 = t252 * t186;
t225 = t231 * qJ(5);
t254 = t235 * pkin(4) + t225;
t340 = pkin(9) - pkin(10);
t206 = t340 * t231;
t208 = t340 * t235;
t153 = t206 * t234 - t208 * t230;
t222 = pkin(10) * t290;
t194 = -pkin(9) * t290 + t222;
t195 = qJD(4) * t208;
t82 = qJD(6) * t153 + t194 * t234 + t195 * t230;
t155 = t206 * t230 + t208 * t234;
t83 = -qJD(6) * t155 - t194 * t230 + t195 * t234;
t358 = t83 * mrSges(7,1) - t82 * mrSges(7,2);
t214 = pkin(2) * t232 + pkin(9);
t327 = -pkin(10) + t214;
t179 = t327 * t231;
t180 = t327 * t235;
t109 = t179 * t234 - t180 * t230;
t157 = -t214 * t290 + t235 * t278 + t222;
t264 = t231 * t278;
t158 = qJD(4) * t180 + t264;
t50 = qJD(6) * t109 + t157 * t234 + t158 * t230;
t110 = t179 * t230 + t180 * t234;
t51 = -qJD(6) * t110 - t157 * t230 + t158 * t234;
t357 = t51 * mrSges(7,1) - t50 * mrSges(7,2);
t10 = pkin(10) * t248 + t13;
t238 = -pkin(4) - pkin(5);
t44 = t144 + (-pkin(10) * t186 - t118) * t235 + t238 * t184;
t303 = t186 * t231;
t49 = pkin(10) * t303 + t59;
t18 = -t230 * t49 + t234 * t44;
t8 = pkin(10) * t247 + t136 * t238 - t22;
t2 = qJD(6) * t18 + t10 * t234 + t230 * t8;
t19 = t230 * t44 + t234 * t49;
t3 = -qJD(6) * t19 - t10 * t230 + t234 * t8;
t356 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t354 = qJD(4) - qJD(6);
t353 = (-mrSges(4,2) + (mrSges(6,2) + mrSges(5,3)) * t373) * t278;
t352 = 0.2e1 * m(5);
t351 = 2 * m(6);
t350 = 2 * m(7);
t349 = -2 * mrSges(4,3);
t348 = -2 * Ifges(4,4);
t251 = t230 * t231 + t234 * t235;
t133 = t354 * t251;
t134 = t354 * t252;
t72 = mrSges(7,1) * t134 + mrSges(7,2) * t133;
t347 = 0.2e1 * t72;
t346 = 0.2e1 * t85;
t141 = mrSges(7,1) * t251 - mrSges(7,2) * t252;
t345 = 0.2e1 * t141;
t344 = -0.2e1 * t359;
t259 = t231 * mrSges(6,1) - t235 * mrSges(6,3);
t188 = t259 * qJD(4);
t343 = 0.2e1 * t188;
t342 = 0.2e1 * t216;
t336 = t136 / 0.2e1;
t203 = Ifges(5,2) * t235 + t323;
t334 = -t203 / 0.2e1;
t255 = Ifges(6,3) * t231 + t320;
t92 = Ifges(6,6) * t184 + t186 * t255;
t256 = -Ifges(5,2) * t231 + t322;
t319 = Ifges(5,6) * t184;
t93 = t186 * t256 + t319;
t326 = t92 - t93;
t324 = mrSges(7,3) * t251;
t318 = Ifges(5,6) * t231;
t314 = t134 * mrSges(7,3);
t313 = t359 * t85;
t197 = -t230 * qJ(5) + t234 * t238;
t162 = t234 * qJD(5) + qJD(6) * t197;
t311 = t162 * mrSges(7,2);
t198 = t234 * qJ(5) + t230 * t238;
t163 = -t230 * qJD(5) - qJD(6) * t198;
t310 = t163 * mrSges(7,1);
t307 = qJ(5) * t235;
t304 = t359 * t232;
t302 = t186 * t235;
t300 = t231 * t236;
t298 = t235 * t236;
t120 = -mrSges(5,2) * t184 - mrSges(5,3) * t303;
t123 = -mrSges(6,2) * t303 + mrSges(6,3) * t184;
t297 = t120 + t123;
t121 = mrSges(5,1) * t184 - mrSges(5,3) * t302;
t122 = -mrSges(6,1) * t184 + mrSges(6,2) * t302;
t296 = -t121 + t122;
t295 = Ifges(7,5) * t133 - Ifges(7,6) * t134;
t294 = t363 * t214;
t293 = t363 * pkin(9);
t288 = qJD(5) * t235;
t287 = qJD(6) * t230;
t286 = qJD(6) * t234;
t285 = 0.2e1 * mrSges(7,3);
t283 = pkin(3) + t254;
t31 = t106 * t354 - t251 * t135;
t32 = t133 * t186 + t135 * t252;
t282 = -Ifges(7,5) * t31 - Ifges(7,6) * t32 + Ifges(7,3) * t136;
t281 = m(6) * t288;
t279 = t232 * t316;
t167 = pkin(4) * t290 - qJ(5) * t289 - t231 * qJD(5);
t215 = -pkin(2) * t236 - pkin(3);
t270 = -t290 / 0.2e1;
t268 = 0.2e1 * t280;
t260 = t231 * mrSges(5,1) + t235 * mrSges(5,2);
t200 = -t235 * mrSges(6,1) - t231 * mrSges(6,3);
t253 = pkin(4) * t231 - t307;
t173 = t215 - t254;
t250 = t367 * t279;
t249 = t231 * t238 + t307;
t159 = -pkin(5) * t290 - t167;
t246 = mrSges(6,2) * t289 - t230 * t314 - t286 * t324 + (-t133 * t234 - t252 * t287) * mrSges(7,3);
t56 = -t136 * mrSges(6,1) - mrSges(6,2) * t247;
t245 = t248 * Ifges(6,6) + t369 * t136 - t305 * t328 + t282;
t244 = -mrSges(6,2) * t254 - t318;
t243 = -t198 * t314 - t162 * t324 + mrSges(6,2) * t288 + (-t133 * t197 + t163 * t252) * mrSges(7,3) - t295 + t365;
t142 = -Ifges(7,4) * t252 - Ifges(7,2) * t251;
t143 = -Ifges(7,1) * t252 - Ifges(7,4) * t251;
t190 = t255 * qJD(4);
t191 = t256 * qJD(4);
t202 = -Ifges(6,3) * t235 + t321;
t73 = Ifges(7,4) * t133 - Ifges(7,2) * t134;
t74 = Ifges(7,1) * t133 - Ifges(7,4) * t134;
t242 = t133 * t143 - t134 * t142 - t251 * t73 - t252 * t74 + (t202 - t203) * t290 + t371 * t289 + (-t190 + t191) * t235 + t366 * t231;
t241 = -m(6) * t254 + t200 + t201;
t55 = mrSges(5,1) * t136 + mrSges(5,3) * t247;
t57 = -mrSges(5,2) * t136 - mrSges(5,3) * t248;
t58 = -mrSges(6,2) * t248 + mrSges(6,3) * t136;
t240 = (t57 + t58) * t235 + (-t55 + t56) * t231 + (-t231 * t297 + t235 * t296) * qJD(4) + m(6) * t361 + m(5) * (-t289 * t70 - t290 * t71 + t362);
t107 = t251 * t186;
t17 = -t249 * t135 + (t288 + (t235 * t238 - t225) * qJD(4)) * t186 - t85;
t189 = t260 * qJD(4);
t28 = -t253 * t135 + (qJD(4) * t254 - t288) * t186 + t85;
t38 = -Ifges(6,5) * t247 + Ifges(6,6) * t136 + Ifges(6,3) * t248;
t39 = -Ifges(5,4) * t247 - Ifges(5,2) * t248 + Ifges(5,6) * t136;
t52 = Ifges(7,4) * t107 - Ifges(7,2) * t106 - Ifges(7,6) * t184;
t53 = Ifges(7,1) * t107 - Ifges(7,4) * t106 - Ifges(7,5) * t184;
t6 = Ifges(7,4) * t31 + Ifges(7,2) * t32 - Ifges(7,6) * t136;
t65 = t186 * t249 + t359;
t7 = Ifges(7,1) * t31 + Ifges(7,4) * t32 - Ifges(7,5) * t136;
t81 = t186 * t253 - t359;
t239 = ((-t231 * t71 - t235 * t70) * qJD(4) + t362) * mrSges(5,3) + t361 * mrSges(6,2) + t370 * t336 + (-t191 / 0.2e1 + t190 / 0.2e1) * t303 + (t186 * t202 + t325) * t289 / 0.2e1 + t368 * t231 / 0.2e1 + t367 * t85 - t251 * t6 / 0.2e1 - t136 * (-Ifges(7,5) * t252 - Ifges(7,6) * t251) / 0.2e1 + (-t133 * t18 + t252 * t3) * mrSges(7,3) - t252 * t7 / 0.2e1 + t371 * (t186 * t270 - t305 / 0.2e1) + t273 * t334 - t19 * t314 - t2 * t324 - (t334 + t202 / 0.2e1) * t306 + t93 * t270 + (-Ifges(6,6) * t336 + t39 / 0.2e1 - t38 / 0.2e1) * t235 - t184 * t295 / 0.2e1 + t92 * t290 / 0.2e1 + (-Ifges(5,6) * t290 + t365) * t184 / 0.2e1 + t366 * t302 / 0.2e1 + t65 * t72 - t359 * t189 - t84 * mrSges(4,2) - t106 * t73 / 0.2e1 + t107 * t74 / 0.2e1 + t133 * t53 / 0.2e1 - t134 * t52 / 0.2e1 - Ifges(4,5) * t135 - Ifges(4,6) * t136 + t17 * t141 + t32 * t142 / 0.2e1 + t31 * t143 / 0.2e1 + t81 * t188 + t28 * t200;
t226 = t235 * pkin(5);
t174 = t226 + t283;
t161 = -t173 + t226;
t160 = t167 + t279;
t146 = t159 - t279;
t112 = t260 * t186;
t111 = t259 * t186;
t91 = -mrSges(7,1) * t184 - mrSges(7,3) * t107;
t90 = mrSges(7,2) * t184 - mrSges(7,3) * t106;
t61 = mrSges(7,1) * t106 + mrSges(7,2) * t107;
t45 = mrSges(6,1) * t248 + mrSges(6,3) * t247;
t25 = mrSges(7,2) * t136 + mrSges(7,3) * t32;
t24 = -mrSges(7,1) * t136 - mrSges(7,3) * t31;
t9 = -mrSges(7,1) * t32 + mrSges(7,2) * t31;
t1 = [(t21 * t71 + t22 * t70 - t313) * t352 + 0.2e1 * m(4) * (t156 * t84 + t216 * t280 - t313) + (mrSges(4,1) * t268 + t84 * t349 - (t348 - t318) * t135 + ((2 * Ifges(4,2)) + Ifges(7,3) + t369) * t136 + t245) * t184 + 0.2e1 * ((-t233 ^ 2 + t237 ^ 2) * Ifges(3,4) - pkin(1) * (mrSges(3,1) * t233 + mrSges(3,2) * t237) + (-Ifges(3,2) + Ifges(3,1)) * t233 * t237) * qJD(2) + (mrSges(4,2) * t268 + mrSges(4,3) * t346 - 0.2e1 * Ifges(4,1) * t135 + t368 * t235 + (t38 - t39) * t231 + (t348 + t328 * t235 + (-Ifges(5,6) + Ifges(6,6)) * t231) * t136 + ((-t319 + t326) * t235 + (-t325 - t374) * t231) * qJD(4)) * t186 + (t17 * t65 + t18 * t3 + t19 * t2) * t350 + (t13 * t59 + t16 * t60 + t28 * t81) * t351 + t46 * t344 + t112 * t346 + (mrSges(4,1) * t342 - Ifges(7,5) * t107 + Ifges(7,6) * t106 + t156 * t349) * t136 - (mrSges(4,2) * t342 + mrSges(4,3) * t344 + t326 * t231 + t325 * t235) * t135 + 0.2e1 * t18 * t24 + 0.2e1 * t19 * t25 + t32 * t52 + t31 * t53 + 0.2e1 * t59 * t58 + 0.2e1 * t60 * t56 + 0.2e1 * t17 * t61 + 0.2e1 * t65 * t9 + 0.2e1 * t70 * t55 + 0.2e1 * t71 * t57 + 0.2e1 * t81 * t45 + 0.2e1 * t2 * t90 + 0.2e1 * t3 * t91 - t106 * t6 + t107 * t7 + 0.2e1 * t28 * t111 + 0.2e1 * t21 * t120 + 0.2e1 * t22 * t121 + 0.2e1 * t16 * t122 + 0.2e1 * t13 * t123; t240 * t214 + (m(4) * (t232 * t84 - t236 * t85) + (t135 * t236 - t136 * t232) * mrSges(4,3) + ((t186 * mrSges(4,3) + t112) * t232 + (-t184 * mrSges(4,3) + t231 * t296 + t235 * t297) * t236 + m(5) * (t298 * t71 - t300 * t70 - t304) + m(6) * (t298 * t59 + t300 * t60) + m(4) * (t156 * t236 - t304)) * qJD(3)) * pkin(2) + t239 + (Ifges(3,5) * t237 - Ifges(3,6) * t233 + (-mrSges(3,1) * t237 + mrSges(3,2) * t233) * pkin(7)) * qJD(2) + m(6) * (t160 * t81 + t173 * t28) + m(7) * (t109 * t3 + t110 * t2 + t146 * t65 + t161 * t17 + t18 * t51 + t19 * t50) + t50 * t90 + t51 * t91 + t109 * t24 + t110 * t25 + t146 * t61 + t160 * t111 + t161 * t9 + t173 * t45 + t360 * t215; 0.2e1 * t250 + (-t109 * t133 - t110 * t134 - t251 * t50 + t252 * t51) * t285 + (t215 * t279 + t294) * t352 + (t160 * t173 + t294) * t351 + (t109 * t51 + t110 * t50 + t146 * t161) * t350 + t242 + 0.2e1 * t353 + t146 * t345 + t161 * t347 + t173 * t343 + 0.2e1 * t160 * t200 + 0.2e1 * t215 * t189; t240 * pkin(9) + m(7) * (t153 * t3 + t155 * t2 + t159 * t65 + t17 * t174 + t18 * t83 + t19 * t82) + t239 + m(6) * (t167 * t81 - t28 * t283) + t82 * t90 + t83 * t91 + t153 * t24 + t155 * t25 + t159 * t61 + t167 * t111 + t174 * t9 - t283 * t45 - t360 * pkin(3); (-pkin(3) + t215) * t189 + (t173 - t283) * t188 + (t146 + t159) * t141 + t353 + t250 + t242 + (-(-t51 - t83) * t252 - (t50 + t82) * t251 - (t110 + t155) * t134 + (-t109 - t153) * t133) * mrSges(7,3) + (t167 + t160) * t200 + (t161 + t174) * t72 + m(5) * (-pkin(3) * t279 + t293) + m(6) * (-t160 * t283 + t167 * t173 + t293) + m(7) * (t109 * t83 + t110 * t82 + t146 * t174 + t153 * t51 + t155 * t50 + t159 * t161); (t153 * t83 + t155 * t82 + t159 * t174) * t350 + 0.2e1 * (-m(6) * t283 + t200) * t167 + (-t133 * t153 - t155 * t134 - t251 * t82 + t252 * t83) * t285 + t242 + t159 * t345 + t174 * t347 - 0.2e1 * pkin(3) * t189 - t283 * t343; -t370 * t186 * qJD(4) + m(7) * (t162 * t19 + t163 * t18 + t197 * t3 + t198 * t2) + m(6) * (-pkin(4) * t16 + qJ(5) * t13 + qJD(5) * t59) + t245 + Ifges(5,6) * t306 + t13 * mrSges(6,3) - t16 * mrSges(6,1) - t21 * mrSges(5,2) + t22 * mrSges(5,1) - pkin(4) * t56 + qJ(5) * t58 + qJD(5) * t123 + t162 * t90 + t163 * t91 + t197 * t24 + t198 * t25 - t356; (t214 * t241 + t244) * qJD(4) + t243 + t214 * t281 + m(7) * (t109 * t163 + t110 * t162 + t197 * t51 + t198 * t50) + (-m(6) * t253 - t259 - t260) * t278 - t357; m(7) * (t153 * t163 + t155 * t162 + t197 * t83 + t198 * t82) + t244 * qJD(4) + (qJD(4) * t241 + t281) * pkin(9) + t243 - t358; -0.2e1 * t310 + (t162 * t198 + t163 * t197) * t350 + 0.2e1 * t311 + 0.2e1 * (m(6) * qJ(5) + mrSges(6,3)) * qJD(5); t230 * t25 + t234 * t24 + (-t230 * t91 + t234 * t90) * qJD(6) + m(7) * (t2 * t230 + t234 * t3 + (-t18 * t230 + t19 * t234) * qJD(6)) + m(6) * t16 + t56; m(7) * (t230 * t50 + t234 * t51 + (-t109 * t230 + t110 * t234) * qJD(6)) + (t214 * t289 + t264) * m(6) + t246; m(7) * (t230 * t82 + t234 * t83 + (-t153 * t230 + t155 * t234) * qJD(6)) + m(6) * pkin(9) * t289 + t246; mrSges(7,1) * t287 + m(7) * (t162 * t230 + t163 * t234 + (-t197 * t230 + t198 * t234) * qJD(6)) + mrSges(7,2) * t286; 0; -t282 + t356; t295 + t357; t295 + t358; t310 - t311; (-mrSges(7,1) * t230 - mrSges(7,2) * t234) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
