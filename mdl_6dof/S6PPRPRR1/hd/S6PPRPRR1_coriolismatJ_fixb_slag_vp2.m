% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PPRPRR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_coriolismatJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRPRR1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRPRR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRPRR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:42:19
% EndTime: 2019-03-08 18:42:24
% DurationCPUTime: 3.78s
% Computational Cost: add. (10193->335), mult. (29522->530), div. (0->0), fcn. (35032->14), ass. (0->220)
t216 = sin(qJ(5));
t219 = cos(qJ(5));
t217 = sin(qJ(3));
t312 = sin(pkin(6));
t275 = sin(pkin(12)) * t312;
t339 = cos(qJ(3));
t266 = cos(pkin(12)) * t312;
t311 = sin(pkin(7));
t314 = cos(pkin(7));
t315 = cos(pkin(6));
t370 = t314 * t266 + t315 * t311;
t131 = -t217 * t275 + t339 * t370;
t132 = t217 * t370 + t339 * t275;
t213 = sin(pkin(13));
t313 = cos(pkin(13));
t240 = t213 * t131 + t313 * t132;
t224 = -t311 * t266 + t315 * t314;
t83 = t216 * t240 - t224 * t219;
t84 = t224 * t216 + t219 * t240;
t378 = -t216 * t83 - t219 * t84 + t240;
t362 = -m(7) / 0.2e1;
t265 = t313 * t311;
t268 = t339 * t311;
t164 = t213 * t268 + t217 * t265;
t139 = t219 * t164 + t216 * t314;
t138 = t164 * t216 - t219 * t314;
t294 = t216 * t138;
t377 = -t139 * t219 + t164 - t294;
t215 = sin(qJ(6));
t209 = t215 ^ 2;
t218 = cos(qJ(6));
t211 = t218 ^ 2;
t285 = t209 + t211;
t376 = mrSges(6,2) - t285 * (m(7) * pkin(10) + mrSges(7,3));
t287 = t218 * t219;
t96 = -t313 * t131 + t132 * t213;
t67 = t215 * t240 - t287 * t96;
t316 = t67 * t218;
t295 = t215 * t219;
t66 = t218 * t240 + t295 * t96;
t317 = t66 * t215;
t257 = t316 - t317;
t323 = t216 * t96;
t375 = (pkin(5) * t323 + t257 * pkin(10)) * t362 - (t316 / 0.2e1 - t317 / 0.2e1) * mrSges(7,3);
t359 = mrSges(6,1) / 0.2e1;
t190 = -mrSges(7,1) * t218 + mrSges(7,2) * t215;
t292 = t216 * t190;
t374 = -t292 / 0.2e1;
t320 = t219 * t96;
t33 = (t257 + t320) * t216;
t274 = t217 * t311;
t163 = t213 * t274 - t339 * t265;
t115 = -t163 * t287 + t215 * t164;
t307 = t115 * t218;
t114 = t163 * t295 + t218 * t164;
t308 = t114 * t215;
t255 = t307 - t308;
t69 = (t163 * t219 + t255) * t216;
t373 = (t33 * qJD(1) / 0.4e1 + t69 * qJD(2) / 0.4e1) * m(7);
t372 = mrSges(7,3) * (t211 / 0.2e1 + t209 / 0.2e1);
t328 = Ifges(7,6) * t215;
t330 = Ifges(7,5) * t218;
t249 = t330 / 0.2e1 - t328 / 0.2e1;
t371 = Ifges(6,4) - t249;
t269 = (-0.1e1 + t285) * t216;
t296 = t215 * t216;
t283 = mrSges(7,3) * t296;
t186 = mrSges(7,2) * t219 - t283;
t291 = t216 * t218;
t188 = -mrSges(7,1) * t219 - mrSges(7,3) * t291;
t342 = -t218 / 0.2e1;
t345 = -t215 / 0.2e1;
t369 = t186 * t345 + t188 * t342;
t337 = t216 * pkin(5);
t196 = -pkin(10) * t219 + t337;
t205 = pkin(3) * t213 + pkin(9);
t143 = t196 * t218 + t205 * t296;
t144 = t196 * t215 - t205 * t291;
t368 = -t143 * t215 + t144 * t218;
t366 = -m(7) * pkin(5) - mrSges(6,1) + t190;
t206 = -t313 * pkin(3) - pkin(4);
t360 = m(5) * pkin(3);
t364 = qJD(3) * (m(6) * t206 - t219 * mrSges(6,1) + t216 * mrSges(6,2) - t313 * t360 - mrSges(5,1));
t322 = t218 * mrSges(7,2);
t324 = t215 * mrSges(7,1);
t264 = t322 + t324;
t177 = t264 * t216;
t210 = t216 ^ 2;
t212 = t219 ^ 2;
t363 = m(6) * t212 * t205 + t216 * t177 + t213 * t360 - mrSges(5,2) + (t210 + t212) * mrSges(6,3);
t361 = m(7) / 0.2e1;
t358 = mrSges(7,1) / 0.2e1;
t357 = -mrSges(7,2) / 0.2e1;
t58 = -t215 * t84 + t218 * t96;
t356 = t58 / 0.2e1;
t355 = t83 / 0.2e1;
t109 = t139 * t218 + t163 * t215;
t354 = -t109 / 0.2e1;
t353 = -t138 / 0.2e1;
t352 = t138 / 0.2e1;
t351 = t177 / 0.2e1;
t178 = t219 * t264;
t350 = -t178 / 0.2e1;
t349 = t186 / 0.2e1;
t348 = t188 / 0.2e1;
t321 = t219 * mrSges(6,2);
t192 = t216 * mrSges(6,1) + t321;
t347 = t192 / 0.2e1;
t346 = -t205 / 0.2e1;
t344 = t215 / 0.2e1;
t343 = t216 / 0.2e1;
t341 = t218 / 0.2e1;
t340 = t219 / 0.2e1;
t176 = -mrSges(7,1) * t291 + mrSges(7,2) * t296;
t338 = pkin(5) * t176;
t289 = t218 * t186;
t298 = t215 * t188;
t235 = t350 - t298 / 0.2e1 + t289 / 0.2e1;
t243 = m(7) * t368;
t179 = -t219 * pkin(5) - t216 * pkin(10) + t206;
t134 = t218 * t179 - t205 * t295;
t135 = t215 * t179 + t205 * t287;
t254 = t134 * t215 - t135 * t218;
t187 = -t216 * mrSges(7,2) - mrSges(7,3) * t295;
t288 = t218 * t187;
t189 = t216 * mrSges(7,1) - mrSges(7,3) * t287;
t297 = t215 * t189;
t45 = (t210 - t212) * t205 * t361 + (t254 * t362 + t235) * t219 + (t243 / 0.2e1 + t288 / 0.2e1 - t297 / 0.2e1 + t351) * t216;
t68 = t176 * t340 - t210 * t372 + t216 * t369;
t335 = t45 * qJD(5) + t68 * qJD(6);
t334 = m(7) * qJD(4);
t333 = mrSges(7,3) * t216;
t332 = Ifges(7,4) * t215;
t331 = Ifges(7,4) * t218;
t329 = Ifges(7,5) * t219;
t327 = Ifges(7,6) * t219;
t319 = t58 * t215;
t59 = t215 * t96 + t218 * t84;
t318 = t59 * t218;
t108 = -t139 * t215 + t163 * t218;
t310 = t108 * t215;
t309 = t109 * t218;
t304 = t163 * t216;
t303 = t205 * t177;
t302 = t205 * t216;
t301 = t210 * t205;
t261 = -Ifges(7,2) * t215 + t331;
t241 = t261 * t216;
t165 = t241 - t327;
t299 = t215 * t165;
t263 = Ifges(7,1) * t218 - t332;
t242 = t263 * t216;
t167 = t242 - t329;
t290 = t218 * t167;
t286 = t68 * qJD(3);
t282 = qJD(3) * t361;
t281 = -t334 / 0.2e1;
t280 = -t333 / 0.2e1;
t279 = t321 / 0.2e1;
t271 = t285 * t219;
t262 = Ifges(7,1) * t215 + t331;
t260 = Ifges(7,2) * t218 + t332;
t259 = -t328 + t330;
t258 = Ifges(7,5) * t215 + Ifges(7,6) * t218;
t103 = (mrSges(7,2) * pkin(5) - t331) * t218 + (pkin(5) * mrSges(7,1) + t332 + (-Ifges(7,1) + Ifges(7,2)) * t218) * t215;
t247 = t143 * t358 + t144 * t357;
t30 = -t338 / 0.2e1 + (Ifges(7,3) / 0.2e1 + pkin(10) * t372) * t216 + (0.3e1 / 0.4e1 * t329 + pkin(10) * t348 - t167 / 0.4e1 + (mrSges(7,2) * t346 + (Ifges(7,2) / 0.2e1 - Ifges(7,1) / 0.4e1) * t218) * t216) * t218 + (-0.3e1 / 0.4e1 * t327 + pkin(10) * t349 + t165 / 0.4e1 + (0.3e1 / 0.2e1 * t331 + mrSges(7,1) * t346 + (Ifges(7,1) / 0.2e1 - Ifges(7,2) / 0.4e1) * t215) * t216) * t215 + t247;
t252 = t30 * qJD(3) + t103 * qJD(5);
t251 = t67 * t357 + t66 * t358;
t250 = -t318 + t84 + t319;
t248 = -t114 * mrSges(7,1) / 0.2e1 + t115 * mrSges(7,2) / 0.2e1;
t246 = t322 / 0.2e1 + t324 / 0.2e1;
t245 = t139 - t309 + t310;
t12 = m(7) * t250 * t83;
t17 = (-t250 * t219 - t83 * t269) * t361;
t9 = (t250 * t138 + t245 * t83) * t361;
t238 = t12 * qJD(1) + t9 * qJD(2) + t17 * qJD(4);
t37 = m(7) * t245 * t138;
t43 = (-t138 * t269 - t245 * t219) * t361;
t237 = t9 * qJD(1) + t37 * qJD(2) + t43 * qJD(4);
t236 = t205 * t219 + t254;
t234 = t279 + (-t190 / 0.2e1 + t359) * t216;
t4 = (t108 * t66 + t109 * t67 + t114 * t58 + t115 * t59 + (-t138 * t96 - t163 * t83) * t216) * t361 + m(6) * (t378 * t163 + t377 * t96) / 0.2e1;
t7 = m(7) * (-t83 * t323 + t58 * t66 + t59 * t67) + m(6) * t378 * t96;
t233 = -t7 * qJD(1) - t4 * qJD(2) + t33 * t281;
t232 = t143 * t108 + t144 * t109 + t139 * t302;
t29 = m(6) * t377 * t163 + (t108 * t114 + t109 * t115 - t294 * t163) * m(7);
t231 = -t4 * qJD(1) - t29 * qJD(2) + t69 * t281;
t230 = -t58 * t186 / 0.2e1 + t59 * t348 + t176 * t355;
t229 = -t139 * t177 / 0.2e1 - t108 * t189 / 0.2e1 + t187 * t354;
t221 = (-t308 / 0.2e1 + t307 / 0.2e1) * mrSges(7,3) + (pkin(5) * t304 + t255 * pkin(10)) * t361;
t223 = t236 * t361 - t235;
t10 = (-t192 / 0.2e1 + t234) * t163 + t232 * t362 - t223 * t138 + t221 + t229;
t166 = Ifges(7,6) * t216 + t261 * t219;
t168 = Ifges(7,5) * t216 + t263 * t219;
t27 = m(7) * (t134 * t143 + t135 * t144) + t144 * t186 + t135 * t187 + t143 * t188 + t134 * t189 + t206 * t192 + (t290 / 0.2e1 - t299 / 0.2e1 + t303 + t371 * t219) * t219 + (t168 * t341 + t166 * t345 + t205 * t178 - t371 * t216 + (m(7) * t205 ^ 2 + Ifges(6,1) - Ifges(6,2) - Ifges(7,3)) * t219) * t216;
t220 = t223 * t83 + (t143 * t58 + t144 * t59 + t84 * t302) * t361 + t189 * t356 + t59 * t187 / 0.2e1 + t84 * t351 + t96 * t347;
t3 = -t234 * t96 + t220 + t375;
t228 = t3 * qJD(1) - t10 * qJD(2) + t27 * qJD(3) + t45 * qJD(4);
t222 = (t310 / 0.2e1 - t309 / 0.2e1) * t333 + t108 * t349 + t188 * t354 + t176 * t353;
t22 = t222 + t248;
t32 = t134 * t186 - t135 * t188 + (t258 * t340 + t167 * t345 + t165 * t342 - t205 * t176 + (t260 * t344 + t262 * t342) * t216 + t254 * mrSges(7,3)) * t216;
t5 = (t318 / 0.2e1 - t319 / 0.2e1) * t333 + t230 + t251;
t227 = -t5 * qJD(1) + t22 * qJD(2) + t32 * qJD(3) + t68 * qJD(4);
t147 = t219 * t269;
t225 = t17 * qJD(1) + t43 * qJD(2) + t45 * qJD(3) + t147 * t334;
t140 = t163 * t301;
t137 = -t246 * t219 + t350;
t89 = t96 * t301;
t41 = t246 * t138 + t264 * t352;
t31 = t290 / 0.4e1 - t262 * t296 / 0.2e1 - t299 / 0.4e1 - t260 * t291 / 0.2e1 + t338 / 0.2e1 - t215 * t241 / 0.4e1 + t218 * t242 / 0.4e1 + t303 / 0.2e1 + Ifges(7,3) * t343 + t247 + (-t259 / 0.4e1 + t249) * t219 + (t280 * t285 + t369) * pkin(10);
t26 = t43 * qJD(5) + t69 * t282;
t21 = t222 - t248;
t15 = t246 * t83 + t264 * t355;
t11 = (t236 * t138 + t232) * t361 + t289 * t353 + t304 * t359 + t221 - t229 + (t178 + t298) * t352 + (t347 + t279 + t374) * t163;
t8 = t17 * qJD(5) + t33 * t282;
t6 = t280 * t318 + t283 * t356 - t230 + t251;
t2 = t96 * t374 + mrSges(6,2) * t320 / 0.2e1 + t323 * t359 + t220 - t375;
t1 = t4 * qJD(3) + t9 * qJD(5);
t13 = [t7 * qJD(3) + t12 * qJD(5), t1, t2 * qJD(5) + t6 * qJD(6) + t240 * t364 - t233 + (t67 * t186 + t66 * t188 + m(7) * (t134 * t66 + t135 * t67 - t89) - m(6) * t89 - t131 * mrSges(4,2) - t132 * mrSges(4,1) - t363 * t96) * qJD(3), t8, t2 * qJD(3) + (t366 * t84 + t376 * t83) * qJD(5) + t15 * qJD(6) + t238, t6 * qJD(3) + t15 * qJD(5) + (-mrSges(7,1) * t59 - mrSges(7,2) * t58) * qJD(6); t1, t29 * qJD(3) + t37 * qJD(5), t11 * qJD(5) + t21 * qJD(6) + t164 * t364 - t231 + (-mrSges(4,2) * t268 - mrSges(4,1) * t274 + m(7) * (t134 * t114 + t135 * t115 - t140) + t114 * t188 + t115 * t186 - m(6) * t140 - t363 * t163) * qJD(3), t26, t11 * qJD(3) + (t376 * t138 + t366 * t139) * qJD(5) + t41 * qJD(6) + t237, t21 * qJD(3) + t41 * qJD(5) + (-mrSges(7,1) * t109 - mrSges(7,2) * t108) * qJD(6); t3 * qJD(5) - t5 * qJD(6) + t233, -t10 * qJD(5) + t22 * qJD(6) + t231, qJD(5) * t27 + qJD(6) * t32, t335 - 0.2e1 * t373, t31 * qJD(6) + t228 + (mrSges(6,2) * t302 - Ifges(6,6) * t216 - pkin(5) * t178 + t166 * t341 + t168 * t344 + t258 * t343 + (t243 + t288 - t297) * pkin(10) + (t205 * t366 + t260 * t345 + t262 * t341 + Ifges(6,5)) * t219 + t368 * mrSges(7,3)) * qJD(5), t31 * qJD(5) + (-mrSges(7,1) * t135 - mrSges(7,2) * t134 - t216 * t258) * qJD(6) + t227; t8, t26, t335 + 0.2e1 * t373, m(7) * t147 * qJD(5) (t292 + m(7) * (pkin(10) * t271 - t337) + mrSges(7,3) * t271 - t192) * qJD(5) + t137 * qJD(6) + t225, t137 * qJD(5) + t176 * qJD(6) + t286; -qJD(3) * t3 - t238, qJD(3) * t10 - t237, -qJD(6) * t30 - t228, -t225, -t103 * qJD(6) (pkin(10) * t190 + t259) * qJD(6) - t252; t5 * qJD(3), -t22 * qJD(3), qJD(5) * t30 - t227, -t286, t252, 0;];
Cq  = t13;
