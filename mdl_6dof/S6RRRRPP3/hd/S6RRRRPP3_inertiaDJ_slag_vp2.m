% Calculate time derivative of joint inertia matrix for
% S6RRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:53:56
% EndTime: 2019-03-09 20:54:09
% DurationCPUTime: 5.60s
% Computational Cost: add. (4661->465), mult. (10628->610), div. (0->0), fcn. (9041->6), ass. (0->202)
t357 = Ifges(5,1) + Ifges(7,3);
t356 = Ifges(7,2) + Ifges(6,3);
t209 = sin(qJ(3));
t210 = sin(qJ(2));
t212 = cos(qJ(3));
t213 = cos(qJ(2));
t155 = t209 * t210 - t212 * t213;
t346 = Ifges(5,5) + Ifges(7,5);
t355 = t155 * t346;
t350 = Ifges(7,4) + Ifges(6,5);
t354 = -Ifges(6,4) + t346;
t208 = sin(qJ(4));
t211 = cos(qJ(4));
t353 = t208 ^ 2 + t211 ^ 2;
t299 = Ifges(7,6) * t208;
t305 = Ifges(5,4) * t208;
t352 = t357 * t211 + t299 - t305;
t298 = Ifges(7,6) * t211;
t300 = Ifges(6,6) * t211;
t351 = t356 * t208 + t298 - t300;
t301 = Ifges(6,6) * t208;
t349 = -t211 * t356 + t299 - t301;
t304 = Ifges(5,4) * t211;
t348 = t208 * t357 - t298 + t304;
t156 = t209 * t213 + t212 * t210;
t329 = qJD(2) + qJD(3);
t109 = t329 * t156;
t275 = qJD(4) * t208;
t259 = t156 * t275;
t108 = t329 * t155;
t289 = t108 * t211;
t223 = t259 + t289;
t274 = qJD(4) * t211;
t290 = t108 * t208;
t224 = t156 * t274 - t290;
t345 = (-Ifges(5,4) + Ifges(7,6)) * t224 - t357 * t223 + t346 * t109;
t344 = t356 * t224 + (Ifges(6,6) - Ifges(7,6)) * t223 + t350 * t109;
t343 = t156 * t352 + t355;
t342 = t155 * t350 + t351 * t156;
t172 = -t211 * mrSges(5,1) + t208 * mrSges(5,2);
t341 = -mrSges(4,1) + t172;
t340 = t351 * qJD(4);
t339 = t352 * qJD(4);
t297 = pkin(2) * qJD(3);
t265 = t212 * t297;
t336 = t353 * t265;
t318 = -pkin(8) - pkin(7);
t181 = t318 * t210;
t183 = t318 * t213;
t129 = t181 * t209 - t183 * t212;
t266 = pkin(2) * qJD(2) * t210;
t47 = pkin(3) * t109 + pkin(9) * t108 + t266;
t260 = qJD(2) * t318;
t169 = t210 * t260;
t242 = t213 * t260;
t330 = t212 * t181 + t183 * t209;
t55 = qJD(3) * t330 + t212 * t169 + t209 * t242;
t191 = -pkin(2) * t213 - pkin(1);
t92 = t155 * pkin(3) - t156 * pkin(9) + t191;
t11 = -t129 * t275 + t208 * t47 + t211 * t55 + t92 * t274;
t12 = -t129 * t274 - t208 * t55 + t211 * t47 - t92 * t275;
t335 = t11 * t211 - t12 * t208;
t334 = Ifges(6,1) + Ifges(7,1) + Ifges(5,3);
t333 = t346 * t274 + t275 * t350;
t50 = t211 * t129 + t208 * t92;
t42 = -qJ(5) * t155 - t50;
t116 = t208 * t129;
t49 = t211 * t92 - t116;
t43 = -pkin(4) * t155 - t49;
t7 = -qJ(5) * t109 - qJD(5) * t155 - t11;
t9 = -pkin(4) * t109 - t12;
t332 = t208 * t9 - t211 * t7 + t43 * t274 + t42 * t275;
t30 = mrSges(5,1) * t224 - mrSges(5,2) * t223;
t56 = qJD(3) * t129 + t169 * t209 - t212 * t242;
t331 = m(5) * t56 + t30;
t328 = (-mrSges(4,2) + (mrSges(6,1) + mrSges(5,3)) * t353) * t265;
t327 = 0.2e1 * m(5);
t326 = 2 * m(6);
t325 = 2 * m(7);
t324 = -2 * Ifges(4,4);
t323 = 0.2e1 * t56;
t322 = -0.2e1 * t330;
t235 = -mrSges(7,2) * t211 + mrSges(7,3) * t208;
t158 = t235 * qJD(4);
t321 = 0.2e1 * t158;
t236 = -t208 * mrSges(6,2) - t211 * mrSges(6,3);
t159 = t236 * qJD(4);
t320 = 0.2e1 * t159;
t173 = -mrSges(7,2) * t208 - mrSges(7,3) * t211;
t319 = 0.2e1 * t173;
t312 = pkin(2) * t212;
t309 = -mrSges(6,1) - mrSges(7,1);
t207 = -pkin(4) - qJ(6);
t287 = t156 * t208;
t94 = -mrSges(5,2) * t155 - mrSges(5,3) * t287;
t97 = mrSges(6,1) * t287 - mrSges(6,3) * t155;
t307 = t94 - t97;
t286 = t156 * t211;
t95 = mrSges(5,1) * t155 - mrSges(5,3) * t286;
t99 = mrSges(6,1) * t286 + mrSges(6,2) * t155;
t306 = -t95 + t99;
t303 = Ifges(5,6) * t155;
t302 = Ifges(5,6) * t211;
t294 = t330 * t56;
t293 = -mrSges(6,1) * t289 + t109 * mrSges(6,2);
t292 = qJ(5) * t208;
t291 = qJ(5) * t211;
t288 = t330 * t209;
t284 = t208 * t212;
t282 = t211 * t212;
t189 = pkin(2) * t209 + pkin(9);
t281 = t336 * t189;
t280 = t336 * pkin(9);
t194 = mrSges(7,1) * t274;
t279 = mrSges(6,1) * t274 + t194;
t276 = qJ(5) * qJD(5);
t273 = qJD(5) * t211;
t272 = qJD(6) * t208;
t271 = 0.2e1 * mrSges(7,1);
t233 = -Ifges(5,2) * t208 + t304;
t63 = t156 * t233 + t303;
t270 = -t63 + t342;
t231 = -Ifges(6,2) * t211 + t301;
t68 = Ifges(6,4) * t155 + t156 * t231;
t269 = -t68 + t343;
t267 = m(6) * t273;
t264 = pkin(9) * t274;
t202 = t209 * t297;
t263 = mrSges(7,1) * t275;
t253 = -t275 / 0.2e1;
t251 = -t274 / 0.2e1;
t250 = t274 / 0.2e1;
t249 = -mrSges(7,1) * t289 - t109 * mrSges(7,3);
t248 = 0.2e1 * t266;
t164 = t231 * qJD(4);
t178 = Ifges(5,2) * t211 + t305;
t247 = -qJD(4) * t178 - t164;
t246 = pkin(4) * t275 - qJD(5) * t208;
t239 = pkin(4) * t287 - t330;
t177 = -Ifges(6,2) * t208 - t300;
t238 = -qJD(4) * t177 - t340;
t237 = t208 * mrSges(5,1) + t211 * mrSges(5,2);
t171 = t211 * mrSges(6,2) - t208 * mrSges(6,3);
t232 = -Ifges(6,4) * t211 - Ifges(5,6) * t208;
t227 = -pkin(4) * t211 - t292;
t226 = qJ(6) * t208 - t291;
t225 = t341 * t202;
t170 = -pkin(3) + t227;
t165 = t233 * qJD(4);
t222 = t211 * t165 + t339 * t208 + t348 * t274 + t349 * t275;
t145 = t207 * t211 - pkin(3) - t292;
t139 = -qJ(5) * t274 + t246;
t221 = t189 * t274 + t208 * t265;
t40 = -mrSges(7,1) * t224 + t109 * mrSges(7,2);
t220 = Ifges(6,4) * t259 + t334 * t109 + t350 * t224 - t346 * t289;
t219 = mrSges(6,1) * t273 + t207 * t194 + (-t272 + t273) * mrSges(7,1) + t333;
t124 = qJ(6) * t275 + (-qJ(5) * qJD(4) - qJD(6)) * t211 + t246;
t218 = (qJ(5) * t309 - Ifges(5,6)) * t208 + (-mrSges(6,1) * pkin(4) - Ifges(6,4)) * t211;
t217 = pkin(4) * t224 + qJ(5) * t259 + t56;
t216 = m(6) * t227 + t171 + t172;
t36 = mrSges(5,1) * t109 + mrSges(5,3) * t223;
t37 = -mrSges(5,2) * t109 - mrSges(5,3) * t224;
t39 = mrSges(6,1) * t224 - mrSges(6,3) * t109;
t41 = -mrSges(6,1) * t259 + t293;
t215 = (t37 - t39) * t211 + (-t36 + t41) * t208 + (-t208 * t307 + t211 * t306) * qJD(4) + m(6) * t332 + m(5) * (-t274 * t49 - t275 * t50 + t335);
t14 = (qJ(5) * t108 - qJD(5) * t156) * t211 + t217;
t160 = t237 * qJD(4);
t2 = -pkin(5) * t223 - qJD(6) * t155 + t109 * t207 - t12;
t21 = -Ifges(5,4) * t223 - Ifges(5,2) * t224 + Ifges(5,6) * t109;
t26 = Ifges(6,4) * t109 + Ifges(6,2) * t223 + Ifges(6,6) * t224;
t28 = t116 + (pkin(5) * t156 - t92) * t211 + t207 * t155;
t32 = -pkin(5) * t287 - t42;
t4 = -pkin(5) * t224 - t7;
t44 = t156 * t226 + t239;
t54 = -qJ(5) * t286 + t239;
t6 = -t226 * t108 + (t272 + (qJ(6) * qJD(4) - qJD(5)) * t211) * t156 + t217;
t214 = (t156 * t251 + t290 / 0.2e1) * t178 + t28 * t194 + t68 * t251 + t63 * t253 + t348 * (t156 * t253 - t289 / 0.2e1) + t349 * (t156 * t250 - t290 / 0.2e1) + (t156 * t177 + t342) * t275 / 0.2e1 + t343 * t250 - t344 * t211 / 0.2e1 + t345 * t208 / 0.2e1 + t341 * t56 + ((-t208 * t50 - t211 * t49) * qJD(4) + t335) * mrSges(5,3) + (-Ifges(6,4) * t274 - Ifges(5,6) * t275 + t333) * t155 / 0.2e1 + t332 * mrSges(6,1) + t211 * t21 / 0.2e1 - t208 * t26 / 0.2e1 + t14 * t171 + t6 * t173 + t44 * t158 + t54 * t159 - Ifges(4,6) * t109 - Ifges(4,5) * t108 - t55 * mrSges(4,2) - t330 * t160 + (t2 * t208 + t211 * t4) * mrSges(7,1) + (t340 / 0.2e1 - t165 / 0.2e1) * t287 + (t339 / 0.2e1 - t164 / 0.2e1) * t286 + (t208 * t354 - t211 * t350 + t302) * t109 / 0.2e1 - t32 * t263 + t177 * t289 / 0.2e1;
t204 = t211 * pkin(5);
t203 = t208 * pkin(5);
t201 = pkin(5) * t274;
t190 = -pkin(3) - t312;
t182 = pkin(9) * t211 + t204;
t180 = pkin(9) * t208 + t203;
t168 = t201 + t264;
t167 = (-pkin(5) - pkin(9)) * t275;
t152 = t189 * t211 + t204;
t151 = t189 * t208 + t203;
t146 = t170 - t312;
t133 = t145 - t312;
t132 = t139 + t202;
t131 = t201 + t221;
t130 = t211 * t265 + (-pkin(5) - t189) * t275;
t98 = -mrSges(7,1) * t287 + mrSges(7,2) * t155;
t96 = mrSges(7,1) * t286 - mrSges(7,3) * t155;
t93 = t202 + t124;
t84 = t237 * t156;
t83 = t236 * t156;
t82 = t235 * t156;
t38 = -mrSges(7,1) * t259 + t249;
t31 = mrSges(7,2) * t223 + mrSges(7,3) * t224;
t29 = -mrSges(6,2) * t224 + mrSges(6,3) * t223;
t1 = [(mrSges(4,2) * t248 + mrSges(4,3) * t323 - 0.2e1 * Ifges(4,1) * t108 + (-t26 + t345) * t211 + (-t21 + t344) * t208 + (t324 + t354 * t211 + (-Ifges(5,6) + t350) * t208) * t109 + ((t270 - t303) * t211 + (-t269 - t355) * t208) * qJD(4)) * t156 + (t11 * t50 + t12 * t49 - t294) * t327 + 0.2e1 * m(4) * (t129 * t55 + t191 * t266 - t294) + 0.2e1 * (-pkin(1) * (mrSges(3,1) * t210 + mrSges(3,2) * t213) + (-Ifges(3,2) + Ifges(3,1)) * t210 * t213 + (-t210 ^ 2 + t213 ^ 2) * Ifges(3,4)) * qJD(2) + (t220 + mrSges(4,1) * t248 + ((2 * Ifges(4,2)) + t334) * t109 - (t324 + t232) * t108 - 0.2e1 * t55 * mrSges(4,3)) * t155 + 0.2e1 * t6 * t82 + 0.2e1 * t14 * t83 + 0.2e1 * t11 * t94 + 0.2e1 * t12 * t95 + 0.2e1 * t2 * t96 + 0.2e1 * t7 * t97 + 0.2e1 * t4 * t98 + 0.2e1 * t9 * t99 + 0.2e1 * t54 * t29 + 0.2e1 * t32 * t40 + 0.2e1 * t42 * t39 + 0.2e1 * t43 * t41 + 0.2e1 * t44 * t31 + 0.2e1 * t49 * t36 + 0.2e1 * t50 * t37 + 0.2e1 * t28 * t38 + 0.2e1 * (mrSges(4,1) * t191 - mrSges(4,3) * t129) * t109 + t30 * t322 + t84 * t323 + (t2 * t28 + t32 * t4 + t44 * t6) * t325 + (t14 * t54 + t42 * t7 + t43 * t9) * t326 - (0.2e1 * t191 * mrSges(4,2) + mrSges(4,3) * t322 + t270 * t208 + t269 * t211) * t108; t214 + t215 * t189 + (m(4) * (t209 * t55 - t212 * t56) + (t108 * t212 - t109 * t209) * mrSges(4,3) + ((t156 * mrSges(4,3) + t84) * t209 + (-t155 * mrSges(4,3) + t208 * t306 + t211 * t307) * t212 + m(4) * (t129 * t212 - t288) + m(6) * (-t282 * t42 + t284 * t43) + m(5) * (t282 * t50 - t284 * t49 - t288)) * qJD(3)) * pkin(2) + (Ifges(3,5) * t213 - Ifges(3,6) * t210 + (-mrSges(3,1) * t213 + mrSges(3,2) * t210) * pkin(7)) * qJD(2) + t151 * t38 + t152 * t40 + m(6) * (t132 * t54 + t14 * t146) + t146 * t29 + t130 * t98 + t131 * t96 + t132 * t83 + t133 * t31 + t93 * t82 + m(7) * (t130 * t32 + t131 * t28 + t133 * t6 + t151 * t2 + t152 * t4 + t44 * t93) + t331 * t190; t222 + ((-qJD(4) * t152 + t131) * t271 + t247) * t208 + ((qJD(4) * t151 + t130) * t271 + t238) * t211 + 0.2e1 * t190 * t160 + 0.2e1 * t132 * t171 + t93 * t319 + t133 * t321 + t146 * t320 + 0.2e1 * t225 + 0.2e1 * t328 + (t190 * t202 + t281) * t327 + (t130 * t152 + t131 * t151 + t133 * t93) * t325 + (t132 * t146 + t281) * t326; m(6) * (t139 * t54 + t14 * t170) + t214 + m(7) * (t124 * t44 + t145 * t6 + t167 * t32 + t168 * t28 + t180 * t2 + t182 * t4) + t215 * pkin(9) + t170 * t29 + t180 * t38 + t182 * t40 + t167 * t98 + t168 * t96 + t139 * t83 + t145 * t31 + t124 * t82 - t331 * pkin(3); t238 * t211 + t247 * t208 + t222 + t328 + t225 + (t124 + t93) * t173 + (t139 + t132) * t171 + m(5) * (-pkin(3) * t202 + t280) + m(7) * (t124 * t133 + t130 * t182 + t131 * t180 + t145 * t93 + t151 * t168 + t152 * t167) + m(6) * (t132 * t170 + t139 * t146 + t280) + (t190 - pkin(3)) * t160 + (t170 + t146) * t159 + (t145 + t133) * t158 + ((t130 + t167) * t211 + (t131 + t168) * t208 + ((t151 + t180) * t211 + (-t152 - t182) * t208) * qJD(4)) * mrSges(7,1); (t124 * t145 + t167 * t182 + t168 * t180) * t325 + t124 * t319 + t145 * t321 + t170 * t320 - 0.2e1 * pkin(3) * t160 + 0.2e1 * (m(6) * t170 + t171) * t139 + ((-qJD(4) * t182 + t168) * t271 + t247) * t208 + ((qJD(4) * t180 + t167) * t271 + t238) * t211 + t222; t220 + t207 * t38 - qJD(6) * t96 - pkin(4) * t41 - t11 * mrSges(5,2) + t12 * mrSges(5,1) - t7 * mrSges(6,3) + t9 * mrSges(6,2) - t2 * mrSges(7,3) + t4 * mrSges(7,2) + (-t346 * t208 - t302) * t156 * qJD(4) + (-t97 + t98) * qJD(5) + (-t39 + t40) * qJ(5) + m(7) * (qJ(5) * t4 + qJD(5) * t32 - qJD(6) * t28 + t2 * t207) + m(6) * (-pkin(4) * t9 - qJ(5) * t7 - qJD(5) * t42) - t232 * t108; t130 * mrSges(7,2) - t131 * mrSges(7,3) + t189 * t267 + m(7) * (qJ(5) * t130 + qJD(5) * t152 - qJD(6) * t151 + t131 * t207) + (m(6) * (-pkin(4) * t208 + t291) - t236 - t237) * t265 + (t189 * t216 + t218) * qJD(4) + t219; t167 * mrSges(7,2) - t168 * mrSges(7,3) + pkin(9) * t267 + m(7) * (qJ(5) * t167 + qJD(5) * t182 - qJD(6) * t180 + t168 * t207) + (pkin(9) * t216 + t218) * qJD(4) + t219; 0.2e1 * m(6) * t276 + 0.2e1 * qJD(6) * mrSges(7,3) + 0.2e1 * m(7) * (-qJD(6) * t207 + t276) + 0.2e1 * (mrSges(6,3) + mrSges(7,2)) * qJD(5); m(6) * t9 + m(7) * t2 + t259 * t309 + t249 + t293; m(6) * t221 + m(7) * t131 + t279; m(6) * t264 + m(7) * t168 + t279; -m(7) * qJD(6); 0; m(7) * t4 + t40; m(7) * t130 - t263; m(7) * t167 - t263; m(7) * qJD(5); 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
