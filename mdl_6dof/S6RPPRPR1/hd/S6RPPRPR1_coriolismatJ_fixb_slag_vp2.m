% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
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
% Datum: 2019-03-09 01:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPRPR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:39:14
% EndTime: 2019-03-09 01:39:21
% DurationCPUTime: 4.42s
% Computational Cost: add. (13976->368), mult. (27591->521), div. (0->0), fcn. (31170->10), ass. (0->207)
t213 = sin(pkin(10));
t215 = cos(pkin(10));
t218 = sin(qJ(4));
t326 = cos(qJ(4));
t201 = t213 * t326 + t218 * t215;
t212 = sin(pkin(11));
t214 = cos(pkin(11));
t217 = sin(qJ(6));
t219 = cos(qJ(6));
t245 = t217 * t212 - t214 * t219;
t152 = t245 * t201;
t200 = t212 * t219 + t217 * t214;
t231 = t200 * t201;
t363 = m(7) * (t152 * t245 + t200 * t231);
t365 = -t363 / 0.2e1;
t364 = t363 / 0.2e1;
t334 = -t152 / 0.2e1;
t198 = t213 * t218 - t215 * t326;
t153 = t200 * t198;
t156 = t245 * t198;
t63 = -t153 * t245 + t156 * t200;
t362 = m(7) * t63 * qJD(3);
t341 = -mrSges(7,3) / 0.2e1;
t361 = t152 * t200 * t341;
t329 = t200 / 0.2e1;
t360 = t152 * t329;
t323 = pkin(8) + qJ(5);
t202 = t323 * t212;
t204 = t323 * t214;
t177 = -t202 * t219 - t217 * t204;
t178 = -t217 * t202 + t204 * t219;
t286 = t201 * t212;
t168 = -mrSges(6,2) * t198 - mrSges(6,3) * t286;
t282 = t214 * t168;
t285 = t201 * t214;
t169 = mrSges(6,1) * t198 - mrSges(6,3) * t285;
t284 = t212 * t169;
t237 = -t284 / 0.2e1 + t282 / 0.2e1;
t108 = t198 * mrSges(7,1) + t152 * mrSges(7,3);
t107 = -t198 * mrSges(7,2) - mrSges(7,3) * t231;
t357 = t107 / 0.2e1;
t238 = t108 * t329 + t245 * t357;
t243 = -cos(pkin(9)) * pkin(1) - pkin(3) * t215 - pkin(2);
t151 = pkin(4) * t198 - qJ(5) * t201 + t243;
t205 = sin(pkin(9)) * pkin(1) + qJ(3);
t324 = pkin(7) + t205;
t190 = t324 * t215;
t258 = t324 * t213;
t165 = t190 * t326 - t218 * t258;
t86 = t214 * t151 - t165 * t212;
t87 = t212 * t151 + t214 * t165;
t248 = t212 * t86 - t214 * t87;
t56 = pkin(5) * t198 - pkin(8) * t285 + t86;
t64 = -pkin(8) * t286 + t87;
t31 = -t217 * t64 + t219 * t56;
t32 = t217 * t56 + t219 * t64;
t343 = -m(7) / 0.2e1;
t345 = m(6) / 0.2e1;
t359 = t248 * t345 + (t152 * t177 - t178 * t231 - t200 * t31 - t245 * t32) * t343 - t237 + t238;
t351 = t231 * t245;
t356 = t351 / 0.2e1;
t208 = t212 ^ 2;
t210 = t214 ^ 2;
t269 = -t210 / 0.2e1;
t318 = mrSges(6,3) * t198;
t355 = (-t208 / 0.2e1 + t269) * t318;
t337 = -t231 / 0.2e1;
t303 = t214 * mrSges(6,1);
t307 = t212 * mrSges(6,2);
t203 = -t303 + t307;
t353 = t203 - mrSges(5,1);
t273 = t208 + t210;
t164 = t190 * t218 + t326 * t258;
t297 = qJ(5) * t198;
t325 = pkin(4) * t201;
t171 = t297 + t325;
t92 = t164 * t212 + t214 * t171;
t93 = -t214 * t164 + t212 * t171;
t247 = -t212 * t92 + t214 * t93;
t311 = t156 * mrSges(7,2);
t313 = t153 * mrSges(7,1);
t277 = t313 / 0.2e1 - t311 / 0.2e1;
t283 = t212 * t198;
t118 = -pkin(5) * t283 + t165;
t319 = mrSges(6,2) * t214;
t320 = mrSges(6,1) * t212;
t350 = t118 * t343 - (-t319 / 0.2e1 - t320 / 0.2e1) * t198 + t277;
t191 = t198 * mrSges(5,2);
t290 = t156 * t245;
t294 = t153 * t200;
t349 = t191 + (-t290 - t294) * mrSges(7,3);
t348 = 0.2e1 * t201;
t347 = 2 * qJD(4);
t346 = -m(6) / 0.2e1;
t344 = m(6) / 0.4e1;
t342 = m(7) / 0.2e1;
t338 = t153 / 0.2e1;
t336 = t156 / 0.2e1;
t335 = t152 / 0.2e1;
t195 = Ifges(7,4) * t245;
t176 = Ifges(7,1) * t200 - t195;
t333 = t176 / 0.2e1;
t332 = -t178 / 0.2e1;
t331 = -t245 / 0.2e1;
t330 = t198 / 0.2e1;
t328 = t212 / 0.2e1;
t327 = t214 / 0.2e1;
t265 = t63 * t342;
t322 = qJD(4) * t265;
t317 = Ifges(6,4) * t212;
t316 = Ifges(6,4) * t214;
t315 = Ifges(7,4) * t152;
t314 = Ifges(7,4) * t200;
t309 = t245 * mrSges(7,1);
t308 = t200 * mrSges(7,2);
t306 = t212 * Ifges(6,2);
t305 = t212 * Ifges(6,6);
t302 = t214 * Ifges(6,5);
t88 = -t152 * mrSges(7,1) - t231 * mrSges(7,2);
t8 = t107 * t337 + t108 * t335 + t88 * t330 + (-t152 ^ 2 / 0.2e1 - t231 ^ 2 / 0.2e1) * mrSges(7,3);
t300 = t8 * qJD(1);
t179 = t198 * t201;
t229 = m(7) * (-t152 * t156 - t153 * t231 + t179);
t257 = t273 * t201;
t230 = m(6) * (-t198 * t257 + t179);
t27 = t229 / 0.2e1 + t230 / 0.2e1;
t296 = qJD(1) * t27;
t295 = t153 * t108;
t291 = t156 * t107;
t289 = t164 * t201;
t281 = t214 * t198;
t268 = m(7) / 0.4e1 + t344;
t42 = t364 + (t273 * t344 + t268) * t348;
t278 = t42 * qJD(1);
t275 = -Ifges(7,5) * t231 + Ifges(7,6) * t152;
t274 = -Ifges(7,5) * t245 - Ifges(7,6) * t200;
t271 = qJD(4) * t201;
t170 = t200 * mrSges(7,1) - mrSges(7,2) * t245;
t270 = t170 * qJD(6);
t262 = t170 * t330;
t172 = t308 + t309;
t261 = t172 / 0.2e1 + t203 / 0.2e1;
t259 = m(6) * t273;
t256 = t268 * t348;
t255 = t319 + t320;
t254 = t311 - t313;
t253 = mrSges(7,1) * t231 - mrSges(7,2) * t152;
t117 = pkin(5) * t286 + t164;
t127 = t201 * Ifges(6,6) + (t306 - t316) * t198;
t128 = t201 * Ifges(6,5) + (-t214 * Ifges(6,1) + t317) * t198;
t240 = Ifges(7,5) * t336 + Ifges(7,6) * t338;
t249 = t87 * t212 + t86 * t214;
t62 = pkin(5) * t201 + pkin(8) * t281 + t92;
t80 = pkin(8) * t283 + t93;
t38 = -t217 * t80 + t219 * t62;
t39 = t217 * t62 + t219 * t80;
t76 = Ifges(7,4) * t156 + Ifges(7,2) * t153 + t201 * Ifges(7,6);
t77 = -Ifges(7,2) * t231 + t198 * Ifges(7,6) - t315;
t78 = Ifges(7,1) * t156 + Ifges(7,4) * t153 + t201 * Ifges(7,5);
t150 = Ifges(7,4) * t231;
t79 = -Ifges(7,1) * t152 + t198 * Ifges(7,5) - t150;
t1 = -t243 * t191 + t93 * t168 + t92 * t169 + t118 * t253 + t78 * t334 + t77 * t338 + t76 * t337 + t117 * t254 + t79 * t336 + t38 * t108 + t39 * t107 + (t32 * t153 - t31 * t156) * mrSges(7,3) + m(7) * (t117 * t118 + t31 * t38 + t32 * t39) + m(6) * (t164 * t165 + t86 * t92 + t87 * t93) + (t165 * t255 + t86 * mrSges(6,1) + t243 * mrSges(5,1) - t87 * mrSges(6,2) + Ifges(7,5) * t334 + Ifges(7,6) * t337 - t32 * mrSges(7,2) + t31 * mrSges(7,1) - t212 * t127 / 0.2e1 + t128 * t327 + (t302 / 0.2e1 - t305 / 0.2e1 - Ifges(5,4)) * t201) * t201 + (-t164 * t255 + t249 * mrSges(6,3) + (Ifges(6,3) - Ifges(5,1) + Ifges(5,2) + Ifges(7,3) + Ifges(6,1) * t269 + (t316 - t306 / 0.2e1) * t212) * t201 + t240 + (Ifges(5,4) - t302 + t305) * t198) * t198;
t250 = t153 * t31 + t156 * t32;
t6 = (-t152 * t39 - t231 * t38 + t250) * t343 - t291 / 0.2e1 - t295 / 0.2e1 + (t153 * t335 - t231 * t336) * mrSges(7,3) + ((t165 + t248) * t346 + t237 + t350) * t198 + (-m(7) * t117 / 0.4e1 - m(6) * (t164 + t247) / 0.4e1) * t348;
t252 = t1 * qJD(1) - t6 * qJD(2);
t89 = Ifges(7,2) * t152 - t150;
t90 = -Ifges(7,1) * t231 + t315;
t5 = t275 * t330 + t31 * t107 - t32 * t108 + t117 * t88 - (-t32 * mrSges(7,3) + t90 / 0.2e1 - t77 / 0.2e1) * t152 - (-t31 * mrSges(7,3) + t79 / 0.2e1 + t89 / 0.2e1) * t231;
t251 = t5 * qJD(1) + t8 * qJD(2);
t7 = t291 + t295 + (mrSges(5,3) * t198 - t282 + t284) * t198 + ((mrSges(5,3) + t255) * t201 + t253) * t201 + m(7) * (t117 * t201 + t250) + m(6) * (t198 * t248 + t289) + m(5) * (-t165 * t198 + t289) + (m(4) * t205 + mrSges(4,3)) * (t213 ^ 2 + t215 ^ 2);
t246 = -t7 * qJD(1) - t27 * qJD(2);
t47 = 0.2e1 * t334 * mrSges(7,1) + 0.2e1 * t337 * mrSges(7,2);
t244 = qJD(1) * t47 + qJD(4) * t170;
t207 = -pkin(5) * t214 - pkin(4);
t242 = (-t273 * t297 - t325) * t345 + (t153 * t177 + t156 * t178 + t201 * t207) * t342;
t11 = (t356 - t360) * mrSges(7,3) + t238 + t277;
t234 = t11 * qJD(1);
t224 = (t212 * t93 + t214 * t92) * t346 + (t200 * t39 - t245 * t38) * t343;
t13 = 0.2e1 * t355 + (-mrSges(5,1) + t308 / 0.2e1 + t309 / 0.2e1 + t307 / 0.2e1 - t303 / 0.2e1 + t261) * t201 + t224 + t242 + t349;
t233 = -t13 * qJD(1) + qJD(2) * t265;
t17 = -t231 * t107 + m(7) * (t152 * t31 - t231 * t32) + t152 * t108 + (-m(6) * t249 - t212 * t168 - t214 * t169) * t201;
t232 = t17 * qJD(1);
t41 = t365 + (-t259 / 0.4e1 + t268) * t348;
t55 = (t200 ^ 2 + t245 ^ 2) * mrSges(7,3) + t273 * mrSges(6,3) + m(7) * (-t177 * t200 - t178 * t245) + qJ(5) * t259;
t222 = t165 * t345 - t350;
t9 = (-t351 / 0.2e1 + t360) * mrSges(7,3) + t222 + t359;
t228 = -t9 * qJD(1) - t41 * qJD(2) + t55 * qJD(4);
t26 = t229 + t230;
t226 = t6 * qJD(1) - t26 * qJD(2) - t362 / 0.2e1;
t21 = t262 - t277;
t173 = -Ifges(7,2) * t200 - t195;
t174 = -Ifges(7,2) * t245 + t314;
t175 = -Ifges(7,1) * t245 - t314;
t25 = t207 * t170 + (t175 / 0.2e1 - t174 / 0.2e1) * t200 - (t333 + t173 / 0.2e1) * t245;
t220 = -(t79 / 0.4e1 + t89 / 0.4e1) * t245 + (t90 / 0.4e1 - t77 / 0.4e1) * t200 - (t177 * t341 + t176 / 0.4e1 + t173 / 0.4e1) * t231 - (mrSges(7,3) * t332 + t175 / 0.4e1 - t174 / 0.4e1) * t152 + t117 * t170 / 0.2e1 + t177 * t357 + t108 * t332 + t198 * t274 / 0.4e1 + t207 * t88 / 0.2e1;
t223 = -Ifges(7,3) * t201 / 0.2e1 - t38 * mrSges(7,1) / 0.2e1 + t39 * mrSges(7,2) / 0.2e1 - t240;
t4 = t220 + t223;
t225 = -t4 * qJD(1) - t21 * qJD(2) - t25 * qJD(4);
t44 = -t257 * t345 + t256 + t365;
t43 = t364 + t201 * t259 / 0.2e1 + t256;
t22 = t262 + t277;
t14 = (-mrSges(7,2) * t201 + mrSges(7,3) * t153) * t329 + (mrSges(7,1) * t201 - mrSges(7,3) * t156) * t331 + (-mrSges(6,2) * t201 + mrSges(6,3) * t283) * t328 + (t201 * mrSges(6,1) + mrSges(6,3) * t281) * t327 + t261 * t201 + (-t290 / 0.2e1 - t294 / 0.2e1) * mrSges(7,3) + t355 - t224 + t242;
t12 = t341 * t351 - t238 + t277 - t361;
t10 = mrSges(7,3) * t356 + t222 - t359 + t361;
t3 = t220 - t223;
t2 = t27 * qJD(3) - t6 * qJD(4) + t8 * qJD(6);
t15 = [qJD(3) * t7 + qJD(4) * t1 + qJD(5) * t17 + qJD(6) * t5, t2, t14 * qJD(4) + t44 * qJD(5) + t12 * qJD(6) - t246 + t362, t14 * qJD(3) + t10 * qJD(5) + t3 * qJD(6) + ((t118 * t207 + t177 * t38 + t178 * t39) * t342 + (-pkin(4) * t165 + qJ(5) * t247) * t345) * t347 + (t177 * mrSges(7,1) - t178 * mrSges(7,2) + Ifges(6,5) * t328 + Ifges(7,5) * t329 + Ifges(6,6) * t327 + Ifges(7,6) * t331 - qJ(5) * t255 - Ifges(5,6)) * t271 + t252 + (t164 * mrSges(5,2) + t118 * t172 + t127 * t327 + t128 * t328 + t156 * t333 + t174 * t338 + t207 * t254 + t78 * t329 + t76 * t331 + (t153 * t178 - t156 * t177 - t200 * t38 - t245 * t39) * mrSges(7,3) + (pkin(4) * t255 - Ifges(5,5) - t214 * (Ifges(6,1) * t212 + t316) / 0.2e1 + (Ifges(6,2) * t214 + t317) * t328) * t198 + t353 * t165 + t247 * mrSges(6,3)) * qJD(4), t44 * qJD(3) + t10 * qJD(4) + t232, t12 * qJD(3) + t3 * qJD(4) + (-mrSges(7,1) * t32 - mrSges(7,2) * t31 + t275) * qJD(6) + t251; t2, t26 * qJD(4), t296 + t322 (-t273 * t318 + t349) * qJD(4) + t43 * qJD(5) + t22 * qJD(6) + (t172 + t353) * t271 + t242 * t347 - t226, t43 * qJD(4), t22 * qJD(4) - t88 * qJD(6) + t300; -qJD(4) * t13 - qJD(5) * t42 - qJD(6) * t11 + t246, -t296 + t322, 0, t233, -t278, -t234 - t270; qJD(3) * t13 - qJD(5) * t9 + qJD(6) * t4 - t252, -t41 * qJD(5) + t21 * qJD(6) + t226, -t233, qJD(5) * t55 + qJD(6) * t25, t228 (-mrSges(7,1) * t178 - mrSges(7,2) * t177 + t274) * qJD(6) - t225; t42 * qJD(3) + t9 * qJD(4) + t47 * qJD(6) - t232, t41 * qJD(4), t278, -t228 + t270, 0, t244; qJD(3) * t11 - qJD(4) * t4 - qJD(5) * t47 - t251, -t21 * qJD(4) - t300, t234, -t170 * qJD(5) + t225, -t244, 0;];
Cq  = t15;
