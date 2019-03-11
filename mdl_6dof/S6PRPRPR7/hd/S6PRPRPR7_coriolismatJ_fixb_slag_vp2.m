% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
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
% Datum: 2019-03-08 19:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRPRPR7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR7_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR7_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR7_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:51:17
% EndTime: 2019-03-08 19:51:21
% DurationCPUTime: 2.55s
% Computational Cost: add. (3798->362), mult. (8840->512), div. (0->0), fcn. (7794->8), ass. (0->210)
t210 = cos(pkin(6));
t212 = sin(qJ(4));
t215 = cos(qJ(4));
t209 = sin(pkin(6));
t216 = cos(qJ(2));
t293 = t209 * t216;
t139 = t210 * t215 - t212 * t293;
t320 = m(7) * t139;
t138 = t210 * t212 + t215 * t293;
t211 = sin(qJ(6));
t214 = cos(qJ(6));
t213 = sin(qJ(2));
t294 = t209 * t213;
t81 = t138 * t211 + t214 * t294;
t297 = t81 * t211;
t80 = t138 * t214 - t211 * t294;
t298 = t80 * t214;
t347 = -t138 + t297 + t298;
t206 = t212 ^ 2;
t208 = t215 ^ 2;
t346 = t206 + t208;
t198 = qJ(5) * t212;
t172 = pkin(4) * t215 + t198;
t345 = m(6) * t172;
t319 = m(7) * t212;
t218 = -pkin(2) - pkin(8);
t344 = -pkin(5) + t218;
t323 = -t214 / 0.2e1;
t326 = -t211 / 0.2e1;
t343 = Ifges(7,5) * t326 + Ifges(7,6) * t323 + Ifges(5,4) + Ifges(6,6);
t307 = t214 * mrSges(7,2);
t315 = t211 * mrSges(7,1);
t170 = t307 + t315;
t342 = t170 + mrSges(6,3);
t197 = t215 * qJ(5);
t278 = -t212 * pkin(4) + t197;
t171 = -t212 * mrSges(6,2) - t215 * mrSges(6,3);
t205 = t211 ^ 2;
t207 = t214 ^ 2;
t277 = t205 + t207;
t282 = t214 * t215;
t103 = (-t211 * t216 - t213 * t282) * t209;
t285 = t214 * t103;
t286 = t213 * t215;
t104 = (-t211 * t286 + t214 * t216) * t209;
t292 = t211 * t104;
t249 = t285 + t292;
t149 = pkin(9) * t215 + t172;
t165 = t344 * t212;
t76 = -t149 * t211 + t165 * t214;
t306 = t214 * t76;
t77 = t149 * t214 + t165 * t211;
t313 = t211 * t77;
t341 = -t313 - t306;
t164 = qJ(3) - t278;
t287 = t212 * t214;
t301 = t215 * mrSges(7,2);
t152 = mrSges(7,3) * t287 - t301;
t283 = t214 * t152;
t289 = t211 * t212;
t304 = t215 * mrSges(7,1);
t150 = -mrSges(7,3) * t289 + t304;
t291 = t211 * t150;
t141 = pkin(9) * t212 + t164;
t166 = t344 * t215;
t68 = -t141 * t211 - t166 * t214;
t69 = t141 * t214 - t166 * t211;
t340 = -m(6) * t164 + m(7) * (t211 * t68 - t214 * t69) - t171 - t283 + t291;
t339 = 2 * qJD(4);
t338 = m(6) / 0.2e1;
t337 = m(7) / 0.2e1;
t336 = -mrSges(7,1) / 0.2e1;
t335 = mrSges(7,2) / 0.2e1;
t263 = (0.1e1 - t277) * t215;
t333 = t263 * t319;
t144 = t212 * t170;
t332 = t144 / 0.2e1;
t331 = -t150 / 0.2e1;
t288 = t211 * t215;
t312 = t212 * mrSges(7,1);
t151 = -mrSges(7,3) * t288 - t312;
t330 = t151 / 0.2e1;
t329 = t152 / 0.2e1;
t310 = t212 * mrSges(7,2);
t153 = mrSges(7,3) * t282 + t310;
t328 = t153 / 0.2e1;
t308 = t214 * mrSges(7,1);
t314 = t211 * mrSges(7,2);
t167 = t308 - t314;
t327 = t167 / 0.2e1;
t325 = t211 / 0.2e1;
t324 = t212 / 0.2e1;
t322 = t214 / 0.2e1;
t16 = (t139 * t263 + t212 * t347) * t337;
t318 = t16 * qJD(4);
t317 = Ifges(7,4) * t211;
t316 = Ifges(7,4) * t214;
t311 = t212 * mrSges(5,2);
t309 = t212 * mrSges(7,3);
t305 = t215 * mrSges(5,1);
t303 = t215 * mrSges(5,2);
t302 = t215 * mrSges(6,2);
t300 = t215 * Ifges(7,5);
t299 = t215 * Ifges(7,6);
t296 = qJ(5) * t138;
t178 = t209 ^ 2 * t213 * t216;
t273 = t212 * t294;
t87 = t139 * t273;
t13 = m(7) * (t103 * t80 + t104 * t81 + t87) + m(5) * (-t138 * t209 * t286 + t178 + t87) + m(6) * (t178 + (-t138 * t215 + t139 * t212) * t294);
t295 = t13 * qJD(1);
t290 = t211 * t153;
t284 = t214 * t151;
t237 = t291 / 0.2e1 - t283 / 0.2e1;
t266 = t205 / 0.2e1 + t207 / 0.2e1;
t247 = t266 * t309;
t221 = (t247 + t237) * t215 + t144 * t324;
t242 = t307 / 0.2e1 + t315 / 0.2e1;
t22 = t221 + t242;
t281 = t22 * qJD(2);
t28 = (t301 / 0.2e1 - t152 / 0.2e1) * t214 + (t304 / 0.2e1 + t150 / 0.2e1) * t211 + t247;
t280 = t28 * qJD(2);
t279 = t346 * t218 * t294;
t275 = qJD(4) * t212;
t274 = t170 * qJD(6);
t269 = t294 / 0.2e1;
t174 = -Ifges(7,2) * t211 + t316;
t259 = Ifges(7,1) * t211 + t316;
t268 = -t259 / 0.4e1 - t174 / 0.4e1;
t176 = Ifges(7,1) * t214 - t317;
t258 = Ifges(7,2) * t214 + t317;
t267 = t176 / 0.4e1 - t258 / 0.4e1;
t265 = m(7) * t277;
t217 = -pkin(4) - pkin(9);
t264 = t277 * t217;
t261 = -t277 * mrSges(7,3) - mrSges(5,1);
t260 = -t212 * mrSges(5,1) - t303;
t220 = (-t298 / 0.2e1 - t297 / 0.2e1) * t309 + t139 * t332 + t80 * t329 + t81 * t331;
t244 = t103 * t336 + t104 * t335;
t7 = t220 + t244;
t120 = t258 * t212 + t299;
t192 = Ifges(7,4) * t287;
t122 = Ifges(7,1) * t289 + t192 + t300;
t145 = -Ifges(7,2) * t289 + t192;
t146 = t212 * t176;
t191 = Ifges(7,5) * t287;
t8 = t68 * t152 - t69 * t150 + t215 * t191 / 0.2e1 + t165 * t144 + ((t122 / 0.2e1 + t145 / 0.2e1 - t68 * mrSges(7,3)) * t214 + (t146 / 0.2e1 - t120 / 0.2e1 - t69 * mrSges(7,3) - t299 / 0.2e1) * t211) * t212;
t257 = t7 * qJD(1) + t8 * qJD(2);
t253 = t211 * t80 - t214 * t81;
t177 = t206 * t294;
t223 = (-t249 * t215 + t177) * t337 + 0.2e1 * (m(5) / 0.4e1 + m(6) / 0.4e1) * (t208 * t294 + t177);
t239 = m(7) * t253;
t20 = (-m(5) / 0.2e1 - m(6) / 0.2e1) * t294 + t239 / 0.2e1 + t223;
t245 = mrSges(4,3) - t260;
t26 = (m(5) + m(4)) * qJ(3) + t245 - t340;
t252 = qJD(1) * t20 - qJD(2) * t26;
t27 = t340 * t215;
t233 = m(7) * t249;
t228 = -t233 / 0.2e1;
t231 = t253 * t337;
t32 = t215 * t231 + t228;
t251 = qJD(1) * t32 + qJD(2) * t27;
t12 = t347 * t320;
t250 = t12 * qJD(1) + t16 * qJD(3);
t248 = t266 * t217 * mrSges(7,3);
t246 = t77 * t335 + t76 * t336;
t243 = -t314 / 0.2e1 + t308 / 0.2e1;
t168 = t212 * mrSges(6,3) - t302;
t169 = t305 - t311;
t241 = t169 / 0.2e1 + t168 / 0.2e1 + t345 / 0.2e1;
t238 = qJ(5) * t332 + t165 * t327;
t236 = -t290 / 0.2e1 - t284 / 0.2e1;
t235 = -t120 / 0.4e1 + t146 / 0.4e1 + t217 * t329;
t234 = -t122 / 0.4e1 - t145 / 0.4e1 + t217 * t331;
t142 = t212 * t167;
t143 = t167 * t215;
t222 = (t211 * t69 + t214 * t68 + t166) * t337 + t152 * t325 + t150 * t322 - t143 / 0.2e1;
t219 = t222 * t139 + (-t138 * t165 + t76 * t80 + t77 * t81) * t337 + t138 * t142 / 0.2e1 + t80 * t330 + t81 * t328;
t1 = (t285 / 0.2e1 + t292 / 0.2e1) * mrSges(7,3) + t217 * t228 + ((-mrSges(5,1) / 0.2e1 + mrSges(6,2) / 0.2e1) * t215 + (-t170 / 0.2e1 + mrSges(5,2) / 0.2e1 - mrSges(6,3) / 0.2e1) * t212 - m(7) * t198 / 0.2e1 - t345 / 0.2e1 + t241) * t294 + t219;
t10 = ((t165 + t341) * t337 - t142 / 0.2e1 + t236) * t215 + t222 * t212;
t121 = -t212 * Ifges(7,6) + t258 * t215;
t123 = -t212 * Ifges(7,5) + t259 * t215;
t3 = -t165 * t143 - t166 * t142 + qJ(3) * t169 + t172 * t171 + t76 * t150 + t68 * t151 + t77 * t152 + t69 * t153 + m(7) * (t165 * t166 + t68 * t76 + t69 * t77) + (t168 + t345) * t164 + (t120 * t322 + t122 * t325 - t215 * t343) * t215 + (t121 * t322 + t123 * t325 + (-Ifges(7,3) - Ifges(6,2) + Ifges(6,3) - Ifges(5,1) + Ifges(5,2)) * t215 + t343 * t212) * t212;
t229 = t1 * qJD(1) + t3 * qJD(2) + t10 * qJD(3);
t227 = -t167 / 0.2e1 + t243;
t226 = t16 * qJD(1) + t10 * qJD(2) + qJD(3) * t333;
t17 = t227 * t139;
t47 = qJ(5) * t167 + (-t174 / 0.2e1 - t259 / 0.2e1) * t214 + (t258 / 0.2e1 - t176 / 0.2e1) * t211;
t5 = (Ifges(7,3) / 0.2e1 - t248) * t212 + (-0.3e1 / 0.4e1 * t299 + t267 * t212 + t235) * t214 + (-0.3e1 / 0.4e1 * t300 + t268 * t212 + t234) * t211 + t238 + t246;
t78 = t227 * t212;
t225 = t17 * qJD(1) - t5 * qJD(2) + t78 * qJD(3) - t47 * qJD(4);
t137 = (m(6) + m(7)) * qJ(5) + t342;
t31 = (t330 + t312 / 0.2e1) * t214 + (t328 - t310 / 0.2e1) * t211 + 0.2e1 * (t313 / 0.4e1 + t306 / 0.4e1 - t165 / 0.4e1) * m(7);
t50 = 0.2e1 * (t205 / 0.4e1 + t207 / 0.4e1 - 0.1e1 / 0.4e1) * t320;
t95 = (-0.1e1 / 0.2e1 + t266) * t319;
t224 = qJD(1) * t50 + qJD(2) * t31 + qJD(3) * t95 - qJD(4) * t137;
t190 = qJ(3) * t293;
t85 = t319 / 0.2e1 + (m(6) + t265 / 0.2e1) * t212;
t79 = t167 * t324 + t243 * t212;
t34 = t320 / 0.2e1 + 0.2e1 * (t338 + t265 / 0.4e1) * t139;
t30 = t233 / 0.2e1 + (-m(6) * t294 + t231) * t215;
t29 = t242 * t215 - t237 - t277 * t309 / 0.2e1;
t25 = t289 * t335 + t287 * t336 + (m(6) * t218 - mrSges(6,1)) * t212 - t236 + (t165 - t341) * t337;
t23 = t221 - t242;
t19 = m(4) * t294 - t239 / 0.2e1 + t223 + (m(5) + m(6)) * t269;
t18 = (t243 + t327) * t139;
t9 = t10 * qJD(4);
t6 = t220 - t244;
t4 = Ifges(7,5) * t288 / 0.2e1 + Ifges(7,6) * t282 / 0.2e1 + (-t299 / 0.4e1 + t235) * t214 + (-t300 / 0.4e1 + t234) * t211 + t238 - t246 + (-Ifges(7,3) / 0.2e1 + t268 * t211 + t267 * t214 - t248) * t212;
t2 = (qJ(5) * t273 + t249 * t217) * t337 + t241 * t294 + t219 - t249 * mrSges(7,3) / 0.2e1 - (t311 + t302) * t294 / 0.2e1 + (t342 * t212 + t305 + t345) * t269;
t11 = [t13 * qJD(2) + t12 * qJD(4), t19 * qJD(3) + t2 * qJD(4) + t30 * qJD(5) + t6 * qJD(6) + t295 + (t103 * t150 + t104 * t152 + ((-mrSges(3,2) + t171 + t245) * t216 + (-t212 * t142 - mrSges(3,1) + mrSges(4,2) - (mrSges(6,1) + mrSges(5,3)) * t346) * t213) * t209 + 0.2e1 * (t103 * t68 + t104 * t69 + t165 * t273) * t337 + m(5) * (t190 + t279) + 0.2e1 * (t164 * t293 + t279) * t338 + m(4) * (-pkin(2) * t294 + t190)) * qJD(2), qJD(2) * t19 + t318, t2 * qJD(2) + t34 * qJD(5) + t18 * qJD(6) + ((t139 * t264 - t296) * t337 + (-pkin(4) * t139 - t296) * t338) * t339 + t250 + ((mrSges(5,2) - t342) * t138 + (mrSges(6,2) + t261) * t139) * qJD(4), qJD(2) * t30 + qJD(4) * t34, t6 * qJD(2) + t18 * qJD(4) + (-mrSges(7,1) * t81 - mrSges(7,2) * t80) * qJD(6); -qJD(3) * t20 + qJD(4) * t1 + qJD(5) * t32 + qJD(6) * t7 - t295, qJD(3) * t26 + qJD(4) * t3 + qJD(5) * t27 + qJD(6) * t8, qJD(6) * t23 - t252 + t9, t25 * qJD(5) + t4 * qJD(6) + (pkin(4) * mrSges(6,1) + Ifges(7,5) * t323 + Ifges(7,6) * t325 + Ifges(6,4) - Ifges(5,5)) * t275 + t229 + (t123 * t322 + t121 * t326 + t166 * t170 + (t174 * t322 + t176 * t325 + Ifges(6,5) - Ifges(5,6)) * t215 + (m(6) * t278 - t171 + t260) * t218 + (m(7) * t166 - t215 * mrSges(6,1) - t143) * qJ(5) + (-m(7) * t341 + t284 + t290) * t217 + t341 * mrSges(7,3)) * qJD(4), qJD(4) * t25 + qJD(6) * t29 + t251, t23 * qJD(3) + t4 * qJD(4) + t29 * qJD(5) + (-mrSges(7,1) * t69 - mrSges(7,2) * t68 - Ifges(7,6) * t289 + t191) * qJD(6) + t257; qJD(2) * t20 + t318, qJD(6) * t22 + t252 + t9, qJD(4) * t333 (t215 * t170 - t171 - t303) * qJD(4) + t85 * qJD(5) + t79 * qJD(6) + t261 * t275 + (t278 * t338 + (t212 * t264 + t197) * t337) * t339 + t226, t85 * qJD(4), t79 * qJD(4) + t215 * t274 + t281; -qJD(2) * t1 - qJD(5) * t50 - qJD(6) * t17 - t250, -qJD(5) * t31 + qJD(6) * t5 - t229, -t95 * qJD(5) - t78 * qJD(6) - t226, qJD(5) * t137 + qJD(6) * t47, -t224 ((-mrSges(7,2) * t217 - Ifges(7,6)) * t214 + (-mrSges(7,1) * t217 - Ifges(7,5)) * t211) * qJD(6) - t225; -qJD(2) * t32 + qJD(4) * t50, qJD(4) * t31 - qJD(6) * t28 - t251, t95 * qJD(4), t224, 0, -t274 - t280; -t7 * qJD(2) + t17 * qJD(4), -qJD(3) * t22 - qJD(4) * t5 + qJD(5) * t28 - t257, t78 * qJD(4) - t281, t225, t280, 0;];
Cq  = t11;
