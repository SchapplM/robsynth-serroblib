% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
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
% Datum: 2019-03-08 21:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRRPPR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR3_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR3_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:08:15
% EndTime: 2019-03-08 21:08:22
% DurationCPUTime: 3.16s
% Computational Cost: add. (4115->383), mult. (9329->538), div. (0->0), fcn. (8173->8), ass. (0->211)
t211 = sin(qJ(6));
t203 = t211 ^ 2;
t214 = cos(qJ(6));
t205 = t214 ^ 2;
t276 = t203 + t205;
t326 = Ifges(7,4) * t214;
t251 = t211 * Ifges(7,1) + t326;
t335 = -t214 / 0.2e1;
t360 = t251 * t335;
t212 = sin(qJ(3));
t215 = cos(qJ(3));
t353 = t212 ^ 2 + t215 ^ 2;
t207 = sin(pkin(6));
t216 = cos(qJ(2));
t296 = t207 * t216;
t268 = -t296 / 0.2e1;
t358 = mrSges(4,1) + mrSges(5,1);
t357 = mrSges(4,2) - mrSges(5,3);
t356 = pkin(8) - qJ(5);
t195 = t214 * mrSges(7,1);
t322 = t211 * mrSges(7,2);
t354 = -t195 + t322;
t355 = t354 - mrSges(6,1);
t324 = Ifges(7,6) * t211;
t325 = Ifges(7,5) * t214;
t238 = t325 / 0.2e1 - t324 / 0.2e1;
t352 = Ifges(5,5) - Ifges(4,4) - Ifges(6,4) + t238;
t315 = t214 * mrSges(7,2);
t323 = t211 * mrSges(7,1);
t253 = t315 + t323;
t140 = t253 * t215;
t163 = t356 * t212;
t167 = t356 * t215;
t300 = t167 * t215;
t351 = -m(6) * (t163 * t212 + t300) + t215 * t140;
t193 = t215 * qJ(4);
t217 = -pkin(3) - pkin(4);
t145 = t212 * t217 + t193;
t285 = t214 * t215;
t320 = t212 * mrSges(7,1);
t149 = mrSges(7,3) * t285 + t320;
t286 = t214 * t149;
t291 = t211 * t215;
t319 = t212 * mrSges(7,2);
t147 = mrSges(7,3) * t291 - t319;
t293 = t211 * t147;
t328 = mrSges(7,3) * t215;
t350 = -t286 / 0.2e1 - t293 / 0.2e1 + t276 * t328 / 0.2e1;
t213 = sin(qJ(2));
t297 = t207 * t215;
t307 = cos(pkin(6));
t136 = t307 * t212 + t213 * t297;
t349 = 0.2e1 * t136;
t348 = 2 * qJD(3);
t347 = m(5) / 0.2e1;
t346 = m(6) / 0.2e1;
t345 = m(7) / 0.2e1;
t344 = -mrSges(7,1) / 0.2e1;
t343 = -mrSges(7,2) / 0.2e1;
t298 = t207 * t213;
t135 = t212 * t298 - t307 * t215;
t80 = t135 * t214 + t211 * t296;
t342 = -t80 / 0.2e1;
t341 = t136 / 0.2e1;
t312 = t215 * mrSges(7,1);
t318 = t212 * mrSges(7,3);
t148 = -t214 * t318 + t312;
t340 = -t148 / 0.2e1;
t339 = t253 / 0.2e1;
t202 = -pkin(9) + t217;
t338 = t202 / 0.2e1;
t337 = -t211 / 0.2e1;
t336 = t211 / 0.2e1;
t334 = t214 / 0.2e1;
t333 = t215 / 0.2e1;
t332 = m(6) * t145;
t208 = qJ(4) + pkin(5);
t331 = m(7) * t208;
t329 = mrSges(6,3) - mrSges(5,2);
t327 = Ifges(7,4) * t211;
t284 = t214 * t216;
t79 = -t135 * t211 + t207 * t284;
t321 = t211 * t79;
t317 = t212 * Ifges(7,5);
t316 = t212 * Ifges(7,6);
t314 = t214 * Ifges(7,2);
t313 = t214 * t80;
t311 = t215 * mrSges(7,2);
t94 = pkin(5) * t215 + t202 * t212 + t193;
t63 = -t167 * t211 + t214 * t94;
t310 = t63 * t211;
t64 = t167 * t214 + t211 * t94;
t309 = t64 * t214;
t244 = -t313 + t321;
t9 = m(7) * (-t135 - t244) * t136;
t308 = t9 * qJD(1);
t290 = t212 * t216;
t103 = (-t211 * t290 - t213 * t214) * t207;
t104 = (-t211 * t213 + t212 * t284) * t207;
t299 = t207 ^ 2 * t213;
t302 = t135 * t212;
t272 = t207 * t290;
t86 = t135 * t272;
t271 = t215 * t296;
t87 = t136 * t271;
t10 = m(6) * (t86 + (t136 * t297 - t299) * t216) + m(7) * (t103 * t79 + t104 * t80 + t87) + m(5) * (-t216 * t299 + t86 + t87) + m(4) * (t87 + (t207 * t302 - t299) * t216);
t306 = t10 * qJD(1);
t261 = t276 * t212;
t264 = -t205 / 0.2e1 - t203 / 0.2e1;
t219 = -t264 * t318 - t145 * t346 + (-t202 * t261 - t208 * t215) * t345 + t354 * t333;
t146 = -t211 * t318 - t311;
t223 = t332 / 0.2e1 + (t211 * t64 + t214 * t63) * t345 + t146 * t336 + t148 * t334;
t196 = t215 * mrSges(6,1);
t277 = t212 * mrSges(6,2) + t196;
t13 = t219 - t223 - t277;
t305 = t13 * qJD(2);
t304 = t135 * qJ(4);
t301 = t136 * t215;
t295 = t211 * t103;
t250 = -Ifges(7,2) * t211 + t326;
t231 = t250 * t215;
t114 = -t231 + t316;
t294 = t211 * t114;
t292 = t211 * t149;
t289 = t214 * t104;
t187 = Ifges(7,4) * t291;
t116 = -Ifges(7,1) * t285 + t187 + t317;
t288 = t214 * t116;
t287 = t214 * t147;
t255 = t319 / 0.2e1 - t147 / 0.2e1;
t256 = t320 / 0.2e1 + t149 / 0.2e1;
t257 = mrSges(7,3) * t264;
t22 = -t255 * t211 + t256 * t214 + t215 * t257;
t282 = t22 * qJD(2);
t26 = t256 * t211 + t255 * t214;
t281 = t26 * qJD(2);
t65 = (-0.2e1 * t264 * m(7) + m(6)) * t212;
t280 = t65 * qJD(2);
t278 = Ifges(7,5) * t291 + Ifges(7,6) * t285;
t158 = -t215 * pkin(3) - t212 * qJ(4) - pkin(2);
t273 = Ifges(7,2) / 0.4e1 - Ifges(7,1) / 0.4e1;
t269 = mrSges(6,1) / 0.2e1 - t354 / 0.2e1;
t267 = t296 / 0.2e1;
t263 = -m(5) * qJ(4) - mrSges(5,3);
t262 = m(5) * t158 - t215 * mrSges(5,1) - t212 * mrSges(5,3);
t141 = t215 * pkin(4) - t158;
t258 = t276 * t345;
t252 = Ifges(7,1) * t214 - t327;
t249 = t324 - t325;
t139 = t253 * t212;
t166 = pkin(3) * t212 - t193;
t93 = pkin(5) * t212 + pkin(9) * t215 + t141;
t61 = -t163 * t211 + t214 * t93;
t62 = t163 * t214 + t211 * t93;
t246 = t211 * t61 - t214 * t62;
t218 = -m(7) * (-t135 * t167 + t63 * t79 + t64 * t80 + (-t163 - t246) * t136) / 0.2e1 - t135 * t140 / 0.2e1 + t79 * t340 + t146 * t342 + m(5) * t166 * t267 + t268 * t332;
t159 = qJ(4) * t271;
t220 = (t295 / 0.2e1 - t289 / 0.2e1) * mrSges(7,3) + (-pkin(3) * t272 + t159) * t347 + (t217 * t272 + t159) * t346;
t241 = t289 - t295;
t232 = m(7) * t241;
t227 = t232 / 0.2e1;
t234 = t292 / 0.2e1 - t287 / 0.2e1;
t2 = t196 * t268 + (-t139 / 0.2e1 + t234) * t136 + t202 * t227 + (t331 / 0.2e1 + t269) * t271 + t218 + t220;
t113 = t215 * Ifges(7,6) + t250 * t212;
t115 = t215 * Ifges(7,5) + t252 * t212;
t168 = mrSges(6,1) * t212 - mrSges(6,2) * t215;
t3 = t163 * t140 + t167 * t139 + t145 * t168 + t61 * t148 + t63 * t149 + t62 * t146 + t64 * t147 + t262 * t166 + m(7) * (-t163 * t167 + t61 * t63 + t62 * t64) + (-pkin(2) * mrSges(4,2) - t158 * mrSges(5,3) + t113 * t336 + t115 * t335 - t352 * t215) * t215 + (-pkin(2) * mrSges(4,1) + t158 * mrSges(5,1) + t288 / 0.2e1 - t294 / 0.2e1 + (Ifges(7,3) + Ifges(6,2) + Ifges(5,1) - Ifges(5,3) - Ifges(4,2) + Ifges(4,1) - Ifges(6,1)) * t215 + t352 * t212) * t212 + (t277 + t332) * t141;
t248 = -t2 * qJD(1) + t3 * qJD(2);
t138 = t354 * t215;
t222 = (-t321 / 0.2e1 + t313 / 0.2e1) * t328 + t138 * t341 + t79 * t147 / 0.2e1 + t149 * t342;
t237 = t103 * t344 + t104 * mrSges(7,2) / 0.2e1;
t7 = t222 + t237;
t8 = t61 * t147 - t62 * t149 + t167 * t138 + t212 * t278 / 0.2e1 + (t215 * t360 - t246 * mrSges(7,3) + t114 * t334 + (t314 * t215 + t116 + t187) * t336) * t215;
t247 = t7 * qJD(1) + t8 * qJD(2);
t245 = t309 - t310;
t17 = (-t287 + t292) * t212 + t353 * mrSges(6,3) + m(7) * (t246 * t212 - t300) + t351;
t221 = (-t301 - t302) * t346 + (t244 * t212 - t301) * t345;
t224 = (t103 * t214 + t104 * t211) * t345 - m(6) * t298 / 0.2e1;
t18 = -t221 + t224;
t243 = qJD(1) * t18 - qJD(2) * t17;
t230 = m(6) * t141 + t168 - t262;
t20 = (m(7) * (t211 * t62 + t214 * t61) + t293 + t286 + t230) * t212;
t228 = (t211 * t80 + t214 * t79) * t345;
t29 = -t232 / 0.2e1 + t212 * t228;
t242 = qJD(1) * t29 + qJD(2) * t20;
t239 = t63 * mrSges(7,1) / 0.2e1 + t64 * t343;
t236 = t315 / 0.2e1 + t323 / 0.2e1;
t235 = t167 * t339 - t208 * t138 / 0.2e1;
t169 = -t314 - t327;
t233 = t169 * t336 - t360;
t15 = (t339 - t236) * t136;
t4 = (Ifges(7,3) / 0.2e1 + t202 * t257) * t215 + (0.3e1 / 0.4e1 * t317 + t116 / 0.4e1 + t187 / 0.4e1 + t149 * t338 + (-t169 / 0.4e1 + t273 * t214) * t215) * t214 + (-0.3e1 / 0.4e1 * t316 - t114 / 0.4e1 + t147 * t338 + (0.3e1 / 0.4e1 * t326 + t251 / 0.4e1 - t273 * t211) * t215) * t211 + t235 + t239;
t43 = -t208 * t253 + t250 * t334 + t252 * t336 + t233;
t226 = t15 * qJD(1) + t4 * qJD(2) - t43 * qJD(3);
t24 = (t146 / 0.2e1 + t311 / 0.2e1) * t214 + (t340 + t312 / 0.2e1) * t211 + 0.2e1 * (-t310 / 0.4e1 + t309 / 0.4e1 - t167 / 0.4e1) * m(7);
t46 = (0.1e1 / 0.4e1 - t203 / 0.4e1 - t205 / 0.4e1) * m(7) * t349;
t90 = -m(6) * qJ(4) + t263 - t331 + t355;
t225 = qJD(1) * t46 - qJD(2) * t24 - qJD(3) * t90;
t70 = t212 * t258 - t261 * t345;
t31 = t136 * t258 + (t346 + m(7) / 0.4e1 + t347) * t349;
t27 = t236 * t212 - t234;
t25 = t227 + (t228 + (m(5) + m(6)) * t296) * t212;
t23 = (t195 / 0.2e1 - t322 / 0.2e1) * t212 + t350;
t21 = t285 * t343 + t291 * t344 + t146 * t334 + t148 * t337 + t245 * t345 + (m(5) * pkin(8) - t329) * t215 + (t345 + m(6)) * t167;
t19 = t221 + t224;
t16 = -t236 * t136 - t253 * t341;
t14 = t219 + t223;
t6 = t222 - t237;
t5 = -t288 / 0.4e1 + t294 / 0.4e1 - t214 * (Ifges(7,2) * t285 + t187) / 0.4e1 - t211 * t231 / 0.4e1 + Ifges(7,3) * t333 - t235 + t239 - t251 * t291 / 0.2e1 + (t249 / 0.4e1 + t238) * t212 + t350 * t202 + (t252 + t169) * t285 / 0.4e1;
t1 = t277 * t267 - t136 * t292 / 0.2e1 + ((-mrSges(4,1) / 0.2e1 - mrSges(5,1) / 0.2e1 + mrSges(6,2) / 0.2e1) * t212 + (-mrSges(4,2) / 0.2e1 + mrSges(5,3) / 0.2e1 + t269) * t215) * t296 + (t241 * t202 + t208 * t271) * t345 - t218 + t220 + (t139 + t287) * t341 + (t212 * t358 + t215 * t357) * t268;
t11 = [t10 * qJD(2) + t9 * qJD(3), t1 * qJD(3) + t25 * qJD(4) + t19 * qJD(5) + t6 * qJD(6) + t306 + (t103 * t149 + t104 * t147 + 0.2e1 * (t103 * t61 + t104 * t62) * t345 + 0.2e1 * (t347 + m(4) / 0.2e1) * t353 * pkin(8) * t296 + ((-m(4) * pkin(2) - mrSges(4,1) * t215 + mrSges(4,2) * t212 - mrSges(3,1) - t230) * t213 + (m(7) * t300 - mrSges(3,2) + t353 * (mrSges(4,3) - t329) - t351) * t216) * t207) * qJD(2), t308 + t1 * qJD(2) + t31 * qJD(4) + t16 * qJD(6) + ((t136 * t217 - t304) * t346 + (t276 * t202 * t136 - t135 * t208) * t345 + (-pkin(3) * t136 - t304) * t347) * t348 + ((t355 + t357) * t135 + (-t276 * mrSges(7,3) + mrSges(6,2) - t358) * t136) * qJD(3), qJD(2) * t25 + qJD(3) * t31, qJD(2) * t19, t6 * qJD(2) + t16 * qJD(3) + (-mrSges(7,1) * t80 - mrSges(7,2) * t79) * qJD(6); -qJD(3) * t2 + qJD(4) * t29 - qJD(5) * t18 + qJD(6) * t7 - t306, qJD(3) * t3 + qJD(4) * t20 + qJD(5) * t17 + qJD(6) * t8, t21 * qJD(4) + t14 * qJD(5) + t5 * qJD(6) + ((-qJ(4) * t163 + t167 * t217) * t346 + (-t163 * t208 + t245 * t202) * t345) * t348 + t248 + (t167 * mrSges(6,2) + t208 * t139 + (-t113 / 0.2e1 + t202 * t146 - t64 * mrSges(7,3)) * t214 + (-t115 / 0.2e1 - t202 * t148 + t63 * mrSges(7,3)) * t211 + (Ifges(5,4) + Ifges(6,6) + Ifges(4,5) + Ifges(7,5) * t337 + Ifges(7,6) * t335 - pkin(3) * mrSges(5,2) - t217 * mrSges(6,3) + (-m(5) * pkin(3) - t358) * pkin(8)) * t215 + (Ifges(5,6) - Ifges(6,5) - Ifges(4,6) + t329 * qJ(4) + (mrSges(4,2) + t263) * pkin(8) - t233) * t212 + t355 * t163) * qJD(3), qJD(3) * t21 + qJD(5) * t70 + qJD(6) * t23 + t242, qJD(3) * t14 + qJD(4) * t70 + qJD(6) * t27 - t243, t5 * qJD(3) + t23 * qJD(4) + t27 * qJD(5) + (-mrSges(7,1) * t62 - mrSges(7,2) * t61 + t278) * qJD(6) + t247; qJD(2) * t2 + qJD(4) * t46 - qJD(6) * t15 - t308, -qJD(4) * t24 + qJD(5) * t13 - qJD(6) * t4 - t248, -qJD(4) * t90 + qJD(6) * t43, t225, t305 (t202 * t354 + t249) * qJD(6) - t226; -qJD(2) * t29 - qJD(3) * t46, qJD(3) * t24 - qJD(5) * t65 - qJD(6) * t22 - t242, -t225, 0, -t280, qJD(6) * t354 - t282; qJD(2) * t18, -qJD(3) * t13 + qJD(4) * t65 - qJD(6) * t26 + t243, -t305, t280, 0, -qJD(6) * t253 - t281; -t7 * qJD(2) + t15 * qJD(3), qJD(3) * t4 + qJD(4) * t22 + qJD(5) * t26 - t247, t226, t282, t281, 0;];
Cq  = t11;
