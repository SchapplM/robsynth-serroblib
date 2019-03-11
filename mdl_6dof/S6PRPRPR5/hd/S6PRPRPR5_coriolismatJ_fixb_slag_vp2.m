% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRPRPR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR5_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR5_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:43:13
% EndTime: 2019-03-08 19:43:19
% DurationCPUTime: 3.60s
% Computational Cost: add. (7602->388), mult. (17306->550), div. (0->0), fcn. (19091->10), ass. (0->218)
t207 = sin(pkin(11));
t209 = cos(pkin(11));
t208 = sin(pkin(6));
t212 = sin(qJ(2));
t291 = t208 * t212;
t305 = cos(pkin(6));
t180 = t207 * t305 + t209 * t291;
t211 = sin(qJ(4));
t226 = t207 * t291 - t209 * t305;
t328 = cos(qJ(4));
t110 = t180 * t211 + t226 * t328;
t210 = sin(qJ(6));
t213 = cos(qJ(6));
t214 = cos(qJ(2));
t290 = t208 * t214;
t89 = -t210 * t110 + t213 * t290;
t340 = -t89 / 0.2e1;
t190 = t207 * t328 + t211 * t209;
t366 = -t190 / 0.2e1;
t205 = t210 ^ 2;
t206 = t213 ^ 2;
t276 = t205 + t206;
t365 = mrSges(7,3) * t276;
t310 = t213 * mrSges(7,2);
t316 = t210 * mrSges(7,1);
t194 = t310 + t316;
t333 = -t194 / 0.2e1;
t354 = -t207 * t211 + t328 * t209;
t227 = t354 * t333;
t363 = m(6) + m(5);
t111 = t180 * t328 - t211 * t226;
t326 = m(7) * t111;
t324 = mrSges(5,3) + mrSges(6,1);
t364 = t190 * t324;
t361 = m(7) * t276;
t325 = -mrSges(5,1) + mrSges(6,2);
t267 = t205 / 0.2e1 + t206 / 0.2e1;
t360 = mrSges(7,3) * t267;
t359 = t194 + mrSges(6,3);
t356 = t209 * t180 + t207 * t226;
t335 = t354 / 0.2e1;
t355 = t335 * t365;
t318 = Ifges(7,6) * t213;
t320 = Ifges(7,5) * t210;
t241 = t320 / 0.2e1 + t318 / 0.2e1;
t353 = -Ifges(6,6) - Ifges(5,4) + t241;
t343 = m(7) / 0.2e1;
t344 = m(6) / 0.2e1;
t88 = t110 * t213 + t210 * t290;
t352 = ((t210 * t88 + t213 * t89) * t343 + t290 * t344) * t190;
t311 = t213 * mrSges(7,1);
t315 = t210 * mrSges(7,2);
t193 = -t311 + t315;
t132 = t193 * t354;
t351 = -t324 * t354 + t132;
t138 = -mrSges(7,3) * t190 * t210 + mrSges(7,1) * t354;
t309 = t213 * mrSges(7,3);
t140 = -mrSges(7,2) * t354 + t190 * t309;
t183 = t354 * mrSges(6,3);
t141 = -t190 * mrSges(6,2) - t183;
t184 = t354 * mrSges(5,2);
t143 = t190 * mrSges(5,1) + t184;
t303 = qJ(5) * t354;
t142 = pkin(4) * t190 - t303;
t327 = m(6) * t142;
t274 = -t327 / 0.2e1;
t350 = t88 * t138 / 0.2e1 + t140 * t340 + t110 * t132 / 0.2e1 + (t274 - t141 / 0.2e1 - t143 / 0.2e1) * t290;
t349 = 0.2e1 * t190;
t348 = 2 * qJD(4);
t347 = m(4) / 0.2e1;
t346 = m(5) / 0.2e1;
t345 = -m(6) / 0.2e1;
t342 = mrSges(7,1) / 0.2e1;
t341 = -mrSges(7,2) / 0.2e1;
t339 = pkin(4) + pkin(9);
t338 = -qJ(5) / 0.2e1;
t158 = t190 * t290;
t337 = t158 / 0.2e1;
t332 = -t210 / 0.2e1;
t331 = t210 / 0.2e1;
t330 = -t213 / 0.2e1;
t329 = t213 / 0.2e1;
t323 = pkin(8) + qJ(3);
t322 = Ifges(7,4) * t210;
t321 = Ifges(7,4) * t213;
t319 = Ifges(7,6) * t190;
t317 = t190 * Ifges(7,5);
t100 = t190 * t339 - t303;
t192 = t323 * t209;
t266 = t323 * t207;
t147 = t192 * t328 - t211 * t266;
t107 = pkin(5) * t354 + t147;
t58 = t100 * t213 + t107 * t210;
t314 = t210 * t58;
t313 = t210 * t89;
t294 = t354 * t213;
t179 = Ifges(7,4) * t294;
t295 = t354 * t210;
t94 = -Ifges(7,1) * t295 - t179 + t317;
t312 = t210 * t94;
t57 = -t100 * t210 + t107 * t213;
t308 = t213 * t57;
t307 = t213 * t88;
t255 = Ifges(7,2) * t213 + t322;
t92 = -t255 * t354 + t319;
t306 = t213 * t92;
t304 = qJ(5) * t110;
t302 = t107 * t193;
t159 = t354 * t290;
t75 = t111 * t159;
t250 = t307 - t313;
t12 = (-t110 + t250) * t326;
t300 = t12 * qJD(1);
t129 = t158 * t213 - t210 * t291;
t130 = t158 * t210 + t213 * t291;
t293 = t208 ^ 2 * t212;
t13 = m(7) * (t129 * t88 - t130 * t89 + t75) + m(4) * (t356 * t208 - t293) * t214 + t363 * (t110 * t158 - t214 * t293 + t75);
t299 = t13 * qJD(1);
t297 = t159 * qJ(5);
t263 = t276 * t339;
t217 = t274 + (-t190 * t263 + t303) * t343 - t227;
t219 = t327 / 0.2e1 + (-t210 * t57 + t213 * t58) * t343 + t138 * t332 + t140 * t329;
t16 = t183 - t184 + (-t360 + t325) * t190 + t217 - t219;
t296 = t16 * qJD(2);
t86 = t354 * t111;
t288 = t210 * t130;
t275 = mrSges(7,3) * t295;
t137 = t190 * mrSges(7,1) + t275;
t287 = t210 * t137;
t139 = -t190 * mrSges(7,2) - mrSges(7,3) * t294;
t286 = t210 * t139;
t285 = t213 * t129;
t284 = t213 * t137;
t283 = t213 * t139;
t282 = t339 * t137;
t281 = t339 * t139;
t228 = (t310 / 0.2e1 + t316 / 0.2e1) * t190;
t236 = t287 / 0.2e1 - t283 / 0.2e1;
t22 = -t354 * t360 + t228 + t236;
t280 = t22 * qJD(2);
t239 = t315 / 0.2e1 - t311 / 0.2e1;
t229 = t239 * t190;
t235 = t286 / 0.2e1 + t284 / 0.2e1;
t25 = -t229 + t235;
t279 = t25 * qJD(2);
t49 = (m(7) * t267 + t344) * t349;
t278 = t49 * qJD(2);
t277 = t207 ^ 2 + t209 ^ 2;
t273 = -Ifges(7,2) / 0.4e1 + Ifges(7,1) / 0.4e1;
t271 = -t309 / 0.2e1;
t201 = -pkin(3) * t209 - pkin(2);
t265 = t277 * mrSges(4,3);
t264 = t277 * qJ(3);
t146 = t192 * t211 + t328 * t266;
t259 = t361 / 0.4e1;
t196 = Ifges(7,1) * t213 - t322;
t256 = Ifges(7,1) * t210 + t321;
t105 = pkin(5) * t190 + t146;
t244 = -qJ(5) * t190 + t201;
t131 = -pkin(4) * t354 + t244;
t133 = t193 * t190;
t144 = mrSges(6,2) * t354 - mrSges(6,3) * t190;
t91 = -t339 * t354 + t244;
t53 = t105 * t213 - t210 * t91;
t54 = t105 * t210 + t213 * t91;
t93 = Ifges(7,6) * t354 + t190 * t255;
t95 = Ifges(7,5) * t354 + t190 * t256;
t1 = t105 * t132 + t107 * t133 + t57 * t137 + t53 * t138 + t58 * t139 + t54 * t140 + t142 * t144 + t201 * t143 + m(7) * (-t105 * t107 + t53 * t57 + t54 * t58) + (t306 / 0.2e1 + t312 / 0.2e1 + t353 * t190) * t190 + (t141 + t327) * t131 - (t93 * t329 + t95 * t331 + (Ifges(6,3) - Ifges(5,1) - Ifges(6,2) + Ifges(5,2) - Ifges(7,3)) * t190 + t353 * t354) * t354;
t247 = t285 + t288;
t218 = (-pkin(4) * t158 + t297) * t345 - m(7) * (-t247 * t339 + t297) / 0.2e1;
t224 = t133 / 0.2e1 + t235;
t231 = -t107 * t110 + t57 * t88 - t58 * t89;
t252 = t210 * t54 + t213 * t53;
t242 = -t105 + t252;
t2 = (-mrSges(6,2) / 0.2e1 + mrSges(5,1) / 0.2e1) * t158 + (t285 / 0.2e1 + t288 / 0.2e1) * mrSges(7,3) + (t333 - mrSges(6,3) / 0.2e1 + mrSges(5,2) / 0.2e1) * t159 + t224 * t111 + (t111 * t242 + t231) * t343 + t218 + t350;
t254 = t2 * qJD(1) + t1 * qJD(2);
t178 = Ifges(7,5) * t294;
t4 = t53 * t139 - t54 * t137 + t178 * t366 - ((-t53 * mrSges(7,3) + t107 * mrSges(7,2) + t94 / 0.2e1 - t179 / 0.2e1) * t213 + (-t54 * mrSges(7,3) - t319 / 0.2e1 + t107 * mrSges(7,1) - t92 / 0.2e1 - (-t322 / 0.2e1 + (Ifges(7,1) / 0.2e1 - Ifges(7,2) / 0.2e1) * t213) * t354) * t210) * t354;
t238 = -t88 * t139 / 0.2e1 + t137 * t340;
t240 = t129 * t342 + t130 * t341;
t7 = -(t111 * t333 + (t307 / 0.2e1 - t313 / 0.2e1) * mrSges(7,3)) * t354 + t238 + t240;
t253 = -t7 * qJD(1) + t4 * qJD(2);
t251 = t308 + t314;
t10 = t265 - t351 * t354 + (t284 + t286 + t364) * t190 + m(7) * (t107 * t354 + t190 * t252) + m(4) * t264 + t363 * (t146 * t190 + t147 * t354);
t216 = (t190 * t250 + t86) * t343 + t356 * t347 + (t344 + t346) * (t110 * t190 + t86);
t232 = m(7) * (-t210 * t129 + t213 * t130);
t15 = (t345 - m(5) / 0.2e1 - m(4) / 0.2e1) * t291 - t232 / 0.2e1 + t216;
t249 = -qJD(1) * t15 - qJD(2) * t10;
t18 = (m(7) * (t210 * t53 - t213 * t54) - t283 + t287 - m(6) * t131 - t144) * t190;
t220 = m(6) * t337 + t247 * t343;
t29 = -t220 + t352;
t248 = qJD(1) * t29 + qJD(2) * t18;
t246 = t146 * t158 + t147 * t159;
t243 = t341 * t58 + t342 * t57;
t234 = t138 * t330 + t140 * t332;
t195 = -Ifges(7,2) * t210 + t321;
t233 = t195 * t329 + t196 * t331;
t20 = (t193 / 0.2e1 - t239) * t111;
t5 = t302 / 0.2e1 - (-Ifges(7,3) / 0.2e1 - t339 * t360) * t354 + (0.3e1 / 0.4e1 * t317 - t179 / 0.4e1 + t94 / 0.4e1 - t282 / 0.2e1 - (mrSges(7,1) * t338 + t195 / 0.4e1 + t273 * t210) * t354) * t210 + (0.3e1 / 0.4e1 * t319 + t92 / 0.4e1 + t281 / 0.2e1 - (0.3e1 / 0.4e1 * t322 + mrSges(7,2) * t338 - t196 / 0.4e1 - t273 * t213) * t354) * t213 + t243;
t74 = -qJ(5) * t193 + t255 * t331 + t256 * t330 - t233;
t225 = t20 * qJD(1) + t5 * qJD(2) - t74 * qJD(4);
t182 = (m(6) + m(7)) * qJ(5) + t359;
t24 = -t239 * t354 + 0.2e1 * (t107 / 0.4e1 - t314 / 0.4e1 - t308 / 0.4e1) * m(7) + t234;
t40 = 0.2e1 * (t205 / 0.4e1 + t206 / 0.4e1 - 0.1e1 / 0.4e1) * t326;
t223 = qJD(1) * t40 - qJD(2) * t24 - qJD(4) * t182;
t50 = (t259 + m(6) / 0.4e1) * t349 + (m(6) + t361) * t366;
t38 = t326 / 0.2e1 + 0.2e1 * (t344 + t259) * t111;
t28 = t352 + t220;
t26 = -t229 - t235;
t23 = t228 - t236 + t355;
t21 = (-t193 / 0.2e1 - t239) * t111;
t19 = m(6) * t147 - (-mrSges(6,1) + t239) * t354 - t234 + (t251 + t107) * t343;
t17 = -t190 * t360 + t217 + t219;
t14 = t232 / 0.2e1 + t216 + (m(4) + t363) * t291 / 0.2e1;
t8 = -t88 * t354 * t271 + t111 * t227 + t275 * t340 - t238 + t240;
t6 = t190 * (-t318 - t320) / 0.4e1 - t302 / 0.2e1 - t306 / 0.4e1 - t210 * (Ifges(7,2) * t295 - t179) / 0.4e1 - t312 / 0.4e1 + qJ(5) * t227 - t281 * t329 - t282 * t332 + Ifges(7,3) * t335 + t241 * t190 + t243 - (t196 / 0.2e1 - t255 / 0.4e1) * t294 - t355 * t339 + (t256 + t195) * t295 / 0.4e1;
t3 = t129 * t271 - mrSges(7,3) * t288 / 0.2e1 + t231 * t343 + mrSges(6,2) * t337 - t158 * mrSges(5,1) / 0.2e1 + (t242 * t343 + t224) * t111 - t218 + t350 + (-mrSges(5,2) / 0.2e1 + t359 / 0.2e1) * t159;
t9 = [t13 * qJD(2) + t12 * qJD(4), t299 + (t129 * t137 + t130 * t139 + t158 * t364) * qJD(2) + t14 * qJD(3) + t3 * qJD(4) + t28 * qJD(5) + t8 * qJD(6) - t351 * qJD(2) * t159 + ((-mrSges(3,2) + t265) * t214 + (-mrSges(4,1) * t209 - mrSges(5,1) * t354 + mrSges(4,2) * t207 + mrSges(5,2) * t190 - mrSges(3,1) + t144) * t212) * qJD(2) * t208 + 0.2e1 * ((t107 * t159 + t129 * t53 + t130 * t54) * t343 + (t131 * t291 + t246) * t344 + (t201 * t291 + t246) * t346 + (-pkin(2) * t212 + t214 * t264) * t208 * t347) * qJD(2), qJD(2) * t14, t300 + t3 * qJD(2) + t38 * qJD(5) + t21 * qJD(6) + ((-t111 * t263 - t304) * t343 + (-pkin(4) * t111 - t304) * t344) * t348 + ((mrSges(5,2) - t359) * t110 + (t325 - t365) * t111) * qJD(4), qJD(2) * t28 + qJD(4) * t38, t8 * qJD(2) + t21 * qJD(4) + (mrSges(7,1) * t89 - mrSges(7,2) * t88) * qJD(6); qJD(3) * t15 + qJD(4) * t2 + qJD(5) * t29 - qJD(6) * t7 - t299, qJD(3) * t10 + qJD(4) * t1 + qJD(5) * t18 + qJD(6) * t4, qJD(4) * t17 + qJD(5) * t50 + qJD(6) * t26 - t249, t17 * qJD(3) + t19 * qJD(5) + t6 * qJD(6) + ((-qJ(5) * t105 - t251 * t339) * t343 + (-pkin(4) * t147 - qJ(5) * t146) * t344) * t348 + t254 + (qJ(5) * t133 - t105 * t194 + (t95 / 0.2e1 - t339 * t138 - t57 * mrSges(7,3)) * t213 + (-t93 / 0.2e1 - t339 * t140 - t58 * mrSges(7,3)) * t210 - (pkin(4) * mrSges(6,1) + Ifges(7,5) * t330 + Ifges(7,6) * t331 + Ifges(6,4) - Ifges(5,5)) * t354 + (-qJ(5) * mrSges(6,1) + Ifges(6,5) - Ifges(5,6) + t233) * t190 + t325 * t147 + (-mrSges(6,3) + mrSges(5,2)) * t146) * qJD(4), qJD(3) * t50 + qJD(4) * t19 + qJD(6) * t23 + t248, t26 * qJD(3) + t6 * qJD(4) + t23 * qJD(5) + (-mrSges(7,1) * t54 - mrSges(7,2) * t53 + Ifges(7,6) * t295 - t178) * qJD(6) + t253; -qJD(2) * t15, -qJD(4) * t16 - qJD(5) * t49 - qJD(6) * t25 + t249, 0, -t296, -t278, t193 * qJD(6) - t279; -qJD(2) * t2 - qJD(5) * t40 - qJD(6) * t20 - t300, qJD(3) * t16 + qJD(5) * t24 - qJD(6) * t5 - t254, t296, qJD(5) * t182 + qJD(6) * t74, -t223 ((mrSges(7,2) * t339 - Ifges(7,6)) * t213 + (mrSges(7,1) * t339 - Ifges(7,5)) * t210) * qJD(6) - t225; -qJD(2) * t29 + qJD(4) * t40, qJD(3) * t49 - qJD(4) * t24 - qJD(6) * t22 - t248, t278, t223, 0, -t194 * qJD(6) - t280; t7 * qJD(2) + t20 * qJD(4), qJD(3) * t25 + qJD(4) * t5 + qJD(5) * t22 - t253, t279, t225, t280, 0;];
Cq  = t9;
