% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-03-08 19:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRPRPR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_coriolismatJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:24:54
% EndTime: 2019-03-08 19:25:00
% DurationCPUTime: 3.88s
% Computational Cost: add. (9187->336), mult. (21665->492), div. (0->0), fcn. (24655->12), ass. (0->207)
t212 = sin(qJ(4));
t215 = cos(qJ(4));
t298 = sin(pkin(12));
t300 = cos(pkin(12));
t183 = t212 * t298 - t215 * t300;
t214 = cos(qJ(6));
t308 = t214 * mrSges(7,2);
t211 = sin(qJ(6));
t313 = t211 * mrSges(7,1);
t188 = -t308 - t313;
t140 = t183 * t188;
t185 = -t212 * t300 - t215 * t298;
t291 = t185 * t211;
t146 = -t183 * mrSges(7,2) + mrSges(7,3) * t291;
t282 = t214 * t146;
t290 = t185 * t214;
t148 = t183 * mrSges(7,1) + mrSges(7,3) * t290;
t286 = t211 * t148;
t235 = t282 / 0.2e1 - t286 / 0.2e1;
t301 = cos(pkin(11));
t359 = pkin(2) * t301;
t201 = -pkin(3) - t359;
t186 = -t215 * pkin(4) + t201;
t129 = t183 * pkin(5) + t185 * pkin(9) + t186;
t209 = sin(pkin(11));
t360 = pkin(2) * t209;
t199 = pkin(8) + t360;
t258 = (-qJ(5) - t199) * t212;
t192 = t215 * t199;
t273 = t215 * qJ(5) + t192;
t351 = t298 * t258 + t273 * t300;
t71 = t129 * t214 - t211 * t351;
t72 = t129 * t211 + t214 * t351;
t246 = t211 * t71 - t214 * t72;
t348 = m(7) / 0.2e1;
t367 = t140 / 0.2e1 - t235 + (t351 + t246) * t348;
t213 = sin(qJ(2));
t299 = sin(pkin(6));
t254 = t301 * t299;
t331 = cos(qJ(2));
t255 = t331 * t299;
t169 = t209 * t255 + t213 * t254;
t210 = cos(pkin(6));
t154 = -t169 * t212 + t210 * t215;
t155 = t169 * t215 + t210 * t212;
t226 = t154 * t298 + t155 * t300;
t261 = t213 * t299;
t168 = t209 * t261 - t254 * t331;
t67 = t168 * t211 + t214 * t226;
t307 = t214 * t67;
t66 = t168 * t214 - t211 * t226;
t311 = t211 * t66;
t247 = -t307 + t311;
t366 = t247 + t226;
t365 = m(6) * t186;
t205 = t211 ^ 2;
t207 = t214 ^ 2;
t272 = t205 + t207;
t364 = (-0.1e1 + t272) * t185;
t136 = -t300 * t258 + t273 * t298;
t141 = t188 * t185;
t320 = t183 * mrSges(7,3);
t145 = mrSges(7,2) * t185 + t211 * t320;
t147 = -mrSges(7,1) * t185 + t214 * t320;
t292 = t168 * t212;
t317 = t185 * mrSges(6,3);
t343 = t226 / 0.2e1;
t345 = t66 / 0.2e1;
t349 = m(6) / 0.2e1;
t330 = t212 * pkin(4);
t144 = -pkin(5) * t185 + pkin(9) * t183 + t330;
t84 = t136 * t211 + t144 * t214;
t85 = -t136 * t214 + t144 * t211;
t95 = -t300 * t154 + t155 * t298;
t363 = (t343 - t226 / 0.2e1) * t317 + pkin(4) * t292 * t349 + (t136 * t226 + t66 * t84 + t67 * t85) * t348 + t147 * t345 + t67 * t145 / 0.2e1 + t141 * t343 + t367 * t95;
t187 = -mrSges(7,1) * t214 + mrSges(7,2) * t211;
t139 = t185 * t187;
t361 = -t139 / 0.2e1;
t265 = t298 * pkin(4);
t198 = t265 + pkin(9);
t358 = t198 * t272;
t357 = t141 - t317;
t203 = Ifges(7,4) * t214;
t356 = t211 * Ifges(7,1) + t203;
t245 = -t211 * t84 + t214 * t85;
t354 = t272 * mrSges(7,3);
t251 = Ifges(7,2) * t211 - t203;
t347 = m(6) * pkin(4);
t352 = t298 * t347 - mrSges(6,2);
t266 = t300 * pkin(4);
t200 = -t266 - pkin(5);
t350 = m(7) * t200 - t300 * t347 - mrSges(6,1) + t187;
t346 = -Ifges(7,3) / 0.2e1;
t341 = t168 / 0.2e1;
t340 = -t183 / 0.2e1;
t339 = t183 / 0.2e1;
t338 = t188 / 0.2e1;
t305 = t215 * mrSges(5,2);
t189 = t212 * mrSges(5,1) + t305;
t337 = t189 / 0.2e1;
t335 = t198 / 0.2e1;
t334 = t211 / 0.2e1;
t333 = -t214 / 0.2e1;
t332 = t214 / 0.2e1;
t328 = mrSges(7,3) * t185;
t327 = Ifges(7,4) * t211;
t326 = Ifges(7,5) * t214;
t324 = Ifges(7,6) * t211;
t107 = t185 * t168;
t323 = t107 * t95;
t321 = t183 * mrSges(6,3);
t319 = t183 * Ifges(7,5);
t318 = t183 * Ifges(7,6);
t316 = t185 * t95;
t309 = t212 * mrSges(5,2);
t108 = t183 * t168;
t87 = -t108 * t211 + t169 * t214;
t304 = t87 * t211;
t88 = t108 * t214 + t169 * t211;
t303 = t88 * t214;
t302 = t95 * t188;
t244 = -t303 + t304;
t295 = t107 * t183;
t28 = (t185 * t244 + t295) * t348 + (-t108 * t185 + t295) * t349;
t297 = qJD(1) * t28;
t296 = t107 * t136;
t293 = t136 * t185;
t119 = t185 * t251 + t318;
t289 = t211 * t119;
t288 = t211 * t146;
t287 = t211 * t147;
t190 = t214 * Ifges(7,2) + t327;
t285 = t211 * t190;
t253 = Ifges(7,1) * t214 - t327;
t121 = -t185 * t253 + t319;
t284 = t214 * t121;
t283 = t214 * t145;
t281 = t214 * t148;
t280 = t214 * t356;
t224 = (-t183 * t298 + t185 * t300) * t347;
t225 = m(7) * (-t183 * t358 - t185 * t200);
t218 = t225 / 0.2e1 + t361 + t224 / 0.2e1 - t272 * t320 / 0.2e1;
t271 = t347 / 0.2e1;
t221 = (t211 * t85 + t214 * t84) * t348 + t145 * t334 + t147 * t332 + t212 * t271;
t260 = t185 * mrSges(6,1) + t183 * mrSges(6,2);
t24 = -t218 + t221 - t260;
t279 = t24 * qJD(2);
t256 = mrSges(7,3) * (-t205 / 0.2e1 - t207 / 0.2e1);
t29 = t139 * t339 + (t281 / 0.2e1 + t288 / 0.2e1 + t185 * t256) * t185;
t278 = t29 * qJD(2);
t237 = t308 / 0.2e1 + t313 / 0.2e1;
t229 = t237 * t183;
t33 = t229 - t235;
t277 = t33 * qJD(2);
t50 = m(7) * t183 * t364;
t49 = t50 / 0.2e1;
t276 = t49 * qJD(2);
t11 = t367 * t183 + (-t141 / 0.2e1 + t287 / 0.2e1 - t283 / 0.2e1 + (-t136 - t245) * t348) * t185;
t270 = t11 * qJD(4) + t49 * qJD(5) + t29 * qJD(6);
t269 = -Ifges(7,2) / 0.4e1 + Ifges(7,1) / 0.4e1;
t250 = -t324 + t326;
t249 = Ifges(7,5) * t211 + Ifges(7,6) * t214;
t130 = t168 * t169;
t12 = m(7) * (t66 * t87 + t67 * t88 + t323) + m(5) * (t130 + (t154 * t212 - t155 * t215) * t168) + m(6) * (t108 * t226 + t130 + t323);
t243 = t12 * qJD(1) + t28 * qJD(3);
t13 = (t366 * t183 + t95 * t364) * t348;
t15 = m(7) * t366 * t95;
t242 = t15 * qJD(1) + t13 * qJD(3);
t241 = t84 * mrSges(7,1) / 0.2e1 - t85 * mrSges(7,2) / 0.2e1;
t240 = -t87 * mrSges(7,1) / 0.2e1 + t88 * mrSges(7,2) / 0.2e1;
t238 = t324 / 0.2e1 - t326 / 0.2e1;
t236 = t136 * t338 + t200 * t361;
t234 = t185 * t249;
t233 = t237 * t95;
t216 = -t244 * t198 * t348 + mrSges(5,1) * t292 / 0.2e1 + t305 * t341 + (-mrSges(6,2) / 0.2e1 + t298 * t271) * t108 + (t200 * t348 - mrSges(6,1) / 0.2e1 + t187 / 0.2e1 - t300 * t271) * t107 + (-t304 / 0.2e1 + t303 / 0.2e1) * mrSges(7,3);
t2 = t168 * t337 - t260 * t341 - t216 + t363;
t118 = -t185 * Ifges(7,6) + t183 * t251;
t120 = -t185 * Ifges(7,5) - t183 * t253;
t152 = mrSges(6,1) * t183 - mrSges(6,2) * t185;
t230 = t238 * t183;
t3 = t351 * t141 + t136 * t140 + t72 * t145 + t85 * t146 + t71 * t147 + t84 * t148 - t186 * t260 + t201 * t189 + (-Ifges(5,4) * t212 + pkin(4) * t152) * t212 + t330 * t365 + m(7) * (t351 * t136 + t71 * t84 + t72 * t85) + (-t284 / 0.2e1 + t289 / 0.2e1 + Ifges(6,4) * t183 + t230) * t183 + (Ifges(5,4) * t215 + (Ifges(5,1) - Ifges(5,2)) * t212) * t215 + (t120 * t333 + t118 * t334 + (-Ifges(6,4) - t238) * t185 + (-Ifges(7,3) + Ifges(6,1) - Ifges(6,2)) * t183) * t185;
t232 = t2 * qJD(1) + t3 * qJD(2) + t11 * qJD(3);
t220 = (-t311 / 0.2e1 + t307 / 0.2e1) * t328 + t146 * t345 - t67 * t148 / 0.2e1 + t95 * t139 / 0.2e1;
t6 = t220 + t240;
t7 = t136 * t139 + t71 * t146 - t72 * t148 + (t249 * t339 + t119 * t332 + t121 * t334 + (t285 / 0.2e1 + t356 * t333) * t185 - t246 * mrSges(7,3)) * t185;
t231 = t6 * qJD(1) + t7 * qJD(2) + t29 * qJD(3);
t228 = t13 * qJD(1) + t11 * qJD(2) + t50 * qJD(3);
t219 = (t247 * t183 - t316) * t348 + (-t183 * t226 - t316) * t349;
t223 = -m(6) * t169 / 0.2e1 - m(7) * (t211 * t88 + t214 * t87) / 0.2e1;
t18 = t219 + t223;
t21 = -t357 * t185 + (-t282 + t286 + t321) * t183 + m(7) * (t183 * t246 - t293) + m(6) * (-t183 * t351 - t293);
t227 = -qJD(1) * t18 - qJD(2) * t21 - qJD(3) * t49;
t19 = t302 / 0.2e1 + t233;
t8 = (0.3e1 / 0.4e1 * t318 + t119 / 0.4e1 + t146 * t335) * t211 + (t346 + t198 * t256 + (-t203 / 0.4e1 - t356 / 0.4e1 - t269 * t211) * t211) * t185 + (-0.3e1 / 0.4e1 * t319 - t121 / 0.4e1 + t148 * t335 + (-0.3e1 / 0.4e1 * t327 - t190 / 0.4e1 + t269 * t214) * t185) * t214 + t236 + t241;
t93 = -t200 * t188 + t280 / 0.2e1 + t253 * t334 - t285 / 0.2e1 - t251 * t332;
t99 = (t338 + t237) * t183;
t222 = t19 * qJD(1) + t8 * qJD(2) + t99 * qJD(3) - t93 * qJD(4);
t100 = t188 * t340 + t229;
t34 = t229 + t235;
t26 = t218 + t221;
t20 = -t302 / 0.2e1 + t233;
t17 = t219 - t223;
t9 = t183 * t250 / 0.4e1 + t284 / 0.4e1 - t289 / 0.4e1 + t185 * t346 + t230 - t236 + t241 - (t281 + t288) * t198 / 0.2e1 + (t190 / 0.2e1 - t253 / 0.4e1) * t290 + t272 * t328 * t335 + (0.2e1 * t356 - t251) * t291 / 0.4e1;
t5 = t220 - t240;
t4 = t28 * qJD(2) + t13 * qJD(4);
t1 = t216 + (t337 - t260 / 0.2e1) * t168 + t363;
t10 = [t12 * qJD(2) + t15 * qJD(4) (-mrSges(3,2) * t255 - mrSges(3,1) * t261 + m(7) * (t71 * t87 + t72 * t88 + t296) + t87 * t148 + t88 * t146 + m(6) * (t108 * t351 + t296) - t108 * t321 + (-m(4) * t359 + m(5) * t201 - mrSges(5,1) * t215 - mrSges(4,1) + t152 + t309 + t365) * t169 + t357 * t107 + (-m(4) * t360 + mrSges(4,2) + (m(5) * t199 + mrSges(5,3)) * (-t212 ^ 2 - t215 ^ 2)) * t168) * qJD(2) + t1 * qJD(4) + t17 * qJD(5) + t5 * qJD(6) + t243, t4, t1 * qJD(2) + t20 * qJD(6) + t242 + (-t155 * mrSges(5,1) - t154 * mrSges(5,2) + t350 * t226 - (m(7) * t358 + t352 + t354) * t95) * qJD(4), qJD(2) * t17, t5 * qJD(2) + t20 * qJD(4) + (-mrSges(7,1) * t67 - mrSges(7,2) * t66) * qJD(6); qJD(4) * t2 + qJD(5) * t18 + qJD(6) * t6 - t243, qJD(4) * t3 + qJD(5) * t21 + qJD(6) * t7, t270 - t297 (t266 * t321 + t265 * t317 + Ifges(5,5) * t215 - Ifges(5,6) * t212 - Ifges(6,5) * t183 + Ifges(6,6) * t185 - t234 / 0.2e1 + t120 * t334 + t118 * t332 + t200 * t140 + t285 * t339 + t280 * t340 + t199 * t309 - mrSges(5,1) * t192 - t352 * t136 + (m(7) * t245 + t283 - t287) * t198 + t350 * t351 + t245 * mrSges(7,3)) * qJD(4) + t26 * qJD(5) + t9 * qJD(6) + t232, qJD(4) * t26 + qJD(6) * t34 - t227, t9 * qJD(4) + t34 * qJD(5) + (-mrSges(7,1) * t72 - mrSges(7,2) * t71 + t234) * qJD(6) + t231; t4, t270 + t297, t50 * qJD(4) (-t183 * t354 - t139 - t189 + t224 + t225 + t260) * qJD(4) + t100 * qJD(6) + t228, t276, t100 * qJD(4) - t139 * qJD(6) + t278; -qJD(2) * t2 - qJD(6) * t19 - t242, -qJD(5) * t24 - qJD(6) * t8 - t232, -qJD(6) * t99 - t228, t93 * qJD(6), -t279 (t187 * t198 + t250) * qJD(6) - t222; -qJD(2) * t18, qJD(4) * t24 - qJD(6) * t33 + t227, -t276, t279, 0, t188 * qJD(6) - t277; -t6 * qJD(2) + t19 * qJD(4), qJD(4) * t8 + qJD(5) * t33 - t231, t99 * qJD(4) - t278, t222, t277, 0;];
Cq  = t10;
