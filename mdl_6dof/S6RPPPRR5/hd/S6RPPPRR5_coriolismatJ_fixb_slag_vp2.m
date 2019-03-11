% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
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
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPPRR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR5_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR5_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:37:15
% EndTime: 2019-03-09 01:37:20
% DurationCPUTime: 2.64s
% Computational Cost: add. (4862->339), mult. (9270->521), div. (0->0), fcn. (7909->6), ass. (0->211)
t210 = sin(qJ(5));
t209 = sin(qJ(6));
t211 = cos(qJ(6));
t258 = mrSges(7,1) * t211 - mrSges(7,2) * t209;
t149 = t258 * t210;
t319 = Ifges(7,2) * t211;
t320 = Ifges(7,4) * t209;
t173 = t319 + t320;
t197 = Ifges(7,4) * t211;
t321 = Ifges(7,1) * t209;
t175 = t197 + t321;
t335 = t211 / 0.2e1;
t338 = -t209 / 0.2e1;
t364 = t173 * t338 + t175 * t335;
t201 = t209 ^ 2;
t203 = t211 ^ 2;
t363 = mrSges(7,3) * (t203 / 0.2e1 + t201 / 0.2e1);
t205 = sin(pkin(9));
t206 = cos(pkin(9));
t212 = cos(qJ(5));
t281 = t211 * t212;
t147 = t205 * t281 - t209 * t206;
t300 = t147 * t211;
t287 = t209 * t212;
t145 = -t205 * t287 - t211 * t206;
t302 = t145 * t209;
t229 = (t300 / 0.2e1 - t302 / 0.2e1) * mrSges(7,3);
t278 = t201 + t203;
t362 = t278 * mrSges(7,3);
t196 = Ifges(7,5) * t211;
t318 = Ifges(7,6) * t209;
t361 = Ifges(6,4) - t196 / 0.2e1 + t318 / 0.2e1;
t360 = -t209 * Ifges(7,2) + t197;
t314 = t211 * mrSges(7,2);
t317 = t209 * mrSges(7,1);
t172 = t314 + t317;
t242 = t314 / 0.2e1 + t317 / 0.2e1;
t336 = t210 / 0.2e1;
t359 = t172 * t336 + t242 * t210;
t358 = t210 * mrSges(6,1) + t212 * mrSges(6,2) - t149;
t207 = qJ(2) + pkin(3);
t208 = pkin(1) + qJ(3);
t159 = t205 * t207 - t206 * t208;
t153 = pkin(7) + t159;
t329 = t210 * pkin(5);
t330 = pkin(8) * t212;
t179 = t329 - t330;
t288 = t209 * t210;
t66 = t153 * t288 + t179 * t211;
t285 = t210 * t211;
t67 = -t153 * t285 + t179 * t209;
t357 = -t209 * t66 + t211 * t67;
t327 = qJD(4) / 0.2e1;
t148 = t209 * t205 + t206 * t281;
t299 = t148 * t211;
t146 = -t211 * t205 + t206 * t287;
t301 = t146 * t209;
t248 = t299 + t301;
t249 = -t300 + t302;
t291 = t206 * t212;
t34 = (t249 * t206 + (-t248 + 0.2e1 * t291) * t205) * t210;
t347 = t34 / 0.2e1;
t202 = t210 ^ 2;
t204 = t212 ^ 2;
t275 = -0.1e1 + t278;
t247 = -t275 * t202 - t204;
t40 = t206 * t247 + t212 * t248;
t200 = t206 ^ 2;
t47 = (t200 * t212 - t206 * t248) * t210;
t150 = t210 * t172;
t151 = t212 * t172;
t164 = -mrSges(7,1) * t212 - mrSges(7,3) * t285;
t274 = mrSges(7,3) * t288;
t162 = mrSges(7,2) * t212 - t274;
t283 = t211 * t162;
t227 = -t151 / 0.2e1 + t164 * t338 + t283 / 0.2e1;
t158 = t205 * t208 + t206 * t207;
t152 = -pkin(4) - t158;
t328 = t210 * pkin(8);
t331 = pkin(5) * t212;
t90 = t152 - t328 - t331;
t58 = t153 * t281 + t209 * t90;
t313 = t211 * t58;
t57 = -t153 * t287 + t211 * t90;
t253 = -t209 * t57 + t313;
t232 = 0.2e1 * t153 * t212 - t253;
t352 = -m(7) / 0.2e1;
t217 = t232 * t352 + t227;
t236 = t249 * pkin(8);
t165 = t210 * mrSges(7,1) - mrSges(7,3) * t281;
t163 = -t210 * mrSges(7,2) - mrSges(7,3) * t287;
t342 = -t163 / 0.2e1;
t239 = t146 * t165 / 0.2e1 + t148 * t342;
t350 = m(7) * pkin(5);
t252 = t258 / 0.2e1 + mrSges(6,1) + t350 / 0.2e1;
t254 = -t146 * t66 + t148 * t67;
t339 = -t206 / 0.2e1;
t351 = m(7) / 0.2e1;
t7 = (-t205 * mrSges(6,2) + t150 * t339) * t212 + t229 + t254 * t352 - t236 * t351 + (-t205 * t252 + t206 * t217) * t210 + t239;
t356 = -m(7) * (t47 * qJD(2) + qJD(3) * t347 + t40 * t327) + t7 * qJD(1);
t39 = t205 * t247 - t212 * t249;
t199 = t205 ^ 2;
t46 = (t199 * t212 + t205 * t249) * t210;
t228 = (-t299 / 0.2e1 - t301 / 0.2e1) * mrSges(7,3);
t235 = t248 * pkin(8);
t240 = -t145 * t165 / 0.2e1 + t147 * t342;
t255 = t145 * t66 + t147 * t67;
t6 = (t206 * mrSges(6,2) - t205 * t150 / 0.2e1) * t212 + t228 + t255 * t352 - t235 * t351 + (t205 * t217 + t206 * t252) * t210 + t240;
t355 = -m(7) * (qJD(2) * t347 + t46 * qJD(3) + t39 * t327) + t6 * qJD(1);
t277 = t202 + t204;
t354 = t277 * mrSges(6,3) + t210 * t150 - mrSges(5,2);
t332 = t212 / 0.2e1;
t353 = t164 * t288 / 0.2e1 + t150 * t332 + t151 * t336 + (t232 * t351 - t283 / 0.2e1) * t210;
t349 = -mrSges(7,1) / 0.2e1;
t348 = -mrSges(7,2) / 0.2e1;
t295 = t205 * t206;
t346 = m(6) * (-0.1e1 + t277) * t295;
t345 = t145 / 0.2e1;
t344 = t153 / 0.2e1;
t343 = -t162 / 0.2e1;
t341 = -t164 / 0.2e1;
t340 = t205 / 0.2e1;
t337 = t209 / 0.2e1;
t334 = -t212 / 0.2e1;
t333 = -t212 / 0.4e1;
t322 = m(7) * qJD(5);
t269 = t322 / 0.2e1;
t326 = t34 * t269;
t325 = t39 * t269;
t324 = t40 * t269;
t323 = m(7) * qJD(1);
t311 = t212 * Ifges(7,5);
t310 = t212 * Ifges(7,6);
t171 = -t212 * mrSges(6,1) + t210 * mrSges(6,2);
t307 = mrSges(5,1) - t171;
t306 = -mrSges(6,1) - t258;
t218 = -m(6) * (t277 * t199 + t200) / 0.2e1 + (t145 ^ 2 + t147 ^ 2 + t199 * t202) * t352;
t219 = m(6) * (-t277 * t200 - t199) / 0.2e1 + (-t146 ^ 2 - t148 ^ 2 - t200 * t202) * t351;
t23 = -m(4) + 0.2e1 * (-t200 / 0.2e1 - t199 / 0.2e1) * m(5) + t218 + t219;
t305 = qJD(1) * t23;
t216 = (t149 * t340 - t229) * t210 + t162 * t345 + t147 * t341;
t243 = t146 * t349 + t148 * t348;
t10 = t216 + t243;
t304 = t10 * qJD(1);
t215 = (t206 * t149 / 0.2e1 + t228) * t210 + t146 * t343 + t148 * t341;
t244 = mrSges(7,1) * t345 + t147 * t348;
t12 = t215 - t244;
t303 = t12 * qJD(1);
t298 = t153 * t210;
t241 = m(7) * t357;
t282 = t211 * t163;
t289 = t209 * t165;
t16 = m(7) * (t202 - t204) * t344 + (t253 * t351 + t227) * t212 + (t282 / 0.2e1 - t289 / 0.2e1 + t241 / 0.2e1 + t150 / 0.2e1) * t210;
t297 = t16 * qJD(1);
t296 = t202 * t153;
t294 = t205 * t210;
t293 = t205 * t212;
t292 = t206 * t210;
t126 = t210 * t360 - t310;
t290 = t209 * t126;
t176 = t211 * Ifges(7,1) - t320;
t128 = t210 * t176 - t311;
t284 = t211 * t128;
t237 = t162 * t337 + t164 * t335;
t25 = t149 * t332 + t202 * t363 + t237 * t210;
t280 = t25 * qJD(1);
t276 = qJD(5) * t212;
t273 = -pkin(5) * t149 / 0.2e1;
t272 = -t323 / 0.2e1;
t271 = t323 / 0.2e1;
t270 = qJD(4) * t352;
t262 = t196 - t318;
t256 = Ifges(7,5) * t209 + Ifges(7,6) * t211;
t127 = Ifges(7,6) * t210 + t212 * t360;
t129 = Ifges(7,5) * t210 + t176 * t212;
t1 = t67 * t162 + t58 * t163 + t66 * t164 + t57 * t165 + m(7) * (t57 * t66 + t58 * t67) + (t152 * mrSges(6,2) - t290 / 0.2e1 + t284 / 0.2e1 + t153 * t150 + t361 * t212) * t212 + (t152 * mrSges(6,1) + t127 * t338 + t129 * t335 + t153 * t151 - t361 * t210 + (m(7) * t153 ^ 2 + Ifges(6,1) - Ifges(6,2) - Ifges(7,3)) * t212) * t210;
t251 = t1 * qJD(1) + t16 * qJD(4);
t3 = t58 * t164 - t149 * t298 + (mrSges(7,3) * t313 + t126 * t335 + t128 * t337 + t364 * t210 + t256 * t334) * t210 + (-t162 - t274) * t57;
t250 = -t3 * qJD(1) - t25 * qJD(4);
t246 = t66 * t349 + t67 * mrSges(7,2) / 0.2e1;
t14 = t146 * t164 - t148 * t162 + mrSges(4,3) + t307 * t205 - t354 * t206 + m(7) * (t146 * t57 - t148 * t58 - t206 * t296) + m(6) * (-t277 * t206 * t153 - t152 * t205) + m(5) * (t158 * t205 - t159 * t206) + m(4) * t208;
t50 = (-t248 + t291) * t210;
t234 = -t14 * qJD(1) + t50 * t270;
t113 = t205 * t296;
t15 = t145 * t164 + t147 * t162 + mrSges(4,1) + mrSges(3,3) + t307 * t206 + (m(4) + m(3)) * qJ(2) + t354 * t205 + m(7) * (t145 * t57 + t147 * t58 + t113) + m(6) * (t153 * t204 * t205 - t152 * t206 + t113) + m(5) * (t158 * t206 + t159 * t205);
t49 = (-t249 - t293) * t210;
t233 = -t15 * qJD(1) + t49 * t270;
t230 = -t172 / 0.2e1 + t242;
t226 = -t145 * t146 + t147 * t148 + t202 * t295;
t87 = t275 * t212 * t210;
t223 = -t87 * qJD(4) - t40 * qJD(2) / 0.2e1 - t39 * qJD(3) / 0.2e1;
t221 = t210 * t230;
t28 = t205 * t221;
t29 = t206 * t221;
t213 = t172 * t344 - pkin(8) * t363 + (-t173 / 0.4e1 + t176 / 0.4e1 - t319 / 0.4e1) * t211 + (-t360 / 0.4e1 - t175 / 0.4e1 - t197 / 0.2e1 - t321 / 0.4e1) * t209;
t4 = t196 * t333 + t273 + (pkin(8) * t341 + t128 / 0.4e1 - t311 / 0.2e1) * t211 + (0.3e1 / 0.4e1 * t310 + pkin(8) * t343 - t126 / 0.4e1) * t209 + (-Ifges(7,3) / 0.2e1 + t213) * t210 + t246;
t43 = -pkin(5) * t172 + (t175 / 0.2e1 + t360 / 0.2e1) * t211 + (t176 / 0.2e1 - t173 / 0.2e1) * t209;
t77 = t230 * t212;
t220 = t4 * qJD(1) - t29 * qJD(2) - t28 * qJD(3) + t77 * qJD(4) + t43 * qJD(5);
t214 = qJD(5) * (-t212 * t258 + m(7) * (-t278 * t328 - t331) - t210 * t362 + t171);
t78 = t172 * t334 - t242 * t212;
t31 = t359 * t206;
t30 = t359 * t205;
t24 = -t218 + t219;
t13 = t215 + t244;
t11 = t216 - t243;
t9 = -mrSges(6,2) * t293 / 0.2e1 - mrSges(6,1) * t294 / 0.2e1 + t229 - t239 + (-pkin(5) * t294 - t236 + t254) * t351 + t358 * t340 + t353 * t206;
t8 = mrSges(6,1) * t292 / 0.2e1 + mrSges(6,2) * t291 / 0.2e1 + t228 - t240 + (pkin(5) * t292 - t235 + t255) * t351 + t358 * t339 + t353 * t205;
t5 = t262 * t333 + t284 / 0.4e1 - t290 / 0.4e1 + t273 - Ifges(7,6) * t287 / 0.2e1 + Ifges(7,5) * t281 / 0.2e1 + Ifges(7,3) * t336 - t237 * pkin(8) + t213 * t210 - t246;
t2 = t16 * qJD(5) - t25 * qJD(6) + 0.2e1 * (t49 * qJD(2) / 0.4e1 + t50 * qJD(3) / 0.4e1) * m(7);
t17 = [qJD(2) * t15 + qJD(3) * t14 + qJD(5) * t1 - qJD(6) * t3, 0.2e1 * (t346 / 0.2e1 + t226 * t351) * qJD(2) + t24 * qJD(3) + t9 * qJD(5) + t13 * qJD(6) - t233, t24 * qJD(2) + 0.2e1 * (t226 * t352 - t346 / 0.2e1) * qJD(3) + t8 * qJD(5) + t11 * qJD(6) - t234, t2, t9 * qJD(2) + t8 * qJD(3) + t5 * qJD(6) + (Ifges(6,5) + (t306 - t350) * t153 + t364) * t276 + t251 + (mrSges(6,2) * t298 - Ifges(6,6) * t210 - pkin(5) * t151 + t127 * t335 + t129 * t337 + t256 * t336 + (t241 + t282 - t289) * pkin(8) + t357 * mrSges(7,3)) * qJD(5), t13 * qJD(2) + t11 * qJD(3) + t5 * qJD(5) + (-mrSges(7,1) * t58 - mrSges(7,2) * t57 - t210 * t256) * qJD(6) + t250; t23 * qJD(3) - t7 * qJD(5) + t12 * qJD(6) + t233, t47 * t322, t305 + t326, t49 * t272 + t324, t31 * qJD(6) + t206 * t214 - t356, t303 + t31 * qJD(5) + (-mrSges(7,1) * t148 + mrSges(7,2) * t146) * qJD(6); -t23 * qJD(2) - t6 * qJD(5) + t10 * qJD(6) + t234, -t305 + t326, t46 * t322, t50 * t272 + t325, t30 * qJD(6) + t205 * t214 - t355, t304 + t30 * qJD(5) + (-mrSges(7,1) * t147 - mrSges(7,2) * t145) * qJD(6); t2, t49 * t271 + t324, t50 * t271 + t325, t87 * t322, t297 + t78 * qJD(6) + t306 * qJD(5) * t210 + (-mrSges(6,2) + t362) * t276 + ((t278 * t330 - t329) * qJD(5) - t223) * m(7), t78 * qJD(5) - t149 * qJD(6) - t280; qJD(2) * t7 + qJD(3) * t6 + qJD(6) * t4 - t251, -t29 * qJD(6) + t356, -t28 * qJD(6) + t355, m(7) * t223 + t77 * qJD(6) - t297, t43 * qJD(6) (-pkin(8) * t258 + t262) * qJD(6) + t220; -qJD(2) * t12 - qJD(3) * t10 - qJD(5) * t4 - t250, t29 * qJD(5) - t303, t28 * qJD(5) - t304, -t77 * qJD(5) + t280, -t220, 0;];
Cq  = t17;
