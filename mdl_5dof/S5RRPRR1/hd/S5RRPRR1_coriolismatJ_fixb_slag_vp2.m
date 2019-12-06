% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:26
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(4,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_coriolismatJ_fixb_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR1_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:23:56
% EndTime: 2019-12-05 18:24:05
% DurationCPUTime: 2.92s
% Computational Cost: add. (5242->298), mult. (11043->397), div. (0->0), fcn. (10577->6), ass. (0->191)
t319 = sin(qJ(2));
t187 = t319 * pkin(1);
t173 = t319 * pkin(2) + t187;
t360 = m(5) * t173;
t318 = sin(qJ(4));
t320 = cos(qJ(4));
t321 = cos(qJ(2));
t165 = t318 * t319 - t320 * t321;
t191 = sin(qJ(5));
t189 = t191 ^ 2;
t192 = cos(qJ(5));
t190 = t192 ^ 2;
t340 = t189 + t190;
t237 = t340 * t165;
t193 = pkin(2) + pkin(1);
t245 = t318 * t193;
t175 = t245 + pkin(4);
t166 = -t318 * t321 - t320 * t319;
t248 = t320 * t166;
t169 = -mrSges(6,1) * t192 + mrSges(6,2) * t191;
t282 = t166 * t169;
t253 = -t282 / 0.2e1 - mrSges(6,3) * t237 / 0.2e1;
t317 = m(5) * t193;
t334 = m(6) / 0.2e1;
t359 = (-t175 * t237 + t193 * t248) * t334 + (-t318 * t165 + t248) * t317 / 0.2e1 + t253;
t228 = Ifges(6,5) * t191 + Ifges(6,6) * t192;
t216 = t166 * t228;
t322 = t192 / 0.2e1;
t183 = Ifges(6,4) * t192;
t172 = Ifges(6,1) * t191 + t183;
t271 = t192 * t172;
t313 = Ifges(6,4) * t191;
t171 = Ifges(6,2) * t192 + t313;
t275 = t191 * t171;
t344 = -t271 / 0.2e1 + t275 / 0.2e1;
t349 = t191 / 0.2e1;
t230 = Ifges(6,2) * t191 - t183;
t75 = -Ifges(6,6) * t166 + t230 * t165;
t231 = Ifges(6,1) * t192 - t313;
t77 = -Ifges(6,5) * t166 - t231 * t165;
t214 = -t216 / 0.2e1 + Ifges(5,6) * t166 + t77 * t349 + t75 * t322 + (-Ifges(5,5) + t344) * t165;
t246 = t319 * qJ(3);
t213 = -t319 * pkin(3) - t246;
t179 = t321 * qJ(3);
t269 = t321 * pkin(3) + t179;
t335 = t318 * t213 + t320 * t269;
t345 = t335 * t169;
t346 = t335 * mrSges(5,1);
t130 = -t320 * t213 + t318 * t269;
t356 = t130 * mrSges(5,2);
t358 = t214 + t345 - t346 + t356;
t357 = -t345 / 0.2e1 + t346 / 0.2e1 - t356 / 0.2e1;
t355 = t130 * t191;
t354 = t130 * t192;
t182 = Ifges(6,5) * t192;
t311 = Ifges(6,6) * t191;
t229 = t311 - t182;
t353 = t229 * t165;
t352 = t318 * t130;
t289 = t335 * t130;
t325 = -t191 / 0.2e1;
t348 = -Ifges(6,3) * t166 / 0.2e1;
t347 = -m(6) * pkin(4) / 0.2e1;
t304 = t335 * mrSges(5,3);
t249 = t320 * t335;
t284 = t165 * t191;
t115 = mrSges(6,2) * t166 + mrSges(6,3) * t284;
t283 = t165 * t192;
t117 = -mrSges(6,1) * t166 + mrSges(6,3) * t283;
t270 = t115 * t349 + t117 * t322;
t81 = pkin(4) * t283 + t355;
t82 = pkin(4) * t284 - t354;
t343 = t270 + (t191 * t82 + t192 * t81) * t334;
t299 = t192 * mrSges(6,2);
t302 = t191 * mrSges(6,1);
t170 = t299 + t302;
t113 = t170 * t165;
t114 = t170 * t166;
t262 = t321 * pkin(1);
t174 = -t321 * pkin(2) - t262;
t240 = t166 * mrSges(5,1) + t165 * mrSges(5,2);
t139 = t166 * pkin(4) + t174;
t60 = t139 * t192 - t191 * t335;
t62 = t139 * t191 + t192 * t335;
t342 = t130 * t113 + t335 * t114 - t62 * t115 - t60 * t117 - t166 * t304 + t174 * t240;
t341 = -Ifges(5,1) + Ifges(6,3);
t339 = 0.3e1 / 0.4e1 * t313 + t171 / 0.4e1;
t336 = t230 * t322 + t231 * t325;
t333 = m(4) * pkin(1);
t331 = -mrSges(6,1) / 0.2e1;
t330 = mrSges(6,2) / 0.2e1;
t329 = Ifges(6,3) / 0.2e1;
t328 = -t165 / 0.2e1;
t323 = -t192 / 0.2e1;
t316 = m(6) * t193;
t281 = t166 * t191;
t261 = mrSges(6,3) * t281;
t116 = -t165 * mrSges(6,2) + t261;
t280 = t166 * t192;
t118 = t165 * mrSges(6,1) + mrSges(6,3) * t280;
t157 = Ifges(5,4) * t165;
t243 = -t283 / 0.2e1;
t244 = t284 / 0.2e1;
t267 = t319 ^ 2;
t268 = t321 ^ 2;
t138 = pkin(4) * t165 + t173;
t59 = t138 * t192 + t355;
t61 = t138 * t191 - t354;
t76 = Ifges(6,6) * t165 + t230 * t166;
t78 = Ifges(6,5) * t165 - t231 * t166;
t1 = -t77 * t280 / 0.2e1 + t75 * t281 / 0.2e1 + m(6) * (t59 * t60 + t61 * t62 + t289) + t174 * t360 + t76 * t244 + t78 * t243 + t61 * t116 + t59 * t118 - 0.2e1 * t157 * t328 + (Ifges(4,1) + Ifges(3,1) + (-(2 * mrSges(4,1)) - t333) * pkin(1) - Ifges(4,2) - Ifges(3,2)) * t321 * t319 + (t348 + t173 * mrSges(5,1) + t353 / 0.2e1) * t165 + (-t304 - t173 * mrSges(5,2) + (-Ifges(5,1) + Ifges(5,2)) * t328 + (-Ifges(5,2) / 0.2e1 - t341 / 0.2e1) * t165 + (-Ifges(5,4) - t229 / 0.2e1) * t166) * t166 + (pkin(1) * mrSges(4,2) - Ifges(3,4) - Ifges(4,4)) * (-t268 + t267) - t342;
t309 = t1 * qJD(1);
t301 = t191 * t59;
t300 = t191 * t76;
t298 = t192 * t61;
t297 = t192 * t62;
t296 = t192 * t78;
t256 = t182 / 0.2e1;
t220 = t256 - t311 / 0.2e1;
t2 = -t82 * t116 - t81 * t118 - m(6) * (t60 * t81 + t62 * t82 + t289) + (-t300 / 0.2e1 + t296 / 0.2e1 - t157 + t220 * t165) * t165 + (t75 * t325 + t77 * t322 + t304 + (Ifges(5,4) - t220) * t166 + (Ifges(5,2) + t341) * t165) * t166 + t342;
t295 = t2 * qJD(1);
t158 = -t275 / 0.2e1;
t7 = t62 * t118 - t130 * t282 + (t76 * t323 + t78 * t325 + t228 * t328 - mrSges(6,3) * t297 + (t172 * t322 + t158) * t166) * t166 + (-t116 + t261) * t60;
t294 = t7 * qJD(1);
t293 = t81 * t191;
t292 = t82 * t192;
t273 = t192 * t116;
t276 = t191 * t118;
t10 = (t273 - t276 - m(6) * (t191 * t60 - t297) + m(5) * t335 - t165 * mrSges(5,3)) * t165 + (-m(4) * qJ(3) - mrSges(4,3)) * (t268 + t267) + ((m(6) + m(5)) * t130 - t114 - mrSges(5,3) * t166) * t166;
t291 = qJD(1) * t10;
t203 = t360 / 0.2e1 + (t191 * t61 + t192 * t59) * t334 + t270;
t11 = m(4) * t187 + t319 * mrSges(4,1) + t321 * mrSges(4,2) + t203 - t240 - t359;
t290 = t11 * qJD(1);
t206 = t237 * t347 + t253;
t14 = t206 + t240 - t343;
t285 = t14 * qJD(1);
t210 = (t299 / 0.2e1 + t302 / 0.2e1) * t165;
t217 = t276 / 0.2e1 - t273 / 0.2e1;
t18 = t210 + t217;
t279 = t18 * qJD(1);
t277 = t191 * t117;
t274 = t192 * t115;
t266 = mrSges(6,3) * t301;
t265 = mrSges(6,3) * t298;
t264 = mrSges(6,3) * t293;
t263 = mrSges(6,3) * t292;
t260 = -Ifges(6,2) / 0.4e1 + Ifges(6,1) / 0.4e1;
t259 = t175 * t277;
t258 = t175 * t274;
t252 = t191 * t320;
t251 = t192 * t320;
t250 = t193 * t320;
t247 = t320 * t170;
t63 = t271 / 0.2e1 + t158 - t336;
t234 = -t250 / 0.2e1;
t232 = mrSges(6,3) * (t190 / 0.2e1 + t189 / 0.2e1);
t227 = t298 - t301;
t226 = t292 - t293;
t198 = (t169 - mrSges(5,1)) * t245 + (t340 * mrSges(6,3) - mrSges(5,2)) * t250;
t212 = t340 * t320;
t28 = (t212 * t175 - t320 * t245) * t316 + t198;
t194 = (t226 * t175 + (t62 * t251 - t60 * t252 - t249 + t352) * t193) * t334 - t259 / 0.2e1 + t258 / 0.2e1 - t264 / 0.2e1 + t263 / 0.2e1 - t114 * t245 / 0.2e1 + t250 * t273 / 0.2e1 + (-t113 + t276) * t234 - t357;
t196 = t227 * t347 + t266 / 0.2e1 - t265 / 0.2e1 + (t277 / 0.2e1 - t274 / 0.2e1) * pkin(4) + t357;
t4 = t194 + t196;
t224 = t4 * qJD(1) + t28 * qJD(2);
t200 = t336 + t344;
t40 = t193 * t247 + t200;
t195 = t169 * t234 + Ifges(6,2) * t190 / 0.4e1 + (-t190 / 0.4e1 + t189 / 0.4e1) * Ifges(6,1) + t175 * t232 + (-t230 + t172) * t191 / 0.4e1 + t339 * t192;
t207 = t130 * t170 / 0.2e1 - t300 / 0.4e1 + t296 / 0.4e1;
t201 = (-0.3e1 / 0.4e1 * t311 + t182 / 0.4e1 + t256) * t165 + t207;
t218 = t116 * t325 + t118 * t323;
t208 = t218 * t175;
t222 = t61 * t330 + t59 * t331;
t5 = t208 + (t329 + t195) * t166 + t201 + t222;
t223 = t5 * qJD(1) - t40 * qJD(2);
t221 = t82 * t330 + t81 * t331;
t211 = t218 * pkin(4);
t204 = -mrSges(6,2) * t251 / 0.2e1 + t252 * t331;
t29 = (t247 / 0.2e1 + t204) * t193 + t200;
t197 = pkin(4) * t232 + (t183 / 0.4e1 + t172 / 0.4e1 + t260 * t191) * t191 + (-t260 * t192 + t339) * t192;
t8 = t211 + (t329 + t197) * t166 + t201 + t221;
t209 = t8 * qJD(1) - t29 * qJD(2) + t63 * qJD(4);
t202 = Ifges(6,6) * t244 + Ifges(6,5) * t243 + t348 - t353 / 0.4e1 + t207;
t30 = t170 * t234 + t204 * t193 + t63;
t19 = t210 - t217;
t15 = t206 + t343;
t12 = t203 + t359;
t9 = t197 * t166 + t202 + t211 - t221;
t6 = t195 * t166 + t202 + t208 - t222;
t3 = -t196 + t194 + t214;
t13 = [qJD(2) * t1 - qJD(3) * t10 - qJD(4) * t2 - qJD(5) * t7, t309 + (m(6) * (t227 * t175 - t193 * t249) + mrSges(4,2) * t246 - mrSges(4,3) * t262 + (-t249 - t352) * t317 - t259 - t266 + t258 + t265 + t113 * t250 + (Ifges(4,5) + Ifges(3,5)) * t321 + (-Ifges(4,6) - Ifges(3,6)) * t319 + (-mrSges(4,1) - t333) * t179 + (t165 * t250 + t166 * t245) * mrSges(5,3) + t358) * qJD(2) + t12 * qJD(3) + t3 * qJD(4) + t6 * qJD(5), qJD(2) * t12 + qJD(4) * t15 + qJD(5) * t19 - t291, t3 * qJD(2) + t15 * qJD(3) + t9 * qJD(5) - t295 + (t263 - t264 + (m(6) * t226 + t274 - t277) * pkin(4) + t358) * qJD(4), -t294 + t6 * qJD(2) + t19 * qJD(3) + t9 * qJD(4) + (-mrSges(6,1) * t62 - mrSges(6,2) * t60 + t216) * qJD(5); -qJD(3) * t11 + qJD(4) * t4 + qJD(5) * t5 - t309, qJD(4) * t28 - qJD(5) * t40, -t290, (t212 * pkin(4) * t316 + t198) * qJD(4) + t30 * qJD(5) + t224, t30 * qJD(4) + (t169 * t175 - t229) * qJD(5) + t223; qJD(2) * t11 - qJD(4) * t14 - qJD(5) * t18 + t291, t290, 0, -t285, -qJD(5) * t170 - t279; -qJD(2) * t4 + qJD(3) * t14 + qJD(5) * t8 + t295, -qJD(5) * t29 - t224, t285, t63 * qJD(5), (t169 * pkin(4) - t229) * qJD(5) + t209; -qJD(2) * t5 + qJD(3) * t18 - qJD(4) * t8 + t294, qJD(4) * t29 - t223, t279, -t209, 0;];
Cq = t13;
