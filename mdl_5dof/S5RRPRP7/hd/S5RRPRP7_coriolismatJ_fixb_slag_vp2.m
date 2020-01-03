% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRP7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP7_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP7_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP7_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP7_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:59:58
% EndTime: 2019-12-31 20:00:06
% DurationCPUTime: 3.65s
% Computational Cost: add. (6641->370), mult. (13363->490), div. (0->0), fcn. (13748->6), ass. (0->180)
t331 = sin(qJ(2));
t275 = t331 * pkin(2);
t367 = m(4) * t275;
t299 = sin(pkin(8));
t300 = cos(pkin(8));
t332 = cos(qJ(2));
t199 = t299 * t331 - t300 * t332;
t322 = mrSges(4,3) * t199;
t269 = t300 * pkin(2);
t224 = -t269 - pkin(3);
t234 = sin(qJ(4));
t235 = cos(qJ(4));
t256 = -t235 * pkin(4) - t234 * qJ(5);
t194 = t256 + t224;
t204 = -t235 * mrSges(6,1) - t234 * mrSges(6,3);
t366 = m(6) * t194 + t204;
t364 = mrSges(6,2) + mrSges(5,3);
t363 = -mrSges(5,1) - mrSges(6,1);
t362 = Ifges(6,4) + Ifges(5,5);
t361 = Ifges(6,2) + Ifges(5,3);
t360 = Ifges(5,6) - Ifges(6,6);
t294 = t199 * qJ(5);
t274 = t331 * pkin(6);
t206 = -qJ(3) * t331 - t274;
t276 = t332 * pkin(6);
t210 = qJ(3) * t332 + t276;
t356 = t299 * t206 + t300 * t210;
t282 = t235 * t356;
t201 = -t299 * t332 - t300 * t331;
t258 = -pkin(2) * t332 - pkin(1);
t139 = t199 * pkin(3) + t201 * pkin(7) + t258;
t283 = t234 * t139;
t67 = t282 + t283;
t51 = t67 + t294;
t324 = t51 - t67;
t66 = t139 * t235 - t234 * t356;
t52 = -t199 * pkin(4) - t66;
t323 = t52 + t66;
t307 = t235 * mrSges(5,2);
t311 = t234 * mrSges(5,1);
t209 = t307 + t311;
t359 = mrSges(4,3) + t209;
t290 = t201 * t234;
t273 = mrSges(5,3) * t290;
t142 = -t199 * mrSges(5,2) + t273;
t192 = t199 * mrSges(6,3);
t145 = mrSges(6,2) * t290 + t192;
t358 = t142 + t145;
t227 = Ifges(6,5) * t234;
t216 = Ifges(6,1) * t235 + t227;
t293 = t199 * t234;
t246 = mrSges(6,2) * t293 - mrSges(6,3) * t201;
t247 = mrSges(5,2) * t201 + mrSges(5,3) * t293;
t357 = t247 + t246;
t141 = -t201 * pkin(3) + t199 * pkin(7) + t275;
t158 = -t300 * t206 + t210 * t299;
t80 = t141 * t235 + t158 * t234;
t81 = t234 * t141 - t158 * t235;
t354 = -t234 * t80 + t235 * t81;
t53 = -qJ(5) * t201 + t81;
t54 = t201 * pkin(4) - t80;
t353 = t234 * t54 + t235 * t53;
t230 = Ifges(5,4) * t235;
t352 = Ifges(5,2) * t234 - t230;
t225 = m(6) * qJ(5) + mrSges(6,3);
t233 = t235 ^ 2;
t351 = m(5) / 0.2e1;
t350 = -m(6) / 0.2e1;
t349 = m(6) / 0.2e1;
t348 = m(4) * pkin(2);
t347 = -mrSges(5,1) / 0.2e1;
t346 = mrSges(6,1) / 0.2e1;
t345 = mrSges(5,2) / 0.2e1;
t344 = t54 / 0.2e1;
t343 = t67 / 0.2e1;
t342 = pkin(4) * mrSges(6,2);
t341 = qJ(5) / 0.2e1;
t289 = t201 * t235;
t187 = Ifges(6,5) * t289;
t340 = -t187 / 0.2e1;
t339 = -t194 / 0.2e1;
t338 = -t201 / 0.2e1;
t228 = Ifges(5,5) * t235;
t337 = t228 / 0.4e1;
t336 = -t234 / 0.2e1;
t335 = t234 / 0.2e1;
t334 = -t235 / 0.2e1;
t333 = t235 / 0.2e1;
t255 = t234 * pkin(4) - qJ(5) * t235;
t330 = m(6) * t255;
t328 = -Ifges(5,1) - Ifges(6,1);
t326 = -Ifges(5,2) - Ifges(6,3);
t320 = Ifges(5,4) * t234;
t229 = Ifges(6,4) * t235;
t319 = Ifges(6,5) * t235;
t317 = Ifges(5,6) * t234;
t226 = Ifges(6,6) * t234;
t316 = qJ(5) * mrSges(6,2);
t306 = t235 * mrSges(6,3);
t310 = t234 * mrSges(6,1);
t208 = -t306 + t310;
t138 = t208 * t201;
t143 = t199 * mrSges(5,1) + mrSges(5,3) * t289;
t144 = -t199 * mrSges(6,1) - mrSges(6,2) * t289;
t296 = t158 * t201;
t90 = -t201 * t255 + t158;
t5 = (t201 * t359 + t138) * t201 + (t322 - t358 * t235 + (t143 - t144) * t234) * t199 + m(6) * (-t90 * t201 + (-t234 * t52 - t235 * t51) * t199) + m(5) * (-t296 + (t234 * t66 - t235 * t67) * t199) + m(4) * (-t199 * t356 - t296);
t315 = qJD(1) * t5;
t244 = t208 * t199;
t245 = t209 * t199;
t292 = t199 * t235;
t272 = mrSges(6,2) * t292;
t313 = t201 * mrSges(6,1);
t248 = -t272 + t313;
t249 = -t201 * mrSges(5,1) + mrSges(5,3) * t292;
t250 = -Ifges(4,4) + t228 + t229;
t259 = -t201 * mrSges(4,1) - t199 * mrSges(4,2);
t261 = t360 * t199;
t312 = t201 * mrSges(4,3);
t89 = -t199 * t255 + t356;
t1 = pkin(1) * (mrSges(3,1) * t331 + mrSges(3,2) * t332) + (t331 ^ 2 - t332 ^ 2) * Ifges(3,4) - t67 * t247 - t52 * t248 - t66 * t249 - t51 * t246 + t90 * t244 - t356 * t312 - m(5) * (t158 * t356 + t66 * t80 + t67 * t81) - m(6) * (t51 * t53 + t52 * t54 + t89 * t90) + t89 * t138 - t81 * t142 - t80 * t143 - t54 * t144 - t53 * t145 + (-Ifges(3,1) + Ifges(3,2)) * t332 * t331 + ((t340 - t261) * t234 - mrSges(4,1) * t275 + t250 * t199) * t199 + ((t328 * t233 - Ifges(4,1) + Ifges(4,2) + (t326 * t234 + (0.2e1 * Ifges(5,4) - 0.3e1 / 0.2e1 * Ifges(6,5)) * t235) * t234 + t361) * t199 + mrSges(4,2) * t275 + (-t226 - t250 + t317) * t201 + t359 * t356) * t201 + (-t259 - t367) * t258 + t158 * t245;
t314 = t1 * qJD(1);
t137 = t256 * t201;
t205 = -t235 * mrSges(5,1) + t234 * mrSges(5,2);
t291 = t201 * t204;
t4 = -t205 * t296 - t67 * t144 - m(6) * (t137 * t90 + t52 * t67) + t137 * t138 - t90 * t291 + t67 * t143 + ((-t52 * mrSges(6,2) + (Ifges(6,5) - Ifges(5,4)) * t290 - t362 * t199) * t234 + (Ifges(5,4) * t289 - t51 * mrSges(6,2) - t67 * mrSges(5,3) - t187 - t261 + (t326 - t328) * t290) * t235) * t201 + (-m(6) * t51 + t273 - t358) * t66;
t303 = t4 * qJD(1);
t237 = (t234 * t81 + t235 * t80) * t351 + (t234 * t53 - t235 * t54) * t349 + t249 * t333 + t248 * t334 + t367 / 0.2e1 + t357 * t335;
t268 = t299 * pkin(2);
t223 = t268 + pkin(7);
t285 = t233 * t199;
t232 = t234 ^ 2;
t286 = t232 * t199;
t280 = (-t285 - t286) * t223;
t295 = t194 * t201;
t238 = (-t201 * t224 + t280) * t351 + (t280 - t295) * t349 - t291 / 0.2e1 + t205 * t338 + (-t199 * t299 + t201 * t300) * t348 / 0.2e1 + t364 * (-t286 / 0.2e1 - t285 / 0.2e1);
t6 = -t237 + t238 - t259;
t302 = t6 * qJD(1);
t239 = (t234 * t323 + t235 * t324) * t350 + t143 * t335 + t144 * t336 + t358 * t334;
t240 = (t307 / 0.2e1 + t311 / 0.2e1 + t310 / 0.2e1 - t306 / 0.2e1 + t255 * t349) * t199;
t8 = t240 + t239;
t301 = t8 * qJD(1);
t18 = t199 * t145 + m(6) * (t51 * t199 + t289 * t90) - t138 * t289;
t298 = qJD(1) * t18;
t281 = t235 * t204;
t279 = t229 + t226;
t134 = m(6) * t293;
t278 = t134 * qJD(1);
t277 = m(6) * t344;
t271 = t66 / 0.2e1 + t52 / 0.2e1;
t270 = t343 - t51 / 0.2e1;
t267 = -t293 / 0.2e1;
t266 = t293 / 0.2e1;
t265 = -t292 / 0.2e1;
t218 = Ifges(5,1) * t235 - t320;
t212 = Ifges(6,3) * t234 + t319;
t211 = -Ifges(6,3) * t235 + t227;
t213 = Ifges(5,2) * t235 + t320;
t215 = Ifges(6,1) * t234 - t319;
t217 = Ifges(5,1) * t234 + t230;
t251 = Ifges(5,2) / 0.4e1 + Ifges(6,3) / 0.4e1 - Ifges(5,1) / 0.4e1 - Ifges(6,1) / 0.4e1;
t236 = (-t218 / 0.4e1 - t216 / 0.4e1 + t213 / 0.4e1 - t211 / 0.4e1 + t224 * t347 + mrSges(6,1) * t339 + t251 * t235) * t235 + (t217 / 0.4e1 + t215 / 0.4e1 - t352 / 0.4e1 - t212 / 0.4e1 + t224 * t345 + mrSges(6,3) * t339 - t251 * t234 + (Ifges(5,4) - Ifges(6,5) / 0.2e1) * t235) * t234 + t364 * t223 * (t232 / 0.2e1 + t233 / 0.2e1);
t241 = (t194 * t137 + t255 * t90) * t349 + t137 * t204 / 0.2e1 + t158 * t209 / 0.2e1 - t255 * t138 / 0.2e1 + t90 * t208 / 0.2e1;
t242 = (-pkin(4) * t54 + qJ(5) * t53) * t350 - t53 * mrSges(6,3) / 0.2e1 + mrSges(6,1) * t344 + t80 * t347 + t81 * t345;
t2 = t187 * t336 + (t234 * t270 + t235 * t271) * mrSges(6,2) + ((-t234 * t324 + t323 * t235) * t349 + t144 * t333 + t143 * t334 + t358 * t336) * t223 + (t337 + t229 / 0.4e1 + t226 / 0.4e1 + (0.3e1 / 0.4e1 * Ifges(6,4) + 0.3e1 / 0.4e1 * Ifges(5,5) - t342 / 0.2e1) * t235 + (-Ifges(5,6) + 0.3e1 / 0.4e1 * Ifges(6,6) - t316 / 0.2e1) * t234) * t199 + (pkin(4) * t346 + mrSges(6,3) * t341 + Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + t236) * t201 + t241 + t242;
t20 = t255 * t204 + t224 * t209 + (t208 + t330) * t194 + (t215 / 0.2e1 - t212 / 0.2e1 + t217 / 0.2e1 - t352 / 0.2e1) * t235 + (t216 / 0.2e1 + t211 / 0.2e1 + t218 / 0.2e1 - t213 / 0.2e1) * t234;
t254 = -t2 * qJD(1) - t20 * qJD(2);
t243 = (-t234 * t90 + (t199 * t223 + t295) * t235) * t349 - t138 * t336;
t14 = -t272 + (t346 - t281 / 0.2e1) * t201 + t277 - t243;
t146 = t366 * t234;
t253 = qJD(1) * t14 + qJD(2) * t146;
t22 = -t192 + 0.2e1 * (t67 / 0.4e1 - t294 / 0.2e1 - t283 / 0.4e1 - t282 / 0.4e1) * m(6);
t252 = qJD(1) * t22 - qJD(4) * t225;
t191 = (m(6) * t223 + mrSges(6,2)) * t235;
t19 = (t67 + 0.2e1 * t294) * t349 + m(6) * t343 + t145;
t15 = t201 * t281 / 0.2e1 + t313 / 0.2e1 + t277 + t243;
t9 = t240 - t239;
t7 = t237 + t238;
t3 = -t242 + t241 + (t271 * mrSges(6,2) + (t323 * t349 + t144 / 0.2e1 - t143 / 0.2e1) * t223) * t235 + (t340 + t270 * mrSges(6,2) + (t324 * t350 - t145 / 0.2e1 - t142 / 0.2e1) * t223) * t234 + t236 * t201 + Ifges(5,6) * t266 + Ifges(6,6) * t267 + t246 * t341 - pkin(4) * t248 / 0.2e1 + t361 * t338 + t362 * t265 + ((Ifges(6,4) / 0.4e1 + Ifges(5,5) / 0.4e1) * t235 + (Ifges(6,6) / 0.4e1 - Ifges(5,6) / 0.2e1) * t234 + t337 + t279 / 0.4e1) * t199;
t10 = [-qJD(2) * t1 + qJD(3) * t5 - qJD(4) * t4 + qJD(5) * t18, t7 * qJD(3) + t3 * qJD(4) + t15 * qJD(5) - t314 + ((-Ifges(5,6) * t201 + t199 * t352) * t333 + t353 * mrSges(6,2) - Ifges(3,6) * t331 + Ifges(3,5) * t332 - (t299 * t348 - mrSges(4,2)) * t158 + (t362 * t234 + t360 * t235) * t338 + (-t362 * t201 + (-t218 - t216) * t199) * t335 + t354 * mrSges(5,3) + (t217 + t215) * t265 + t213 * t266 + t211 * t267 + (-Ifges(6,6) * t201 - t199 * t212) * t334 + t269 * t322 + t268 * t312 - mrSges(3,1) * t276 + mrSges(3,2) * t274 - t194 * t244 - t224 * t245 + (m(5) * t224 - t300 * t348 - mrSges(4,1) + t205) * t356 - Ifges(4,5) * t199 + Ifges(4,6) * t201 + (m(6) * t353 + m(5) * t354 + t357 * t235 + (t248 - t249) * t234) * t223 + t366 * t89) * qJD(2), qJD(2) * t7 + qJD(4) * t9 + t315, t3 * qJD(2) + t9 * qJD(3) + t19 * qJD(5) - t303 + ((-m(6) * pkin(4) + t363) * t67 + (-mrSges(5,2) + t225) * t66 + ((t316 + t360) * t235 + (t362 - t342) * t234) * t201) * qJD(4), qJD(2) * t15 + qJD(4) * t19 + t298; qJD(3) * t6 + qJD(4) * t2 - qJD(5) * t14 + t314, qJD(4) * t20 - qJD(5) * t146, t302, t191 * qJD(5) - t254 + (t228 + t279 - t317 + (m(6) * t256 + t204 + t205) * t223 + t256 * mrSges(6,2)) * qJD(4), qJD(4) * t191 - t253; -qJD(2) * t6 - qJD(4) * t8 + qJD(5) * t134 - t315, -t302, 0, -t301 + (t306 - t307 - t330) * qJD(4) + (m(6) * qJD(5) + t363 * qJD(4)) * t234, m(6) * qJD(4) * t234 + t278; -qJD(2) * t2 + qJD(3) * t8 - qJD(5) * t22 + t303, t254, t301, t225 * qJD(5), -t252; qJD(2) * t14 - qJD(3) * t134 + qJD(4) * t22 - t298, t253, -t278, t252, 0;];
Cq = t10;
