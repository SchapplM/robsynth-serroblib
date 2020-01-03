% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRR6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR6_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR6_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR6_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR6_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:01:00
% EndTime: 2019-12-31 19:01:07
% DurationCPUTime: 2.98s
% Computational Cost: add. (7643->302), mult. (15042->412), div. (0->0), fcn. (15026->8), ass. (0->197)
t208 = sin(qJ(5));
t205 = t208 ^ 2;
t209 = cos(qJ(5));
t206 = t209 ^ 2;
t353 = t205 + t206;
t329 = sin(qJ(4));
t330 = sin(qJ(3));
t331 = cos(qJ(4));
t332 = cos(qJ(3));
t183 = t329 * t330 - t331 * t332;
t184 = -t329 * t332 - t331 * t330;
t318 = Ifges(6,2) * t209;
t320 = Ifges(6,4) * t208;
t189 = t318 + t320;
t202 = Ifges(6,4) * t209;
t322 = Ifges(6,1) * t208;
t191 = t202 + t322;
t248 = Ifges(6,5) * t208 + Ifges(6,6) * t209;
t234 = t184 * t248;
t292 = t183 * t209;
t264 = -t292 / 0.2e1;
t293 = t183 * t208;
t265 = t293 / 0.2e1;
t333 = t209 / 0.2e1;
t335 = t208 / 0.2e1;
t250 = Ifges(6,2) * t208 - t202;
t96 = -Ifges(6,6) * t184 + t250 * t183;
t192 = Ifges(6,1) * t209 - t320;
t98 = -Ifges(6,5) * t184 - t192 * t183;
t231 = t189 * t265 + t191 * t264 + t96 * t333 + t98 * t335 + Ifges(5,6) * t184 - Ifges(5,5) * t183 - t234 / 0.2e1;
t251 = mrSges(6,1) * t209 - t208 * mrSges(6,2);
t254 = sin(pkin(9)) * pkin(1) + pkin(6);
t230 = t330 * t254;
t221 = -t330 * pkin(7) - t230;
t193 = t332 * t254;
t282 = t332 * pkin(7) + t193;
t350 = t329 * t221 + t331 * t282;
t360 = t350 * t251;
t363 = t350 * mrSges(5,1);
t132 = -t331 * t221 + t329 * t282;
t371 = t132 * mrSges(5,2);
t373 = t231 - t360 - t363 + t371;
t372 = t360 / 0.2e1 + t363 / 0.2e1 - t371 / 0.2e1;
t370 = t132 * t209;
t297 = t132 * t350;
t369 = t208 * t132;
t368 = t329 * t132;
t307 = t209 * mrSges(6,2);
t310 = t208 * mrSges(6,1);
t188 = t307 + t310;
t136 = t183 * t188;
t291 = t184 * t208;
t274 = mrSges(6,3) * t291;
t139 = -mrSges(6,2) * t183 + t274;
t290 = t184 * t209;
t141 = t183 * mrSges(6,1) + mrSges(6,3) * t290;
t269 = -cos(pkin(9)) * pkin(1) - pkin(2);
t186 = -t332 * pkin(3) + t269;
t124 = t183 * pkin(4) + t184 * pkin(8) + t186;
t47 = t208 * t124 + t209 * t350;
t306 = t209 * t47;
t334 = -t209 / 0.2e1;
t347 = m(6) / 0.2e1;
t46 = t124 * t209 - t208 * t350;
t365 = ((t208 * t46 - t306 + t350) * t347 + t139 * t334 + t141 * t335 - t136 / 0.2e1) * t183;
t364 = pkin(4) * t350;
t362 = t350 * mrSges(5,3);
t277 = t331 * pkin(3);
t200 = -t277 - pkin(4);
t361 = t200 * t350;
t135 = t251 * t184;
t336 = -t208 / 0.2e1;
t359 = t189 * t336 + t191 * t333;
t358 = t353 * t183;
t138 = t184 * mrSges(6,2) + mrSges(6,3) * t293;
t140 = -t184 * mrSges(6,1) + mrSges(6,3) * t292;
t233 = t188 * t184;
t261 = -t184 * mrSges(5,1) - t183 * mrSges(5,2);
t311 = t184 * mrSges(5,3);
t356 = t132 * t136 - t47 * t138 - t46 * t140 - t186 * t261 + (t233 - t311) * t350;
t355 = -Ifges(5,1) + Ifges(6,3);
t237 = t307 / 0.2e1 + t310 / 0.2e1;
t76 = (-t188 / 0.2e1 + t237) * t183;
t284 = t76 * qJD(2);
t257 = -mrSges(6,2) * t331 / 0.2e1;
t337 = -t200 / 0.2e1;
t344 = -mrSges(6,1) / 0.2e1;
t346 = pkin(4) / 0.2e1;
t38 = (t337 + t346) * t188 + (pkin(3) * t257 - t191 / 0.2e1 + t250 / 0.2e1) * t209 + (t277 * t344 - t192 / 0.2e1 + t189 / 0.2e1) * t208;
t258 = t192 * t335 - t250 * t333 + t359;
t327 = pkin(4) * t188;
t69 = t258 - t327;
t218 = (t189 / 0.4e1 - t192 / 0.4e1 + t318 / 0.4e1) * t209 + (-t250 / 0.4e1 + t191 / 0.4e1 + t202 / 0.2e1 + t322 / 0.4e1) * t208;
t252 = mrSges(6,3) * (t206 / 0.2e1 + t205 / 0.2e1);
t214 = pkin(8) * t252 + t218;
t99 = Ifges(6,5) * t183 - t184 * t192;
t304 = t209 * t99;
t97 = Ifges(6,6) * t183 + t250 * t184;
t308 = t208 * t97;
t338 = t188 / 0.2e1;
t222 = t132 * t338 - t308 / 0.4e1 + t304 / 0.4e1;
t235 = t139 * t336 + t141 * t334;
t217 = t235 * pkin(8) + t135 * t346 + t222;
t201 = Ifges(6,5) * t209;
t271 = t201 / 0.2e1;
t317 = Ifges(6,6) * t208;
t223 = (-0.3e1 / 0.4e1 * t317 + t201 / 0.4e1 + t271) * t183;
t343 = mrSges(6,2) / 0.2e1;
t325 = t184 * pkin(4);
t326 = t183 * pkin(8);
t142 = -t325 + t326;
t57 = t142 * t209 + t369;
t58 = t208 * t142 - t370;
t239 = t58 * t343 + t57 * t344;
t342 = Ifges(6,3) / 0.2e1;
t9 = t223 + (t342 + t214) * t184 + t217 + t239;
t352 = -t9 * qJD(1) + t38 * qJD(3) - t69 * qJD(4) + t284;
t288 = t200 * t188;
t62 = t258 + t288;
t275 = t329 * pkin(3);
t199 = t275 + pkin(8);
t213 = t199 * t252 + t218;
t216 = t135 * t337 + t235 * t199 + t222;
t276 = t330 * pkin(3);
t137 = t276 + t142;
t53 = t137 * t209 + t369;
t54 = t208 * t137 - t370;
t240 = t54 * t343 + t53 * t344;
t7 = t223 + (t342 + t213) * t184 + t216 + t240;
t351 = t7 * qJD(1) + t62 * qJD(3) - t284;
t249 = t317 - t201;
t316 = Ifges(6,3) * t184;
t349 = Ifges(6,5) * t264 + Ifges(6,6) * t265 - t249 * t183 / 0.4e1 - t316 / 0.2e1;
t345 = m(6) * pkin(3);
t341 = -t183 / 0.2e1;
t340 = t183 / 0.2e1;
t328 = pkin(4) * t136;
t30 = m(6) * (-0.1e1 + t353) * t184 * t183;
t285 = t30 * qJD(2);
t77 = (t237 + t338) * t183;
t324 = t77 * qJD(5) + t285;
t323 = -t76 * qJD(5) - t285;
t309 = t208 * t57;
t305 = t209 * t58;
t303 = t53 * t208;
t302 = t54 * t209;
t17 = t135 * t340 + (t184 * t252 + t235) * t184;
t301 = qJD(1) * t17;
t286 = t209 * t138;
t262 = -t286 / 0.2e1;
t287 = t208 * t140;
t263 = t287 / 0.2e1;
t219 = t237 * t184 + t262 + t263;
t246 = t302 - t303;
t11 = t365 + (t219 + (-t132 - t246) * t347) * t184;
t300 = t11 * qJD(1);
t245 = t305 - t309;
t12 = ((-t132 - t245) * t347 + t219) * t184 + t365;
t299 = t12 * qJD(1);
t289 = t200 * t136;
t280 = qJD(3) + qJD(4);
t279 = mrSges(6,3) * t303;
t278 = mrSges(6,3) * t302;
t273 = t199 * t287;
t272 = t199 * t286;
t270 = -t200 * t184 - t358 * t199;
t268 = t208 * t331;
t267 = t209 * t331;
t266 = t329 * t183;
t255 = -t268 / 0.2e1;
t178 = Ifges(5,4) * t183;
t229 = t330 * mrSges(4,1) + t332 * mrSges(4,2);
t2 = -t316 * t340 + t269 * t229 - t98 * t290 / 0.2e1 + t99 * t264 + t96 * t291 / 0.2e1 + t97 * t265 + t54 * t139 + t53 * t141 + m(6) * (t46 * t53 + t47 * t54 + t297) + m(5) * t186 * t276 - 0.2e1 * t178 * t341 + (Ifges(4,1) - Ifges(4,2)) * t332 * t330 + (-t330 ^ 2 + t332 ^ 2) * Ifges(4,4) + (mrSges(5,1) * t276 + t249 * t340) * t183 + (-t362 - mrSges(5,2) * t276 + (-Ifges(5,1) + Ifges(5,2)) * t341 + (-Ifges(5,2) / 0.2e1 - t355 / 0.2e1) * t183 + (-Ifges(5,4) - t249 / 0.2e1) * t184) * t184 - t356;
t244 = t2 * qJD(1) + t11 * qJD(2);
t238 = t271 - t317 / 0.2e1;
t4 = -m(6) * (t46 * t57 + t47 * t58 + t297) - t58 * t139 - t57 * t141 + (-t308 / 0.2e1 + t304 / 0.2e1 - t178 + t238 * t183) * t183 + (t362 + t96 * t336 + t98 * t333 + (Ifges(5,4) - t238) * t184 + (Ifges(5,2) + t355) * t183) * t184 + t356;
t243 = -t4 * qJD(1) + t12 * qJD(2);
t6 = t47 * t141 + t132 * t135 + (-mrSges(6,3) * t306 + t359 * t184 + t248 * t341 + t97 * t334 + t99 * t336) * t184 + (-t139 + t274) * t46;
t242 = -t6 * qJD(1) - t17 * qJD(2);
t232 = -t358 * mrSges(6,3) + t135 - t261;
t228 = t353 * t331;
t212 = ((-t228 * t184 + t266) * pkin(3) + t270) * t347;
t224 = m(6) * (-t326 * t353 + t325);
t20 = -t224 / 0.2e1 + t212;
t215 = (-mrSges(5,1) - t251) * t275 + (t353 * mrSges(6,3) - mrSges(5,2)) * t277;
t48 = (t228 * t199 + t329 * t200) * t345 + t215;
t210 = (t245 * t199 + t361) * t347 - t289 / 0.2e1 - t273 / 0.2e1 + t272 / 0.2e1 - t233 * t275 / 0.2e1 + ((t47 * t267 - t46 * t268 + t368) * t347 + t141 * t255 + t139 * t267 / 0.2e1) * pkin(3) + (-t309 / 0.2e1 + t305 / 0.2e1) * mrSges(6,3) - t372;
t211 = -m(6) * (t246 * pkin(8) - t364) / 0.2e1 - t328 / 0.2e1 + pkin(8) * t263 + pkin(8) * t262 + t279 / 0.2e1 - t278 / 0.2e1 + t372;
t5 = t210 + t211;
t227 = t5 * qJD(1) + t20 * qJD(2) + t48 * qJD(3);
t39 = t288 / 0.2e1 - t327 / 0.2e1 + (mrSges(6,1) * t255 + t209 * t257) * pkin(3) + t258;
t19 = t224 / 0.2e1 + t212 + t232;
t10 = t214 * t184 + t217 - t239 + t349;
t8 = t213 * t184 + t216 - t240 + t349;
t3 = t210 - t211 + t231;
t1 = qJD(3) * t11 + qJD(4) * t12 - qJD(5) * t17;
t13 = [qJD(3) * t2 - qJD(4) * t4 - qJD(5) * t6, t1, (mrSges(4,2) * t230 + t275 * t311 + t183 * mrSges(5,3) * t277 + Ifges(4,5) * t332 - Ifges(4,6) * t330 - mrSges(4,1) * t193 + m(5) * (-t331 * t350 - t368) * pkin(3) + m(6) * (t199 * t246 + t361) + t278 - t279 + t272 - t273 - t289 + t373) * qJD(3) + t3 * qJD(4) + t8 * qJD(5) + t244, t3 * qJD(3) + (-m(6) * t364 + t245 * mrSges(6,3) + t328 + (m(6) * t245 + t286 - t287) * pkin(8) + t373) * qJD(4) + t10 * qJD(5) + t243, t8 * qJD(3) + t10 * qJD(4) + (-t47 * mrSges(6,1) - t46 * mrSges(6,2) + t234) * qJD(5) + t242; t1, t280 * t30, t19 * qJD(4) + t300 + t324 + (-t229 + t232 + 0.2e1 * t270 * t347 + (-pkin(3) * t266 + t184 * t277) * m(5)) * qJD(3), t299 + t19 * qJD(3) + (t224 + t232) * qJD(4) + t324, qJD(5) * t135 + t280 * t77 - t301; qJD(4) * t5 + qJD(5) * t7 - t244, qJD(4) * t20 - t300 + t323, qJD(4) * t48 + qJD(5) * t62, ((-t329 * pkin(4) + t228 * pkin(8)) * t345 + t215) * qJD(4) + t39 * qJD(5) + t227, t39 * qJD(4) + (-t199 * t251 - t249) * qJD(5) + t351; -qJD(3) * t5 + qJD(5) * t9 - t243, -qJD(3) * t20 - t299 + t323, -qJD(5) * t38 - t227, t69 * qJD(5), (-pkin(8) * t251 - t249) * qJD(5) - t352; -qJD(3) * t7 - qJD(4) * t9 - t242, t280 * t76 + t301, qJD(4) * t38 - t351, t352, 0;];
Cq = t13;
