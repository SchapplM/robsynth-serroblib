% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRR12_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR12_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR12_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR12_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR12_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:12:32
% EndTime: 2019-12-31 19:12:38
% DurationCPUTime: 2.86s
% Computational Cost: add. (7710->315), mult. (14296->419), div. (0->0), fcn. (14313->6), ass. (0->207)
t206 = sin(qJ(4));
t207 = sin(qJ(3));
t335 = cos(qJ(4));
t336 = cos(qJ(3));
t183 = t206 * t336 + t335 * t207;
t333 = m(6) * t183;
t376 = -t333 / 0.2e1;
t182 = t206 * t207 - t335 * t336;
t209 = -pkin(1) - pkin(6);
t185 = (-pkin(7) + t209) * t207;
t271 = t336 * t209;
t234 = -t336 * pkin(7) + t271;
t221 = t335 * t185 + t206 * t234;
t316 = t221 * mrSges(5,3);
t375 = t182 * t316;
t208 = cos(qJ(5));
t200 = Ifges(6,5) * t208;
t205 = sin(qJ(5));
t321 = Ifges(6,6) * t205;
t256 = t321 - t200;
t374 = t183 * t256;
t203 = t205 ^ 2;
t204 = t208 ^ 2;
t358 = t203 + t204;
t264 = t358 * t182;
t323 = Ifges(6,4) * t205;
t191 = t208 * Ifges(6,1) - t323;
t100 = -Ifges(6,5) * t182 - t183 * t191;
t322 = Ifges(6,2) * t208;
t188 = t322 + t323;
t201 = Ifges(6,4) * t208;
t325 = Ifges(6,1) * t205;
t190 = t201 + t325;
t255 = Ifges(6,5) * t205 + Ifges(6,6) * t208;
t238 = t182 * t255;
t290 = t183 * t208;
t268 = -t290 / 0.2e1;
t291 = t183 * t205;
t269 = t291 / 0.2e1;
t337 = t208 / 0.2e1;
t338 = t205 / 0.2e1;
t257 = t205 * Ifges(6,2) - t201;
t98 = -Ifges(6,6) * t182 + t183 * t257;
t235 = t100 * t338 + t188 * t269 + t190 * t268 + t98 * t337 + Ifges(5,6) * t182 - Ifges(5,5) * t183 - t238 / 0.2e1;
t307 = t208 * mrSges(6,1);
t313 = t205 * mrSges(6,2);
t186 = -t307 + t313;
t359 = t221 * t186;
t362 = t221 * mrSges(5,1);
t142 = t185 * t206 - t335 * t234;
t370 = t142 * mrSges(5,2);
t373 = t235 + t359 - t362 + t370;
t372 = -t359 / 0.2e1 + t362 / 0.2e1 - t370 / 0.2e1;
t194 = t207 * pkin(3) + qJ(2);
t371 = m(5) * t194;
t369 = t142 * t205;
t368 = t142 * t206;
t367 = t142 * t208;
t298 = t142 * t221;
t339 = -t205 / 0.2e1;
t366 = t188 * t339 + t190 * t337;
t363 = -Ifges(5,1) + Ifges(6,3);
t276 = t335 * pkin(3);
t197 = -t276 - pkin(4);
t360 = t197 * t221;
t306 = t208 * mrSges(6,2);
t314 = t205 * mrSges(6,1);
t187 = t306 + t314;
t242 = t306 / 0.2e1 + t314 / 0.2e1;
t79 = (-t187 / 0.2e1 + t242) * t182;
t283 = t79 * qJD(2);
t261 = -mrSges(6,2) * t335 / 0.2e1;
t347 = -mrSges(6,1) / 0.2e1;
t39 = (-t197 / 0.2e1 + pkin(4) / 0.2e1) * t187 + (pkin(3) * t261 - t190 / 0.2e1 + t257 / 0.2e1) * t208 + (t276 * t347 - t191 / 0.2e1 + t188 / 0.2e1) * t205;
t262 = t191 * t338 - t257 * t337 + t366;
t330 = pkin(4) * t187;
t67 = t262 - t330;
t218 = (t188 / 0.4e1 - t191 / 0.4e1 + t322 / 0.4e1) * t208 + (-t257 / 0.4e1 + t190 / 0.4e1 + t201 / 0.2e1 + t325 / 0.4e1) * t205;
t259 = mrSges(6,3) * (t204 / 0.2e1 + t203 / 0.2e1);
t213 = pkin(8) * t259 + t218;
t130 = t186 * t182;
t101 = t183 * Ifges(6,5) - t182 * t191;
t286 = t208 * t101;
t99 = t183 * Ifges(6,6) + t182 * t257;
t309 = t205 * t99;
t340 = t187 / 0.2e1;
t223 = t142 * t340 - t309 / 0.4e1 + t286 / 0.4e1;
t292 = t182 * t208;
t134 = mrSges(6,1) * t183 + mrSges(6,3) * t292;
t285 = t208 * t134;
t294 = t182 * t205;
t275 = mrSges(6,3) * t294;
t133 = -mrSges(6,2) * t183 + t275;
t288 = t205 * t133;
t240 = -t288 / 0.2e1 - t285 / 0.2e1;
t217 = t240 * pkin(8) - pkin(4) * t130 / 0.2e1 + t223;
t273 = t200 / 0.2e1;
t224 = (-0.3e1 / 0.4e1 * t321 + t200 / 0.4e1 + t273) * t183;
t346 = mrSges(6,2) / 0.2e1;
t141 = -pkin(4) * t182 + t183 * pkin(8);
t68 = t141 * t208 + t369;
t69 = t141 * t205 - t367;
t247 = t69 * t346 + t68 * t347;
t345 = Ifges(6,3) / 0.2e1;
t8 = t224 + (t345 + t213) * t182 + t217 + t247;
t357 = -t8 * qJD(1) + t39 * qJD(3) - t67 * qJD(4) + t283;
t332 = pkin(3) * t206;
t196 = pkin(8) + t332;
t212 = t196 * t259 + t218;
t216 = t240 * t196 + t197 * t130 / 0.2e1 + t223;
t277 = t336 * pkin(3);
t132 = t277 + t141;
t57 = t132 * t208 + t369;
t58 = t132 * t205 - t367;
t248 = t58 * t346 + t57 * t347;
t5 = t224 + (t345 + t212) * t182 + t216 + t248;
t289 = t197 * t187;
t54 = t262 + t289;
t356 = t5 * qJD(1) + t54 * qJD(3) - t283;
t293 = t182 * t206;
t355 = t335 * t183 + t293;
t159 = t183 * t186;
t354 = -mrSges(6,3) * t264 + t159;
t237 = t183 * t187;
t239 = t182 * t187;
t245 = mrSges(6,2) * t182 + mrSges(6,3) * t291;
t246 = -mrSges(6,1) * t182 + mrSges(6,3) * t290;
t331 = pkin(4) * t183;
t131 = pkin(8) * t182 + t194 + t331;
t55 = t131 * t208 - t205 * t221;
t56 = t131 * t205 + t208 * t221;
t353 = t221 * t239 + t142 * t237 - t375 - t56 * t245 - t55 * t246 - t194 * (-mrSges(5,1) * t182 - mrSges(5,2) * t183);
t320 = Ifges(6,3) * t182;
t352 = Ifges(6,5) * t268 + Ifges(6,6) * t269 - t374 / 0.4e1 - t320 / 0.2e1;
t267 = t133 * t337;
t305 = t208 * t56;
t350 = -m(6) / 0.2e1;
t351 = ((t205 * t55 + t221 - t305) * t350 + t267 + t134 * t339 + t237 / 0.2e1) * t182;
t349 = m(5) * pkin(3);
t348 = m(6) * pkin(3);
t342 = t182 / 0.2e1;
t341 = -t183 / 0.2e1;
t334 = m(6) * (-pkin(8) * t264 - t331);
t304 = t208 * t58;
t311 = t205 * t57;
t253 = t304 - t311;
t11 = (t142 + t253) * t376 + t351;
t303 = t208 * t69;
t310 = t205 * t68;
t252 = t303 - t310;
t13 = (t142 + t252) * t376 + t351;
t329 = -t11 * qJD(3) - t13 * qJD(4);
t32 = (0.1e1 - t358) * t182 * t333;
t284 = t32 * qJD(2);
t80 = (t242 + t340) * t182;
t328 = t80 * qJD(5) + t284;
t327 = -t79 * qJD(5) - t284;
t324 = Ifges(5,4) * t182;
t308 = t207 * mrSges(4,1);
t7 = t56 * t134 - t142 * t130 + (-t208 * t99 / 0.2e1 + t101 * t339 - mrSges(6,3) * t305 + t255 * t341 + t366 * t182) * t182 + (-t133 + t275) * t55;
t302 = t7 * qJD(1);
t266 = -t183 * mrSges(5,1) + t182 * mrSges(5,2);
t222 = -t336 * mrSges(4,2) + t266 - t308;
t22 = t288 + t285 + mrSges(3,3) + (m(4) + m(3)) * qJ(2) + m(6) * (t205 * t56 + t208 * t55) + t371 - t222;
t301 = qJD(1) * t22;
t300 = t11 * qJD(1);
t299 = t13 * qJD(1);
t219 = (t182 * t259 + t240) * t183 + t130 * t342;
t243 = t313 / 0.2e1 - t307 / 0.2e1;
t17 = t219 + t243;
t295 = t17 * qJD(1);
t287 = t208 * t100;
t280 = qJD(3) + qJD(4);
t279 = mrSges(6,3) * t311;
t278 = mrSges(6,3) * t304;
t274 = t334 / 0.2e1;
t272 = t205 * t335;
t260 = -t272 / 0.2e1;
t175 = Ifges(5,4) * t183;
t1 = t183 * (-t320 + t374) / 0.2e1 + m(6) * (t55 * t57 + t56 * t58 + t298) + (-qJ(2) * mrSges(4,2) + Ifges(4,4) * t207) * t207 - t375 + t101 * t268 + t98 * t294 / 0.2e1 + t99 * t269 + t58 * t133 + t57 * t134 + (-Ifges(5,2) * t183 - t324) * t342 + (-0.2e1 * t175 + (-Ifges(5,1) + Ifges(5,2)) * t182) * t341 - (t182 * t256 + t363 * t183 + t287 + t324) * t182 / 0.2e1 + ((-Ifges(4,1) + Ifges(4,2)) * t207 + qJ(2) * mrSges(4,1) - Ifges(4,4) * t336) * t336 - t353 + (t371 - t266) * t277;
t251 = t1 * qJD(1) - t11 * qJD(2);
t244 = t273 - t321 / 0.2e1;
t3 = -t69 * t133 - t68 * t134 - m(6) * (t55 * t68 + t56 * t69 + t298) + (-t309 / 0.2e1 + t286 / 0.2e1 - t175 + t244 * t183) * t183 + (t98 * t339 + t287 / 0.2e1 + t316 + (Ifges(5,4) - t244) * t182 + (Ifges(5,2) + t363) * t183) * t182 + t353;
t250 = -t3 * qJD(1) - t13 * qJD(2);
t236 = t266 + t354;
t233 = t205 * t246;
t232 = t208 * t245;
t231 = t358 * t335;
t230 = -t237 / 0.2e1;
t227 = t197 * t183 - t196 * t264;
t214 = m(6) * ((t183 * t231 + t293) * pkin(3) + t227);
t23 = t274 - t214 / 0.2e1;
t225 = t196 * t232;
t226 = t196 * t233;
t210 = m(6) * (t360 + t252 * t196 + (-t55 * t272 + t335 * t305 + t368) * pkin(3)) / 0.2e1 + t197 * t230 - t226 / 0.2e1 + t225 / 0.2e1 + pkin(3) * t134 * t260 - t239 * t332 / 0.2e1 + t267 * t276 + (-t310 / 0.2e1 + t303 / 0.2e1) * mrSges(6,3) - t372;
t211 = t279 / 0.2e1 - t278 / 0.2e1 + (-t221 * t350 + t230) * pkin(4) + (t253 * t350 + t233 / 0.2e1 - t232 / 0.2e1) * pkin(8) + t372;
t4 = t210 + t211;
t215 = (-mrSges(5,1) + t186) * t332 + (t358 * mrSges(6,3) - mrSges(5,2)) * t276;
t44 = (t196 * t231 + t197 * t206) * t348 + t215;
t229 = t4 * qJD(1) - t23 * qJD(2) + t44 * qJD(3);
t40 = t289 / 0.2e1 - t330 / 0.2e1 + (mrSges(6,1) * t260 + t208 * t261) * pkin(3) + t262;
t20 = t214 / 0.2e1 + t274 + t236;
t16 = t219 - t243;
t9 = t182 * t213 + t217 - t247 + t352;
t6 = t182 * t212 + t216 - t248 + t352;
t2 = -t211 + t210 + t235;
t10 = [qJD(2) * t22 + qJD(3) * t1 - qJD(4) * t3 - qJD(5) * t7, qJD(5) * t16 + t301 + t329, (-t226 + t225 - Ifges(4,6) * t336 + t278 - mrSges(4,2) * t271 - t209 * t308 - t279 + (-t221 * t335 - t368) * t349 - t197 * t237 - Ifges(4,5) * t207 + m(6) * (t196 * t253 + t360) + t355 * mrSges(5,3) * pkin(3) + t373) * qJD(3) + t2 * qJD(4) + t6 * qJD(5) + t251, t2 * qJD(3) + t9 * qJD(5) + t250 + (t252 * mrSges(6,3) + (-m(6) * t221 + t237) * pkin(4) + (m(6) * t252 + mrSges(6,1) * t294 + mrSges(6,2) * t292) * pkin(8) + t373) * qJD(4), -t302 + t16 * qJD(2) + t6 * qJD(3) + t9 * qJD(4) + (-mrSges(6,1) * t56 - mrSges(6,2) * t55 + t238) * qJD(5); qJD(5) * t17 - t301 + t329, t280 * t32, -t300 + (m(6) * t227 - t355 * t349 + t222 + t354) * qJD(3) + t20 * qJD(4) + t328, -t299 + t20 * qJD(3) + (t236 + t334) * qJD(4) + t328, qJD(5) * t159 + t280 * t80 + t295; qJD(4) * t4 + qJD(5) * t5 - t251, -qJD(4) * t23 + t300 + t327, qJD(4) * t44 + qJD(5) * t54, ((-pkin(4) * t206 + pkin(8) * t231) * t348 + t215) * qJD(4) + t40 * qJD(5) + t229, t40 * qJD(4) + (t186 * t196 - t256) * qJD(5) + t356; -qJD(3) * t4 + qJD(5) * t8 - t250, qJD(3) * t23 + t299 + t327, -qJD(5) * t39 - t229, t67 * qJD(5), (pkin(8) * t186 - t256) * qJD(5) - t357; -qJD(2) * t17 - qJD(3) * t5 - qJD(4) * t8 + t302, t280 * t79 - t295, qJD(4) * t39 - t356, t357, 0;];
Cq = t10;
