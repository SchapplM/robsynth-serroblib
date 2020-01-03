% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPP8_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP8_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP8_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP8_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP8_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:07:33
% EndTime: 2019-12-31 21:07:42
% DurationCPUTime: 3.91s
% Computational Cost: add. (4190->468), mult. (9212->596), div. (0->0), fcn. (7042->4), ass. (0->225)
t367 = Ifges(6,4) + Ifges(5,5);
t361 = Ifges(4,5) + Ifges(6,5);
t242 = cos(qJ(3));
t243 = cos(qJ(2));
t302 = t242 * t243;
t209 = mrSges(6,1) * t302;
t241 = sin(qJ(2));
t317 = t241 * mrSges(6,3);
t165 = t209 - t317;
t240 = sin(qJ(3));
t305 = t240 * t243;
t166 = mrSges(5,1) * t305 - t241 * mrSges(5,3);
t318 = t241 * mrSges(6,2);
t167 = -mrSges(6,1) * t305 + t318;
t210 = mrSges(5,1) * t302;
t319 = t241 * mrSges(5,2);
t168 = t210 + t319;
t239 = pkin(3) + qJ(5);
t333 = pkin(7) * t243;
t199 = t241 * pkin(2) - t333;
t303 = t242 * t199;
t335 = pkin(6) * t240;
t71 = -t303 + (-pkin(3) - t335) * t241;
t344 = -t71 / 0.2e1;
t348 = mrSges(6,3) / 0.2e1;
t350 = -m(6) / 0.2e1;
t352 = -m(5) / 0.2e1;
t286 = qJ(5) + t335;
t36 = (pkin(4) * t243 - t199) * t242 + (-pkin(3) - t286) * t241;
t160 = t240 * t199;
t334 = pkin(6) * t242;
t287 = -qJ(4) + t334;
t70 = t287 * t241 - t160;
t46 = -pkin(4) * t305 - t70;
t306 = t240 * t241;
t85 = pkin(6) * t306 + t303;
t304 = t241 * t242;
t86 = -pkin(6) * t304 + t160;
t366 = (t166 / 0.2e1 - t167 / 0.2e1) * qJ(4) + (-pkin(3) * t71 - qJ(4) * t70) * t352 + (qJ(4) * t46 - t239 * t36) * t350 + pkin(3) * t168 / 0.2e1 + t239 * t165 / 0.2e1 + t36 * t348 - t46 * mrSges(6,2) / 0.2e1 + t70 * mrSges(5,3) / 0.2e1 + mrSges(5,2) * t344 - t85 * mrSges(4,1) / 0.2e1 + t86 * mrSges(4,2) / 0.2e1;
t349 = m(6) / 0.2e1;
t365 = 0.2e1 * t349;
t364 = -t241 / 0.2e1;
t338 = t241 / 0.2e1;
t363 = -t243 / 0.2e1;
t362 = t243 / 0.2e1;
t309 = qJ(4) * t240;
t268 = -pkin(3) * t242 - t309;
t173 = -pkin(2) + t268;
t175 = t242 * mrSges(5,2) - t240 * mrSges(5,3);
t360 = m(5) * t173 + t175;
t226 = Ifges(6,6) * t240;
t359 = Ifges(6,3) * t242 + t226;
t185 = -Ifges(6,2) * t242 + t226;
t231 = Ifges(4,4) * t242;
t358 = -Ifges(4,2) * t240 + t231;
t190 = Ifges(4,1) * t240 + t231;
t357 = Ifges(5,4) * t306 + t367 * t304;
t356 = Ifges(5,1) + Ifges(6,1) + Ifges(4,3);
t322 = Ifges(5,6) * t242;
t183 = Ifges(5,3) * t240 - t322;
t321 = Ifges(6,6) * t242;
t184 = Ifges(6,2) * t240 + t321;
t295 = -Ifges(6,4) / 0.2e1 - Ifges(5,5) / 0.2e1;
t346 = -Ifges(4,6) / 0.2e1;
t284 = t346 - t295;
t354 = t284 * t241 + Ifges(4,6) * t364 + t358 * t363 + (t183 + t184) * t362 + t367 * t338;
t353 = -0.2e1 * qJ(4);
t351 = m(5) / 0.2e1;
t347 = -Ifges(5,4) / 0.2e1;
t332 = t243 * pkin(2);
t342 = pkin(4) + pkin(7);
t256 = (-t342 * t241 - pkin(1) - t332) * t240;
t297 = pkin(6) * t302;
t55 = t256 + t297;
t345 = t55 / 0.2e1;
t293 = -t241 * pkin(7) - pkin(1);
t261 = t293 - t332;
t84 = t240 * t261 + t297;
t343 = t84 / 0.2e1;
t298 = pkin(4) * t304;
t159 = t242 * t261;
t83 = pkin(6) * t305 - t159;
t54 = -t83 - t298;
t341 = m(6) * t54;
t174 = -mrSges(6,2) * t242 + mrSges(6,3) * t240;
t141 = t174 * t241;
t340 = -t141 / 0.2e1;
t197 = t342 * t240;
t339 = -t197 / 0.2e1;
t234 = t241 * pkin(6);
t331 = -mrSges(5,1) - mrSges(6,1);
t330 = -mrSges(6,2) - mrSges(5,3);
t329 = m(6) * qJD(4);
t328 = m(6) * qJD(5);
t326 = Ifges(4,4) * t240;
t323 = Ifges(5,6) * t240;
t316 = t243 * mrSges(6,2);
t315 = t243 * mrSges(5,3);
t314 = t243 * Ifges(5,4);
t313 = t243 * Ifges(4,6);
t300 = pkin(3) * t306 + t234;
t108 = -qJ(4) * t304 + t300;
t299 = pkin(3) * t305 + t243 * pkin(6);
t109 = -qJ(4) * t302 + t299;
t144 = t174 * t243;
t178 = -mrSges(5,2) * t240 - mrSges(5,3) * t242;
t146 = t178 * t241;
t147 = t178 * t243;
t179 = mrSges(4,1) * t240 + mrSges(4,2) * t242;
t148 = t243 * t179;
t161 = mrSges(4,2) * t243 - mrSges(4,3) * t306;
t162 = -t241 * mrSges(4,2) - mrSges(4,3) * t305;
t163 = -mrSges(4,1) * t243 - mrSges(4,3) * t304;
t164 = t241 * mrSges(4,1) - mrSges(4,3) * t302;
t225 = t243 * mrSges(6,3);
t296 = mrSges(6,1) * t304;
t169 = t225 + t296;
t170 = mrSges(5,1) * t306 + t315;
t171 = -mrSges(6,1) * t306 - t316;
t172 = mrSges(5,1) * t304 - mrSges(5,2) * t243;
t187 = -Ifges(5,2) * t242 + t323;
t191 = Ifges(4,1) * t242 - t326;
t294 = Ifges(6,5) / 0.2e1 + Ifges(4,5) / 0.2e1;
t285 = t347 + t294;
t255 = Ifges(5,4) * t364 + t187 * t363 + t285 * t241 + (t191 + t359) * t362 + t361 * t338;
t103 = -Ifges(5,5) * t243 + t241 * t183;
t211 = Ifges(6,6) * t304;
t105 = -Ifges(6,4) * t243 + Ifges(6,2) * t306 + t211;
t96 = t241 * t358 - t313;
t280 = t103 / 0.2e1 + t105 / 0.2e1 - t96 / 0.2e1;
t101 = -t243 * Ifges(6,5) + t241 * t359;
t212 = Ifges(5,6) * t306;
t107 = -Ifges(5,2) * t304 + t212 - t314;
t98 = -t243 * Ifges(4,5) + t241 * t191;
t282 = t101 / 0.2e1 - t107 / 0.2e1 + t98 / 0.2e1;
t236 = t243 * pkin(3);
t260 = t159 - t236 - t298;
t34 = t286 * t243 - t260;
t44 = t287 * t243 + t256;
t308 = qJ(4) * t242;
t265 = qJ(5) * t240 - t308;
t65 = t265 * t241 + t300;
t66 = t265 * t243 + t299;
t278 = -pkin(2) * t240 + t334;
t283 = t240 * t293;
t68 = -t283 + (qJ(4) - t278) * t243;
t69 = t236 + t83;
t3 = t66 * t141 + t65 * t144 + t109 * t146 + t108 * t147 + t86 * t161 + t84 * t162 + t85 * t163 - t83 * t164 + t34 * t165 + t68 * t166 + t44 * t167 + t69 * t168 + t36 * t169 + t70 * t170 + t46 * t171 + t71 * t172 + m(4) * (-t83 * t85 + t84 * t86) + m(6) * (t34 * t36 + t44 * t46 + t65 * t66) + m(5) * (t108 * t109 + t68 * t70 + t69 * t71) + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t243 + (-t285 * t243 + t282) * t242 + (-t284 * t243 + t280) * t240) * t243 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t241 + pkin(6) * t148 + t255 * t242 + t354 * t240 + (Ifges(3,1) - Ifges(3,2) + (m(4) * pkin(6) + t179) * pkin(6) - t356) * t243) * t241;
t312 = t3 * qJD(1);
t140 = pkin(3) * t304 + qJ(4) * t306;
t142 = t175 * t241;
t273 = t242 * mrSges(4,1) - t240 * mrSges(4,2);
t143 = t273 * t241;
t272 = mrSges(6,2) * t240 + mrSges(6,3) * t242;
t145 = t272 * t241;
t149 = -Ifges(6,3) * t306 + t211;
t150 = Ifges(5,3) * t304 + t212;
t151 = t241 * t185;
t270 = Ifges(5,2) * t240 + t322;
t152 = t270 * t241;
t188 = Ifges(4,2) * t242 + t326;
t153 = t241 * t188;
t154 = t241 * t190;
t82 = qJ(5) * t304 + t140;
t4 = t82 * t141 - t108 * t142 + t65 * t145 + t140 * t146 + t55 * t169 + t54 * t171 + (-t163 + t172) * t84 + (-t161 + t170) * t83 + m(6) * (t34 * t55 + t44 * t54 + t65 * t82) + m(5) * (t108 * t140 + t68 * t83 + t69 * t84) + (pkin(6) * t143 + (t68 * mrSges(5,1) + t149 / 0.2e1 - t152 / 0.2e1 - t154 / 0.2e1 - t84 * mrSges(4,3) - t44 * mrSges(6,1) + t313 / 0.2e1 + t280) * t242 + (-t69 * mrSges(5,1) - t83 * mrSges(4,3) - t34 * mrSges(6,1) + t150 / 0.2e1 - t151 / 0.2e1 + t153 / 0.2e1 + t294 * t243 - t282) * t240) * t241 + t357 * t363;
t311 = t4 * qJD(1);
t9 = (m(5) * t108 + m(6) * t65 + t141 + t146) * t304 + (-m(5) * t68 + m(6) * t44 - t170 + t171) * t243;
t310 = t9 * qJD(1);
t15 = m(6) * (t243 * t34 + t65 * t306) + t243 * t169 + t141 * t306;
t307 = t15 * qJD(1);
t291 = t272 * t338;
t139 = -t239 * t242 - pkin(2) - t309;
t288 = m(6) * t139 - t272;
t279 = -pkin(3) * mrSges(5,1) - t239 * mrSges(6,1);
t277 = (t272 / 0.2e1 - t175 / 0.2e1) * t242;
t276 = t331 * qJ(4) - Ifges(4,6);
t181 = Ifges(6,3) * t240 - t321;
t275 = t181 / 0.2e1 + t270 / 0.2e1 + t190 / 0.2e1;
t182 = -Ifges(5,3) * t242 - t323;
t274 = t182 / 0.2e1 + t185 / 0.2e1 - t188 / 0.2e1;
t233 = t240 * pkin(3);
t158 = t233 + t265;
t176 = t233 - t308;
t198 = t342 * t242;
t227 = Ifges(6,5) * t242;
t228 = Ifges(5,5) * t240;
t229 = Ifges(4,5) * t242;
t230 = Ifges(6,4) * t240;
t244 = (-t228 / 0.4e1 - t229 / 0.4e1 - t230 / 0.4e1 - t227 / 0.4e1) * t243 + (t108 * t176 + t140 * t173) * t351 + (t139 * t82 + t158 * t65 + (t34 + t54) * t198 + (-t44 + t55) * t197) * t349 - pkin(2) * t143 / 0.2e1 + t108 * t178 / 0.2e1 + t139 * t145 / 0.2e1 + t140 * t175 / 0.2e1 + t158 * t141 / 0.2e1 - t173 * t142 / 0.2e1 + t176 * t146 / 0.2e1 + t171 * t339 + t198 * t169 / 0.2e1 + t65 * t174 / 0.2e1 - t82 * t272 / 0.2e1;
t248 = (t345 - t44 / 0.2e1) * mrSges(6,1) + (t343 + t68 / 0.2e1) * mrSges(5,1) + (t170 / 0.2e1 - t161 / 0.2e1 + (t68 + t84) * t351) * pkin(7) + t103 / 0.4e1 + t105 / 0.4e1 + t149 / 0.4e1 - t152 / 0.4e1 - t154 / 0.4e1 - t96 / 0.4e1;
t249 = (t54 / 0.2e1 + t34 / 0.2e1) * mrSges(6,1) + (-t83 / 0.2e1 + t69 / 0.2e1) * mrSges(5,1) + (t172 / 0.2e1 - t163 / 0.2e1 + (t69 - t83) * t351) * pkin(7) + t101 / 0.4e1 - t107 / 0.4e1 - t150 / 0.4e1 + t151 / 0.4e1 - t153 / 0.4e1 + t98 / 0.4e1;
t251 = pkin(6) * t179 / 0.2e1 + (mrSges(5,1) + mrSges(4,3)) * pkin(7) * (-t242 ^ 2 / 0.2e1 - t240 ^ 2 / 0.2e1);
t253 = -t198 * mrSges(6,1) / 0.2e1 + t191 / 0.4e1 - t188 / 0.4e1 - t187 / 0.4e1 + t185 / 0.4e1 + t182 / 0.4e1 + t359 / 0.4e1;
t254 = mrSges(6,1) * t339 - t190 / 0.4e1 - t358 / 0.4e1 - t270 / 0.4e1 + t184 / 0.4e1 + t183 / 0.4e1 - t181 / 0.4e1;
t2 = ((0.3e1 / 0.4e1 * Ifges(4,6) + t295) * t243 + t254 * t241 + t248) * t240 + ((0.3e1 / 0.4e1 * Ifges(5,4) - t294) * t243 + t253 * t241 + t249) * t242 + (-Ifges(5,1) / 0.2e1 - Ifges(6,1) / 0.2e1 - Ifges(4,3) / 0.2e1 + t251) * t241 + t244 + t366;
t5 = -pkin(2) * t179 - t158 * t272 + t176 * t175 + (m(5) * t176 + t178) * t173 + (m(6) * t158 + t174) * t139 + (t358 / 0.2e1 - t184 / 0.2e1 - t183 / 0.2e1 + t275) * t242 + (t359 / 0.2e1 - t187 / 0.2e1 + t191 / 0.2e1 + t274) * t240;
t267 = -t2 * qJD(1) - t5 * qJD(2);
t22 = (t288 + t360) * t240;
t247 = (t340 - t146 / 0.2e1) * t240 + (-t240 * t108 + (-t173 * t241 - t333) * t242) * t351 + (-t139 * t304 - t198 * t243 - t240 * t65) * t349;
t259 = m(5) * t344 + t36 * t350;
t6 = -t210 - t209 + (t348 - mrSges(5,2) / 0.2e1 + t277) * t241 + t247 + t259;
t266 = -qJD(1) * t6 + qJD(2) * t22;
t252 = (t139 * t306 + t197 * t243 - t242 * t65) * t349 + t242 * t340;
t257 = t46 * t349 + t318 / 0.2e1;
t11 = (t243 * mrSges(6,1) - t291) * t240 + t252 - t257;
t41 = t288 * t242;
t264 = -qJD(1) * t11 + qJD(2) * t41;
t246 = (t283 + (t353 + t278) * t243) * t352 + ((t353 + t334) * t243 + t256) * t350;
t258 = m(5) * t343 + m(6) * t345;
t14 = t246 + t258 + t315 + t316;
t201 = (m(5) + m(6)) * qJ(4) - t330;
t263 = -qJD(1) * t14 + qJD(3) * t201;
t250 = -t225 + ((-t239 - t286) * t243 + t260) * t349;
t17 = -t341 / 0.2e1 + t250;
t200 = m(6) * t239 + mrSges(6,3);
t262 = qJD(1) * t17 + qJD(3) * t200;
t192 = (qJD(1) * t243 - qJD(3)) * m(6);
t110 = -m(6) * t197 - t240 * mrSges(6,1);
t53 = m(6) * t198 + (m(5) * pkin(7) - t331) * t242;
t16 = -t296 + t341 / 0.2e1 + t250;
t13 = -t240 * t291 + t252 + t257;
t12 = t330 * t243 + t331 * t306 - t246 + t258;
t7 = -t317 / 0.2e1 + t319 / 0.2e1 + t241 * t277 + t247 - t259;
t1 = (t254 * t240 + t253 * t242 + t251) * t241 + (t314 / 0.4e1 + t249) * t242 + (t313 / 0.4e1 + t248) * t240 + t244 + t356 * t338 + (t347 + t361 / 0.2e1) * t302 - t366 + (t346 + t367 / 0.2e1) * t305;
t8 = [qJD(2) * t3 + qJD(3) * t4 - qJD(4) * t9 + qJD(5) * t15, t1 * qJD(3) + t7 * qJD(4) + t13 * qJD(5) + t312 + (mrSges(3,2) * t234 - Ifges(3,6) * t241 - pkin(2) * t148 + t139 * t144 + t173 * t147 + t197 * t165 + t198 * t167 - t66 * t272 + (t139 * t66 + t197 * t36 + t198 * t46) * t365 + (-t70 * mrSges(5,1) + t46 * mrSges(6,1) + t86 * mrSges(4,3) - t354) * t242 + (t71 * mrSges(5,1) + t36 * mrSges(6,1) - t85 * mrSges(4,3) + t255) * t240 + ((t162 - t166) * t242 + (-t164 + t168) * t240 + m(5) * (t240 * t71 - t242 * t70) + m(4) * (-t240 * t85 + t242 * t86)) * pkin(7) + (Ifges(3,5) + t275 * t242 + t274 * t240 + (-m(4) * pkin(2) - mrSges(3,1) - t273) * pkin(6)) * t243 + t360 * t109) * qJD(2), t1 * qJD(2) + t12 * qJD(4) + t16 * qJD(5) + t311 + (-t84 * mrSges(4,1) + t83 * mrSges(4,2) + t84 * mrSges(5,2) + t54 * mrSges(6,2) - t83 * mrSges(5,3) - t55 * mrSges(6,3) + (qJ(4) * t54 - t239 * t55) * t365 + 0.2e1 * (-pkin(3) * t84 - qJ(4) * t83) * t351 + (t276 * t242 + (-t279 - t361) * t240) * t241 + t357) * qJD(3), qJD(2) * t7 + qJD(3) * t12 - t310, qJD(2) * t13 + qJD(3) * t16 + t307; qJD(3) * t2 + qJD(4) * t6 + qJD(5) * t11 - t312, qJD(3) * t5 - qJD(4) * t22 - qJD(5) * t41, t53 * qJD(4) + t110 * qJD(5) - t267 + (m(6) * (-qJ(4) * t197 - t198 * t239) + t230 + t227 + t228 + t229 - t197 * mrSges(6,2) - t198 * mrSges(6,3) + (-Ifges(5,4) + t279) * t242 + t276 * t240 + (m(5) * t268 + t175 - t273) * pkin(7)) * qJD(3), qJD(3) * t53 - t266, qJD(3) * t110 - t264; -qJD(2) * t2 - qJD(4) * t14 + qJD(5) * t17 - t311, t267, qJD(4) * t201 + qJD(5) * t200, t263, t262; -t6 * qJD(2) + t14 * qJD(3) + t243 * t328 + t310, t266, -t263 - t328, 0, t192; -t11 * qJD(2) - t17 * qJD(3) - t243 * t329 - t307, t264, -t262 + t329, -t192, 0;];
Cq = t8;
