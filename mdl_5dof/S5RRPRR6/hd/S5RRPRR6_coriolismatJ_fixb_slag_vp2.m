% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% m [6x1]
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
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR6_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR6_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:16:53
% EndTime: 2022-01-20 11:17:01
% DurationCPUTime: 3.58s
% Computational Cost: add. (10609->343), mult. (21889->452), div. (0->0), fcn. (20619->8), ass. (0->210)
t377 = -mrSges(5,1) / 0.2e1;
t225 = sin(qJ(5));
t226 = sin(qJ(4));
t228 = cos(qJ(5));
t229 = cos(qJ(4));
t203 = -t225 * t229 - t228 * t226;
t223 = sin(pkin(9));
t179 = t203 * t223;
t366 = t179 / 0.2e1;
t376 = 0.2e1 * t366;
t227 = sin(qJ(2));
t224 = cos(pkin(9));
t230 = cos(qJ(2));
t289 = t224 * t230;
t186 = (-t226 * t289 + t227 * t229) * pkin(1);
t187 = (t226 * t227 + t229 * t289) * pkin(1);
t109 = t228 * t186 - t225 * t187;
t110 = t225 * t186 + t228 * t187;
t346 = -mrSges(5,2) / 0.2e1;
t347 = m(6) / 0.2e1;
t373 = -mrSges(6,2) / 0.2e1;
t374 = mrSges(6,1) / 0.2e1;
t359 = t109 * t374 + t110 * t373;
t375 = t186 * t377 - t187 * t346 - (t109 * t228 + t110 * t225) * pkin(4) * t347 - t359;
t202 = -t225 * t226 + t228 * t229;
t178 = t202 * t223;
t322 = t178 * mrSges(6,3);
t150 = -t224 * mrSges(6,1) - t322;
t372 = t150 + t322;
t259 = t203 * mrSges(6,1) - t202 * mrSges(6,2);
t371 = qJD(5) * t259;
t295 = t223 * t226;
t272 = mrSges(5,3) * t295;
t199 = t224 * mrSges(5,2) - t272;
t280 = t229 * t199;
t294 = t223 * t229;
t270 = mrSges(5,3) * t294;
t200 = -t224 * mrSges(5,1) - t270;
t284 = t226 * t200;
t357 = -t284 / 0.2e1 + t280 / 0.2e1;
t205 = -t224 * pkin(3) - t223 * pkin(7) - pkin(2);
t196 = t229 * t205;
t274 = pkin(8) * t294;
t134 = -t274 + t196 + (-qJ(3) * t226 - pkin(4)) * t224;
t290 = t224 * t229;
t160 = qJ(3) * t290 + t226 * t205;
t212 = pkin(8) * t295;
t143 = t160 - t212;
t287 = t225 * t143;
t70 = t228 * t134 - t287;
t282 = t228 * t143;
t71 = t225 * t134 + t282;
t291 = t224 * t226;
t159 = -qJ(3) * t291 + t196;
t142 = t159 - t274;
t80 = -t225 * t142 - t282;
t81 = t228 * t142 - t287;
t370 = -m(6) * ((t70 - t81) * t203 + (t71 + t80) * t202) / 0.2e1 - t357;
t112 = -t179 * mrSges(6,1) + t178 * mrSges(6,2);
t275 = pkin(4) * t294;
t368 = t112 * t275 - t224 * (-Ifges(5,5) * t226 - Ifges(5,6) * t229) * t223 / 0.2e1 - (-Ifges(5,5) * t224 + (-0.2e1 * Ifges(5,4) * t226 + (Ifges(5,1) - Ifges(5,2)) * t229) * t223) * t295 / 0.2e1;
t321 = t179 * mrSges(6,3);
t149 = t224 * mrSges(6,2) + t321;
t367 = t149 / 0.2e1;
t113 = Ifges(6,5) * t179 - Ifges(6,6) * t178;
t255 = Ifges(6,4) * t179 * t376 + (-Ifges(6,5) * t366 - t113 / 0.2e1) * t224 + (Ifges(6,6) * t224 / 0.2e1 - Ifges(6,4) * t178 + (Ifges(6,1) - Ifges(6,2)) * t376) * t178;
t60 = t71 * t322;
t111 = t178 * mrSges(6,1) + t179 * mrSges(6,2);
t213 = pkin(4) * t295;
t297 = t223 * qJ(3);
t201 = t213 + t297;
t91 = t201 * t111;
t364 = t255 - t60 + t91;
t265 = -t321 / 0.2e1;
t352 = t202 * t367 + t372 * t203 / 0.2e1;
t257 = t202 * t265 + t352;
t180 = t203 * t224;
t181 = t202 * t224;
t279 = t180 * t374 + t181 * t373;
t355 = t257 + t279;
t363 = t355 * qJD(3);
t361 = qJD(5) * t355 + m(6) * (t202 * t180 - t203 * t181) * qJD(3);
t221 = t223 ^ 2;
t222 = t224 ^ 2;
t317 = t229 * mrSges(5,2);
t318 = t226 * mrSges(5,1);
t254 = t317 + t318;
t360 = (t254 * t223 + t112) * t223 + (t221 + t222) * mrSges(4,3);
t338 = t230 * pkin(1);
t194 = t205 - t338;
t182 = t229 * t194;
t339 = t227 * pkin(1);
t215 = qJ(3) + t339;
t135 = -t215 * t291 + t182;
t136 = t226 * t194 + t215 * t290;
t207 = t221 * t215;
t218 = t221 * qJ(3);
t278 = t207 + t218;
t356 = ((t136 + t160) * t229 + (-t135 - t159) * t226) * t224 + t278;
t343 = m(6) * (t180 * t228 + t181 * t225) * pkin(4);
t353 = t291 * t377 + t290 * t346 + t343 / 0.2e1;
t327 = Ifges(5,4) * t229;
t174 = -Ifges(5,6) * t224 + (-Ifges(5,2) * t226 + t327) * t223;
t193 = (-Ifges(5,1) * t226 - t327) * t223;
t260 = t294 / 0.2e1;
t261 = -t294 / 0.2e1;
t351 = t174 * t261 + t193 * t260 + t368;
t350 = m(4) / 0.2e1;
t348 = m(5) / 0.2e1;
t108 = -t274 + t182 + (-t215 * t226 - pkin(4)) * t224;
t123 = t136 - t212;
t288 = t225 * t123;
t61 = t228 * t108 - t288;
t345 = t61 / 0.2e1;
t344 = -t71 / 0.2e1;
t341 = pkin(4) * t225;
t340 = pkin(4) * t228;
t337 = t61 * mrSges(6,2);
t283 = t228 * t123;
t62 = t225 * t108 + t283;
t336 = t62 * mrSges(6,1);
t122 = t135 - t274;
t65 = -t225 * t122 - t283;
t335 = t65 * mrSges(6,1);
t66 = t228 * t122 - t288;
t334 = t66 * mrSges(6,2);
t333 = t70 * mrSges(6,2);
t332 = t71 * mrSges(6,1);
t331 = t80 * mrSges(6,1);
t330 = t81 * mrSges(6,2);
t325 = pkin(4) * qJD(4);
t296 = t223 * t215;
t183 = t213 + t296;
t153 = t183 * t275;
t189 = (mrSges(5,1) * t229 - mrSges(5,2) * t226) * t223;
t301 = t183 * t111;
t45 = t61 * t321;
t246 = t255 - t45 + t301;
t307 = t135 * t226;
t248 = -t136 * t229 + t307;
t306 = t136 * t200;
t308 = t135 * t199;
t313 = t66 * t149;
t314 = t65 * t150;
t315 = t62 * t178;
t5 = t246 + m(6) * (t61 * t65 + t62 * t66 + t153) + (t248 * mrSges(5,3) + t215 * t189) * t223 - mrSges(6,3) * t315 + t313 + t314 + t308 - t306 + t351;
t316 = t5 * qJD(1);
t7 = t61 * t149 - t372 * t62 + t246;
t312 = t7 * qJD(1);
t311 = t70 * t179;
t241 = t181 * t149 + t180 * t150 + (t280 - t284) * t224 + t360;
t23 = m(6) * (t61 * t180 + t62 * t181 + t183 * t223) + m(5) * (-t248 * t224 + t207) + m(4) * (t222 * t215 + t207) + t241;
t310 = qJD(1) * t23;
t277 = t221 * t338;
t197 = t215 * t277;
t236 = ((-t224 * mrSges(4,1) + t223 * mrSges(4,2) - mrSges(3,1)) * t227 + (-mrSges(3,2) + t360) * t230) * pkin(1);
t242 = t109 * t150 + t110 * t149 + t186 * t200 + t187 * t199;
t276 = t223 * t338;
t298 = t222 * t230;
t13 = t236 + m(6) * (t61 * t109 + t62 * t110 + t183 * t276) + m(5) * (t135 * t186 + t136 * t187 + t197) + m(4) * (t215 * pkin(1) * t298 + t197 + (-pkin(2) - t338) * t339) + t242;
t309 = t13 * qJD(1);
t286 = t225 * t150;
t271 = t228 * t321;
t264 = t321 / 0.2e1;
t258 = t201 * t275;
t154 = t322 * t341;
t256 = -t154 / 0.2e1 - pkin(4) * t286 / 0.2e1 + (t367 + t265) * t340;
t126 = t159 * t199;
t127 = t160 * t200;
t141 = t159 * t272;
t170 = t189 * t297;
t234 = (-t315 / 0.2e1 - t311 / 0.2e1) * mrSges(6,3) - t45 / 0.2e1 - t60 / 0.2e1 + t91 / 0.2e1 + t301 / 0.2e1 + t255;
t54 = t80 * t150;
t55 = t81 * t149;
t231 = t189 * t296 / 0.2e1 + t170 / 0.2e1 + (t80 * t61 + t81 * t62 + t70 * t65 + t71 * t66 + t153 + t258) * t347 + t351 + t234 + t126 / 0.2e1 - t127 / 0.2e1 + t313 / 0.2e1 + t314 / 0.2e1 + t55 / 0.2e1 + t54 / 0.2e1 - t306 / 0.2e1 + t308 / 0.2e1 + t141 / 0.2e1 + (t307 / 0.2e1 + (-t160 / 0.2e1 - t136 / 0.2e1) * t229) * t223 * mrSges(5,3);
t2 = t231 + t375;
t6 = t174 * t260 + t193 * t261 + t160 * t270 - t170 + mrSges(6,3) * t311 - t126 + t127 - t55 - t54 - t141 - m(6) * (t70 * t80 + t71 * t81 + t258) - t364 - t368;
t253 = t2 * qJD(1) - t6 * qJD(2);
t233 = (t345 + t70 / 0.2e1) * t149 + (-t62 / 0.2e1 + t344) * t150 + t234;
t4 = t233 - t359;
t8 = -t71 * t150 + (t149 - t321) * t70 + t364;
t252 = t4 * qJD(1) + t8 * qJD(2);
t27 = m(6) * (t70 * t180 + t71 * t181 + t201 * t223) + m(5) * (t218 + (-t159 * t226 + t160 * t229) * t224) + m(4) * (t222 * qJ(3) + t218) + t241;
t232 = ((qJ(3) + t215) * t222 + t278) * t350 + ((t183 + t201) * t223 + (t62 + t71) * t181 + (t61 + t70) * t180) * t347 + t241;
t235 = (t229 * t186 + t226 * t187) * t348 + (t202 * t109 - t203 * t110) * t347 + t339 * t350;
t9 = -t232 + t235 - m(5) * t356 / 0.2e1;
t251 = qJD(1) * t9 - qJD(2) * t27;
t237 = ((t61 - t66) * t203 + (t62 + t65) * t202) * t347 + t257 + t357;
t12 = (t317 / 0.2e1 + t318 / 0.2e1) * t224 - t343 / 0.2e1 + t237 - t279;
t25 = t202 * t264 + t279 - t352;
t15 = t25 + t353 + t370;
t250 = -t12 * qJD(1) + t15 * qJD(2);
t21 = t257 - t279;
t249 = -t21 * qJD(1) + t25 * qJD(2);
t247 = t256 + t113;
t243 = -Ifges(5,5) * t295 - Ifges(5,6) * t294 + t113 - t154;
t18 = t154 / 0.2e1 + (-t66 / 0.2e1 + t345) * mrSges(6,2) + (t65 / 0.2e1 + t62 / 0.2e1) * mrSges(6,1) + (t286 / 0.2e1 + (-t149 / 0.2e1 + t264) * t228) * pkin(4);
t204 = (mrSges(6,1) * t225 + mrSges(6,2) * t228) * pkin(4);
t24 = (-t70 / 0.2e1 + t81 / 0.2e1) * mrSges(6,2) + (t344 - t80 / 0.2e1) * mrSges(6,1) + t256;
t240 = t18 * qJD(1) - t24 * qJD(2) + t204 * qJD(4);
t211 = qJ(3) * t277;
t198 = t204 * qJD(5);
t22 = -t332 / 0.2e1 - t333 / 0.2e1 + t331 / 0.2e1 - t330 / 0.2e1 + t247;
t16 = t355 + t353 - t370;
t14 = -t336 / 0.2e1 - t337 / 0.2e1 + t335 / 0.2e1 - t334 / 0.2e1 + t247;
t11 = t237 + t279 + t353;
t10 = t356 * t348 + t232 + t235;
t3 = t233 + t359;
t1 = t231 - t375;
t17 = [qJD(2) * t13 + qJD(3) * t23 + qJD(4) * t5 + qJD(5) * t7, t309 + t242 * qJD(2) + t10 * qJD(3) + t1 * qJD(4) + t3 * qJD(5) + qJD(2) * t236 + 0.2e1 * ((t70 * t109 + t71 * t110 + t201 * t276) * t347 + (t159 * t186 + t160 * t187 + t211) * t348 + (t211 + (-pkin(2) * t227 + qJ(3) * t298) * pkin(1)) * t350) * qJD(2), qJD(2) * t10 + qJD(4) * t11 + t310 + t361, t316 + t1 * qJD(2) + t11 * qJD(3) + (-t136 * mrSges(5,1) - t135 * mrSges(5,2) + t243 - t334 + t335) * qJD(4) + t14 * qJD(5) + (m(6) * (t225 * t66 + t228 * t65) - t271) * t325, t312 + t3 * qJD(2) + t363 + t14 * qJD(4) + (t113 - t336 - t337) * qJD(5); -qJD(3) * t9 + qJD(4) * t2 + qJD(5) * t4 - t309, qJD(3) * t27 - qJD(4) * t6 + qJD(5) * t8, qJD(4) * t16 - t251 + t361, t16 * qJD(3) + (-t160 * mrSges(5,1) - t159 * mrSges(5,2) + t243 - t330 + t331) * qJD(4) + t22 * qJD(5) + (-t271 + m(6) * (t225 * t81 + t228 * t80)) * t325 + t253, t363 + t22 * qJD(4) + (t113 - t332 - t333) * qJD(5) + t252; qJD(2) * t9 + qJD(4) * t12 + qJD(5) * t21 - t310, -qJD(4) * t15 - qJD(5) * t25 + t251, 0, (m(6) * (t202 * t341 + t203 * t340) - t254 + t259) * qJD(4) + t371 - t250, qJD(4) * t259 - t249 + t371; -qJD(2) * t2 - qJD(3) * t12 - qJD(5) * t18 - t316, qJD(3) * t15 + qJD(5) * t24 - t253, t250, -t198, -t198 - t240; -qJD(2) * t4 - qJD(3) * t21 + qJD(4) * t18 - t312, qJD(3) * t25 - qJD(4) * t24 - t252, t249, t240, 0;];
Cq = t17;
