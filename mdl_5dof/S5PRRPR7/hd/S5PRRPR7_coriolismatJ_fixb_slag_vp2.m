% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRPR7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR7_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR7_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR7_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR7_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:35:10
% EndTime: 2019-12-05 16:35:24
% DurationCPUTime: 4.54s
% Computational Cost: add. (5862->433), mult. (14842->661), div. (0->0), fcn. (15021->10), ass. (0->247)
t221 = sin(pkin(10));
t223 = cos(pkin(10));
t201 = -mrSges(5,1) * t223 + mrSges(5,2) * t221;
t379 = -mrSges(4,1) + t201;
t365 = m(6) / 0.2e1;
t378 = 0.2e1 * t365;
t368 = m(5) / 0.2e1;
t377 = 0.2e1 * t368;
t376 = Ifges(6,3) / 0.2e1;
t228 = cos(qJ(3));
t375 = -t228 / 0.2e1;
t300 = t221 * t228;
t374 = t300 * t376;
t227 = cos(qJ(5));
t282 = t227 * t228;
t224 = sin(qJ(5));
t225 = sin(qJ(3));
t290 = t224 * t225;
t175 = -t223 * t290 - t282;
t287 = t225 * t227;
t289 = t224 * t228;
t177 = t223 * t287 - t289;
t109 = -mrSges(6,1) * t175 + mrSges(6,2) * t177;
t296 = t223 * t225;
t319 = t228 * mrSges(5,1);
t192 = -mrSges(5,3) * t296 - t319;
t373 = t109 - t192;
t303 = t221 * t224;
t278 = mrSges(6,3) * t303;
t188 = mrSges(6,2) * t223 - t278;
t283 = t227 * t188;
t301 = t221 * t227;
t277 = mrSges(6,3) * t301;
t189 = -t223 * mrSges(6,1) - t277;
t291 = t224 * t189;
t246 = t291 / 0.2e1 - t283 / 0.2e1;
t302 = t221 * t225;
t329 = t175 * mrSges(6,3);
t132 = -mrSges(6,2) * t302 + t329;
t284 = t227 * t132;
t328 = t177 * mrSges(6,3);
t134 = mrSges(6,1) * t302 - t328;
t292 = t224 * t134;
t247 = -t292 / 0.2e1 + t284 / 0.2e1;
t222 = sin(pkin(5));
t226 = sin(qJ(2));
t229 = cos(qJ(2));
t293 = t223 * t229;
t145 = (t221 * t226 + t228 * t293) * t222;
t310 = t145 * t223;
t299 = t222 * t229;
t276 = t221 * t299;
t285 = t226 * t222;
t144 = -t223 * t285 + t228 * t276;
t311 = t144 * t221;
t372 = t310 + t311;
t318 = t228 * mrSges(4,2);
t204 = t225 * mrSges(4,1) + t318;
t176 = -t223 * t289 + t287;
t178 = t223 * t282 + t290;
t371 = -Ifges(6,6) * t176 / 0.2e1 - Ifges(6,5) * t178 / 0.2e1;
t321 = t227 * mrSges(6,2);
t325 = t224 * mrSges(6,1);
t261 = t321 + t325;
t182 = t261 * t221;
t217 = t221 ^ 2;
t218 = t223 ^ 2;
t370 = (t218 + t217) * mrSges(5,3) + t221 * t182;
t369 = -m(5) / 0.2e1;
t367 = m(5) / 0.4e1;
t366 = -m(6) / 0.2e1;
t364 = mrSges(6,1) / 0.2e1;
t363 = -mrSges(6,2) / 0.2e1;
t316 = cos(pkin(5));
t174 = t225 * t316 + t228 * t285;
t127 = t174 * t223 - t276;
t173 = t225 * t285 - t228 * t316;
t68 = -t127 * t224 + t173 * t227;
t361 = t68 / 0.2e1;
t126 = t174 * t221 + t222 * t293;
t360 = -t126 / 0.2e1;
t359 = -t132 / 0.2e1;
t358 = t132 / 0.2e1;
t133 = -mrSges(6,2) * t300 + mrSges(6,3) * t176;
t357 = t133 / 0.2e1;
t356 = t134 / 0.2e1;
t355 = t174 / 0.2e1;
t354 = t175 / 0.2e1;
t353 = t176 / 0.2e1;
t352 = t177 / 0.2e1;
t351 = t178 / 0.2e1;
t350 = -t188 / 0.2e1;
t349 = -t189 / 0.2e1;
t348 = -t221 / 0.2e1;
t347 = t221 / 0.2e1;
t346 = t223 / 0.2e1;
t344 = t227 / 0.2e1;
t200 = -pkin(3) * t228 - qJ(4) * t225 - pkin(2);
t294 = t223 * t228;
t147 = pkin(7) * t294 + t200 * t221;
t131 = -pkin(8) * t228 + t147;
t252 = pkin(4) * t221 - pkin(8) * t223 + pkin(7);
t161 = t252 * t225;
t70 = -t131 * t224 + t161 * t227;
t343 = t70 * mrSges(6,3);
t71 = t131 * t227 + t161 * t224;
t342 = t71 * mrSges(6,3);
t341 = Ifges(5,4) * t221;
t340 = Ifges(5,4) * t223;
t339 = Ifges(6,4) * t177;
t338 = Ifges(6,4) * t224;
t337 = Ifges(6,4) * t227;
t336 = Ifges(5,5) * t223;
t335 = Ifges(5,5) * t225;
t199 = -pkin(4) * t223 - pkin(8) * t221 - pkin(3);
t297 = t223 * t224;
t142 = -qJ(4) * t297 + t199 * t227;
t331 = t142 * mrSges(6,3);
t295 = t223 * t227;
t143 = qJ(4) * t295 + t199 * t224;
t330 = t143 * mrSges(6,3);
t327 = t221 * mrSges(5,1);
t326 = t223 * mrSges(5,2);
t324 = t224 * t68;
t322 = t225 * mrSges(4,2);
t69 = t127 * t227 + t173 * t224;
t320 = t227 * t69;
t317 = t228 * mrSges(5,2);
t315 = qJ(4) * t218;
t314 = t126 * t144;
t313 = t126 * t223;
t305 = t221 * t126;
t255 = t127 * t223 + t305;
t94 = t173 * t297 + t174 * t227;
t95 = -t173 * t295 + t174 * t224;
t14 = 0.4e1 * (t174 - t255) * t367 * t173 + (-t173 * t305 + t68 * t94 + t69 * t95) * m(6);
t312 = t14 * qJD(1);
t286 = t225 * t229;
t275 = t222 * t286;
t106 = -t145 * t224 + t227 * t275;
t107 = t145 * t227 + t224 * t275;
t308 = t173 * t225;
t15 = m(5) * (t127 * t145 + t314) + m(6) * (t106 * t68 + t107 * t69 + t314) + 0.4e1 * (t308 * t367 + m(4) * (t174 * t228 + t308) / 0.4e1 - m(4) * t285 / 0.4e1) * t299;
t309 = t15 * qJD(1);
t307 = t200 * t223;
t203 = pkin(3) * t225 - qJ(4) * t228;
t306 = t203 * t223;
t215 = t217 * qJ(4);
t190 = -mrSges(5,3) * t302 + t317;
t298 = t223 * t190;
t288 = t225 * t182;
t111 = Ifges(6,5) * t175 - Ifges(6,6) * t177;
t280 = pkin(7) * t368;
t274 = -t329 / 0.2e1;
t273 = -t328 / 0.2e1;
t272 = pkin(7) * t221 + pkin(4);
t271 = t302 / 0.2e1;
t269 = t299 / 0.2e1;
t264 = t192 / 0.2e1 - t109 / 0.2e1;
t263 = t225 * t269;
t262 = t326 + t327;
t110 = -mrSges(6,1) * t176 + mrSges(6,2) * t178;
t130 = t228 * t272 - t307;
t135 = mrSges(6,1) * t300 - mrSges(6,3) * t178;
t136 = -t225 * t272 - t306;
t152 = pkin(7) * t302 + t306;
t187 = t221 * t203;
t153 = -pkin(7) * t296 + t187;
t180 = t262 * t225;
t181 = t262 * t228;
t191 = -mrSges(5,2) * t225 - mrSges(5,3) * t300;
t193 = mrSges(5,1) * t225 - mrSges(5,3) * t294;
t146 = -pkin(7) * t300 + t307;
t253 = -t146 * t221 + t147 * t223;
t137 = t187 + (-pkin(7) * t223 + pkin(8)) * t225;
t162 = t252 * t228;
t74 = -t137 * t224 + t162 * t227;
t75 = t137 * t227 + t162 * t224;
t230 = (-t193 / 0.2e1 + t110 / 0.2e1) * t126 + (t181 / 0.2e1 - t298 / 0.2e1 + t264 * t221) * t173 + (pkin(7) * t174 * t225 - t126 * t152 + t127 * t153 + (pkin(7) * t228 - t253) * t173) * t368 + (-t130 * t173 * t221 + t126 * t136 + t68 * t74 + t69 * t75 + t70 * t94 + t71 * t95) * t365 + t127 * t191 / 0.2e1 + t180 * t355 + t135 * t361 + t69 * t357 + t94 * t356 + t95 * t358;
t232 = (-pkin(3) * t275 + qJ(4) * t372) * t369 + (qJ(4) * t311 + t106 * t142 + t107 * t143) * t366 + t106 * t349 + t107 * t350 - t144 * t182 / 0.2e1;
t2 = (t318 / 0.2e1 - t204 / 0.2e1 + (mrSges(4,1) / 0.2e1 - t201 / 0.2e1) * t225) * t299 + t230 + (-t310 / 0.2e1 - t311 / 0.2e1) * mrSges(5,3) + t232;
t160 = t335 + (Ifges(5,1) * t223 - t341) * t228;
t248 = -Ifges(5,6) * t225 + (-t221 * Ifges(5,2) + t340) * t375 + t374 - t371;
t90 = Ifges(6,2) * t175 + Ifges(6,6) * t302 + t339;
t91 = Ifges(6,4) * t178 + Ifges(6,2) * t176 + Ifges(6,6) * t300;
t171 = Ifges(6,4) * t175;
t92 = Ifges(6,1) * t177 + Ifges(6,5) * t302 + t171;
t93 = Ifges(6,1) * t178 + Ifges(6,4) * t176 + Ifges(6,5) * t300;
t3 = -pkin(2) * t204 + t153 * t190 + t147 * t191 + t152 * t192 + t146 * t193 + t92 * t351 + t91 * t354 + t90 * t353 + t93 * t352 + t74 * t134 + t70 * t135 + t136 * t109 + t130 * t110 + t75 * t132 + t71 * t133 + m(5) * (t146 * t152 + t147 * t153) + m(6) * (t130 * t136 + t70 * t74 + t71 * t75) + (t160 * t346 + pkin(7) * t181 + (-Ifges(4,4) + t336 / 0.2e1) * t225 + t248 * t221) * t225 + ((Ifges(6,5) * t177 + Ifges(6,6) * t175) * t347 + pkin(7) * t180 + (-Ifges(4,2) + Ifges(4,1) - Ifges(5,3) + t218 * Ifges(5,1) / 0.2e1 + m(5) * pkin(7) ^ 2 + (-t340 + (Ifges(5,2) / 0.2e1 + t376) * t221) * t221) * t225 + (Ifges(5,6) * t221 + Ifges(4,4) - t336) * t228) * t228;
t260 = t2 * qJD(1) + t3 * qJD(2);
t108 = mrSges(6,1) * t177 + mrSges(6,2) * t175;
t240 = t108 * t360 + t356 * t69 + t359 * t68;
t250 = t106 * t364 + t107 * t363;
t6 = (t352 * t69 + t354 * t68) * mrSges(6,3) + t240 + t250;
t112 = -Ifges(6,2) * t177 + t171;
t113 = Ifges(6,1) * t175 - t339;
t8 = t111 * t271 + t130 * t108 - t71 * t134 + t70 * t132 + (-t342 + t113 / 0.2e1 - t90 / 0.2e1) * t177 + (-t343 + t92 / 0.2e1 + t112 / 0.2e1) * t175;
t259 = -qJD(1) * t6 + qJD(2) * t8;
t258 = -t320 + t324;
t17 = (t373 * t223 + (-t190 - t284 + t292) * t221 + m(6) * (t130 * t223 - t301 * t71 + t303 * t70) + m(5) * (-t146 * t223 - t147 * t221)) * t225;
t242 = (t106 * t227 + t107 * t224) * t365;
t245 = m(5) * (-t127 * t221 + t313);
t21 = t242 + (m(5) * t269 - t245 / 0.2e1 + (-t301 * t69 + t303 * t68 + t313) * t366) * t225;
t257 = qJD(1) * t21 - qJD(2) * t17;
t25 = (mrSges(6,2) * t271 + t359 + t329 / 0.2e1) * t227 + (mrSges(6,1) * t271 + t328 / 0.2e1 + t356) * t224;
t249 = -t321 / 0.2e1 - t325 / 0.2e1;
t243 = t249 * t223;
t38 = t243 + t246;
t256 = qJD(2) * t25 + qJD(3) * t38;
t254 = -t142 * t224 + t143 * t227;
t251 = t363 * t95 + t364 * t94;
t183 = (-Ifges(6,5) * t224 - Ifges(6,6) * t227) * t221;
t157 = -Ifges(6,6) * t223 + (-Ifges(6,2) * t224 + t337) * t221;
t158 = -Ifges(6,5) * t223 + (Ifges(6,1) * t227 - t338) * t221;
t179 = (mrSges(6,1) * t227 - mrSges(6,2) * t224) * t221;
t184 = (-Ifges(6,2) * t227 - t338) * t221;
t185 = (-Ifges(6,1) * t224 - t337) * t221;
t18 = -t223 * t183 / 0.2e1 + t142 * t188 - t143 * t189 + (qJ(4) * t179 + (t185 / 0.2e1 - t157 / 0.2e1 - t330) * t227 + (-t158 / 0.2e1 - t184 / 0.2e1 + t331) * t224) * t221;
t231 = (-t157 / 0.4e1 + t185 / 0.4e1 - t330 / 0.2e1) * t177 + (t184 / 0.4e1 + t158 / 0.4e1 - t331 / 0.2e1) * t175 + t130 * t179 / 0.2e1 + t142 * t358 - t143 * t134 / 0.2e1 - t223 * t111 / 0.4e1 + t70 * t188 / 0.2e1 + t71 * t349;
t234 = t225 * t183 / 0.4e1 + qJ(4) * t108 / 0.2e1 + (-t342 / 0.2e1 + t113 / 0.4e1 - t90 / 0.4e1) * t227 + (t343 / 0.2e1 - t112 / 0.4e1 - t92 / 0.4e1) * t224;
t237 = -t74 * mrSges(6,1) / 0.2e1 + t75 * mrSges(6,2) / 0.2e1 + t371;
t5 = (Ifges(6,3) * t375 + t234) * t221 + t231 + t237;
t239 = t179 * t360 + t68 * t350 + t69 * t189 / 0.2e1;
t9 = (t320 / 0.2e1 - t324 / 0.2e1) * t221 * mrSges(6,3) + t239 + t251;
t244 = -t9 * qJD(1) + t5 * qJD(2) + t18 * qJD(3);
t233 = t253 * t369 + ((-t224 * t70 + t227 * t71) * t223 + (t130 + (qJ(4) * t223 - t254) * t225) * t221) * t366;
t236 = (t224 * t75 + t227 * t74) * t365 + t224 * t357 + t135 * t344;
t11 = t228 * t280 + (t317 / 0.2e1 - t190 / 0.2e1 - t288 / 0.2e1 - t247) * t223 + (t319 / 0.2e1 - t246 * t225 + t264) * t221 + t233 + t236;
t235 = t255 * t369 + (-t223 * t258 + t305) * t366;
t238 = m(5) * t355 + (t224 * t95 + t227 * t94) * t365;
t19 = t235 + t238;
t33 = (t283 - t291) * t223 + m(6) * (t223 * t254 + t215) + m(5) * (t215 + t315) + t370;
t241 = qJD(1) * t19 + qJD(2) * t11 - qJD(3) * t33;
t220 = t228 ^ 2;
t219 = t225 ^ 2;
t202 = t219 * pkin(7) * t299;
t155 = t173 * t215;
t39 = t243 - t246;
t26 = t224 * t273 + t227 * t274 - t249 * t302 + t247;
t22 = m(5) * t263 + t242 + (t245 + m(6) * (t221 * t258 + t313)) * t225 / 0.2e1;
t20 = -t235 + t238;
t12 = t109 * t347 + t298 / 0.2e1 + t192 * t348 + t288 * t346 + (t280 + t326 / 0.2e1 + t327 / 0.2e1) * t228 - t233 + t236 + t246 * t302 + t247 * t223;
t10 = -t69 * t277 / 0.2e1 + t278 * t361 - t239 + t251;
t7 = t273 * t69 + t274 * t68 - t240 + t250;
t4 = t221 * t234 + t231 - t237 + t374;
t1 = t201 * t263 + t230 - t232 + t372 * mrSges(5,3) / 0.2e1 - t204 * t299;
t13 = [qJD(2) * t15 + qJD(3) * t14, t1 * qJD(3) + t22 * qJD(4) + t7 * qJD(5) + t309 + (t106 * t134 + t107 * t132 + t373 * t144 + t145 * t190 + (-t144 * t146 + t145 * t147 + t202) * t377 + (t106 * t70 + t107 * t71 + t130 * t144) * t378 + m(4) * t202 + (t180 * t286 + (-mrSges(3,2) + m(4) * pkin(7) * t220 + (t219 + t220) * mrSges(4,3)) * t229 + (-m(4) * pkin(2) - t228 * mrSges(4,1) - mrSges(3,1) + t322) * t226) * t222) * qJD(2), t1 * qJD(2) + t20 * qJD(4) + t10 * qJD(5) + t312 + (t95 * t188 + t94 * t189 + (mrSges(4,2) - t370) * t173 + (-t173 * t315 - t155) * t377 + (t142 * t94 + t143 * t95 - t155) * t378 + (-pkin(3) * t377 + t379) * t174) * qJD(3), qJD(2) * t22 + qJD(3) * t20, t7 * qJD(2) + t10 * qJD(3) + (-mrSges(6,1) * t69 - mrSges(6,2) * t68) * qJD(5); qJD(3) * t2 - qJD(4) * t21 - qJD(5) * t6 - t309, qJD(3) * t3 + qJD(4) * t17 + qJD(5) * t8, t12 * qJD(4) + t4 * qJD(5) + t260 + (-Ifges(4,6) * t225 + t75 * t188 + t74 * t189 + t158 * t351 - pkin(3) * t181 + t136 * t182 + t157 * t353 + t143 * t133 + t142 * t135 + pkin(7) * t322 + m(6) * (t142 * t74 + t143 * t75) + (t153 * mrSges(5,3) - t248) * t223 + (t335 / 0.2e1 + t160 / 0.2e1 + t93 * t344 - t224 * t91 / 0.2e1 - t152 * mrSges(5,3)) * t221 + ((m(5) * t153 + t191) * t223 + (-m(5) * t152 + m(6) * t136 + t110 - t193) * t221) * qJ(4) + (Ifges(4,5) + (-Ifges(6,3) * t223 + (Ifges(6,5) * t227 - Ifges(6,6) * t224) * t221) * t347 + (Ifges(5,1) * t221 + t340) * t346 + (Ifges(5,2) * t223 + t341) * t348 + (-m(5) * pkin(3) + t379) * pkin(7)) * t228) * qJD(3), qJD(3) * t12 + qJD(5) * t26 - t257, t4 * qJD(3) + t26 * qJD(4) + (-mrSges(6,1) * t71 - mrSges(6,2) * t70 + t111) * qJD(5) + t259; -qJD(2) * t2 - qJD(4) * t19 - qJD(5) * t9 - t312, -qJD(4) * t11 + qJD(5) * t5 - t260, qJD(4) * t33 + qJD(5) * t18, qJD(5) * t39 - t241, t39 * qJD(4) + (-mrSges(6,1) * t143 - mrSges(6,2) * t142 + t183) * qJD(5) + t244; qJD(2) * t21 + qJD(3) * t19, qJD(3) * t11 - qJD(5) * t25 + t257, -qJD(5) * t38 + t241, 0, -qJD(5) * t261 - t256; t6 * qJD(2) + t9 * qJD(3), -qJD(3) * t5 + qJD(4) * t25 - t259, qJD(4) * t38 - t244, t256, 0;];
Cq = t13;
