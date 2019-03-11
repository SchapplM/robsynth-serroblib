% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
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
% Datum: 2019-03-09 02:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPPR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:38:05
% EndTime: 2019-03-09 02:38:14
% DurationCPUTime: 5.18s
% Computational Cost: add. (14611->390), mult. (28365->553), div. (0->0), fcn. (31916->10), ass. (0->219)
t236 = sin(pkin(11));
t237 = cos(pkin(11));
t345 = sin(qJ(6));
t347 = cos(qJ(6));
t220 = -t236 * t347 - t237 * t345;
t318 = sin(pkin(10));
t319 = cos(pkin(10));
t346 = sin(qJ(3));
t348 = cos(qJ(3));
t218 = -t318 * t348 - t319 * t346;
t246 = t220 * t218;
t373 = t345 * t236 - t347 * t237;
t379 = t373 * t218;
t383 = m(7) * (t220 * t246 + t373 * t379);
t385 = -t383 / 0.2e1;
t384 = t383 / 0.2e1;
t360 = t379 / 0.2e1;
t381 = -mrSges(7,1) / 0.2e1;
t216 = t318 * t346 - t319 * t348;
t161 = t220 * t216;
t164 = t373 * t216;
t65 = t161 * t373 - t164 * t220;
t380 = m(7) * t65 * qJD(4);
t310 = t216 * t236;
t231 = sin(pkin(9)) * pkin(1) + pkin(7);
t282 = t346 * t231;
t214 = -qJ(4) * t346 - t282;
t285 = t348 * t231;
t215 = qJ(4) * t348 + t285;
t374 = t318 * t214 + t319 * t215;
t125 = -pkin(5) * t310 + t374;
t332 = t164 * mrSges(7,2);
t334 = t161 * mrSges(7,1);
t262 = t334 / 0.2e1 + t332 / 0.2e1;
t370 = -m(6) / 0.2e1;
t378 = t374 * t370 - m(7) * t125 / 0.2e1 - t262;
t363 = -t246 / 0.2e1;
t234 = t236 ^ 2;
t235 = t237 ^ 2;
t295 = t234 + t235;
t377 = m(6) * t295;
t280 = t318 * pkin(3);
t283 = t319 * pkin(3);
t223 = -mrSges(6,1) * t237 + mrSges(6,2) * t236;
t376 = mrSges(5,1) - t223;
t177 = -t319 * t214 + t215 * t318;
t289 = t346 * pkin(3);
t183 = -t218 * pkin(4) + t216 * qJ(5) + t289;
t98 = t177 * t236 + t237 * t183;
t99 = -t177 * t237 + t236 * t183;
t267 = -t236 * t98 + t237 * t99;
t372 = 0.2e1 * t218;
t371 = m(5) / 0.2e1;
t369 = m(6) / 0.2e1;
t368 = m(7) / 0.2e1;
t365 = -t161 / 0.2e1;
t364 = t246 / 0.2e1;
t362 = t164 / 0.2e1;
t361 = -t379 / 0.2e1;
t224 = t280 + qJ(5);
t344 = pkin(8) + t224;
t206 = t344 * t236;
t207 = t344 * t237;
t175 = -t206 * t345 + t207 * t347;
t359 = -t175 / 0.2e1;
t335 = Ifges(7,4) * t220;
t191 = -Ifges(7,2) * t373 - t335;
t358 = t191 / 0.2e1;
t212 = Ifges(7,4) * t373;
t193 = -Ifges(7,1) * t220 - t212;
t357 = t193 / 0.2e1;
t356 = t216 / 0.2e1;
t354 = t218 / 0.2e1;
t353 = -t373 / 0.2e1;
t352 = -t220 / 0.2e1;
t351 = t220 / 0.2e1;
t350 = t236 / 0.2e1;
t349 = t237 / 0.2e1;
t290 = t65 * t368;
t343 = qJD(3) * t290;
t341 = mrSges(5,3) * t216;
t340 = mrSges(5,3) * t218;
t339 = mrSges(6,3) * t216;
t338 = Ifges(6,4) * t236;
t337 = Ifges(6,4) * t237;
t336 = Ifges(7,4) * t379;
t331 = t379 * mrSges(7,1);
t330 = t373 * mrSges(7,3);
t329 = t220 * mrSges(7,3);
t328 = t236 * mrSges(6,1);
t327 = t236 * Ifges(6,2);
t326 = t236 * Ifges(6,6);
t324 = t237 * mrSges(6,2);
t323 = t237 * Ifges(6,5);
t114 = -mrSges(7,2) * t216 - mrSges(7,3) * t246;
t116 = mrSges(7,1) * t216 - mrSges(7,3) * t379;
t90 = -t246 * mrSges(7,2) + t331;
t7 = t114 * t363 + t116 * t361 + t90 * t356 + (-t379 ^ 2 / 0.2e1 - t246 ^ 2 / 0.2e1) * mrSges(7,3);
t321 = t7 * qJD(1);
t194 = t216 * t218;
t252 = m(7) * (t161 * t246 + t164 * t379 - t194);
t274 = t295 * t218;
t253 = m(6) * (t216 * t274 - t194);
t27 = t252 / 0.2e1 + t253 / 0.2e1;
t317 = qJD(1) * t27;
t316 = t161 * t116;
t315 = t161 * t220;
t314 = t164 * t114;
t313 = t164 * t373;
t311 = t177 * t218;
t309 = t216 * t237;
t308 = t218 * t236;
t307 = t218 * t237;
t186 = -mrSges(6,1) * t218 + mrSges(6,3) * t309;
t306 = t236 * t186;
t187 = t216 * mrSges(6,1) + mrSges(6,3) * t307;
t305 = t236 * t187;
t184 = mrSges(6,2) * t218 + mrSges(6,3) * t310;
t304 = t237 * t184;
t185 = -t216 * mrSges(6,2) + mrSges(6,3) * t308;
t303 = t237 * t185;
t293 = -m(7) / 0.4e1 - m(6) / 0.4e1;
t48 = t385 + (-t377 / 0.4e1 + t293) * t372;
t302 = t48 * qJD(1);
t299 = -Ifges(7,5) * t246 - Ifges(7,6) * t379;
t286 = -cos(pkin(9)) * pkin(1) - pkin(2);
t222 = -pkin(3) * t348 + t286;
t159 = t216 * pkin(4) + t218 * qJ(5) + t222;
t89 = t236 * t159 + t237 * t374;
t296 = -Ifges(7,5) * t373 + Ifges(7,6) * t220;
t188 = -t220 * mrSges(7,1) - mrSges(7,2) * t373;
t294 = t188 * qJD(6);
t288 = t379 * t329;
t279 = -t309 / 0.2e1;
t189 = mrSges(7,1) * t373 - mrSges(7,2) * t220;
t278 = -t189 / 0.2e1 - t223 / 0.2e1;
t88 = t237 * t159 - t236 * t374;
t123 = t288 / 0.2e1;
t232 = -t283 - pkin(4);
t272 = t293 * t372;
t271 = -t324 - t328;
t113 = mrSges(7,2) * t218 - t161 * mrSges(7,3);
t115 = -mrSges(7,1) * t218 - t164 * mrSges(7,3);
t126 = -pkin(5) * t308 + t177;
t136 = -Ifges(6,6) * t218 + (t327 - t337) * t216;
t137 = -Ifges(6,5) * t218 + (-t237 * Ifges(6,1) + t338) * t216;
t180 = t271 * t216;
t181 = t271 * t218;
t208 = t216 * mrSges(5,2);
t254 = mrSges(4,1) * t346 + mrSges(4,2) * t348;
t263 = Ifges(7,5) * t362 + Ifges(7,6) * t365;
t58 = pkin(5) * t216 + pkin(8) * t307 + t88;
t66 = pkin(8) * t308 + t89;
t33 = -t345 * t66 + t347 * t58;
t34 = t345 * t58 + t347 * t66;
t64 = -pkin(5) * t218 + pkin(8) * t309 + t98;
t82 = pkin(8) * t310 + t99;
t38 = -t345 * t82 + t347 * t64;
t39 = t345 * t64 + t347 * t82;
t78 = Ifges(7,4) * t164 - Ifges(7,2) * t161 - Ifges(7,6) * t218;
t79 = -Ifges(7,2) * t246 + t216 * Ifges(7,6) + t336;
t80 = Ifges(7,1) * t164 - Ifges(7,4) * t161 - Ifges(7,5) * t218;
t155 = Ifges(7,4) * t246;
t81 = Ifges(7,1) * t379 + t216 * Ifges(7,5) - t155;
t91 = t332 + t334;
t92 = mrSges(7,1) * t246 + mrSges(7,2) * t379;
t1 = t286 * t254 + m(7) * (t125 * t126 + t33 * t38 + t34 * t39) + (Ifges(4,1) - Ifges(4,2)) * t348 * t346 + (mrSges(5,1) * t289 + (-Ifges(6,3) + Ifges(5,1) - Ifges(5,2) - Ifges(7,3) + t235 * Ifges(6,1) / 0.2e1 + (-t337 + t327 / 0.2e1) * t236) * t218 + t263 + (Ifges(5,4) - t323 + t326) * t216) * t216 + (-t346 ^ 2 + t348 ^ 2) * Ifges(4,4) + (Ifges(7,5) * t361 + Ifges(7,6) * t364 - mrSges(5,2) * t289 + t136 * t350 - t237 * t137 / 0.2e1 + (t323 / 0.2e1 - t326 / 0.2e1 - Ifges(5,4)) * t218) * t218 + t80 * t360 + t81 * t362 + t78 * t363 + t79 * t365 + m(6) * (t177 * t374 + t88 * t98 + t89 * t99) + t374 * t181 + t34 * t113 + t39 * t114 + t33 * t115 + t38 * t116 + t125 * t92 + t126 * t91 + t177 * t180 + t89 * t184 + t99 * t185 + t88 * t186 + t98 * t187 + (m(5) * t289 - t218 * mrSges(5,1) - t208) * t222;
t251 = -t126 * t218 - t161 * t33 + t164 * t34;
t260 = t305 / 0.2e1 - t303 / 0.2e1;
t268 = t236 * t88 - t237 * t89;
t5 = t314 / 0.2e1 + t113 * t360 - t316 / 0.2e1 + t115 * t363 + (-t304 / 0.2e1 + t306 / 0.2e1 - t92 / 0.2e1 - t181 / 0.2e1) * t218 + (t91 / 0.2e1 + t180 / 0.2e1 + t260) * t216 + (t125 * t216 - t246 * t38 + t379 * t39 + t251) * t368 + ((-t177 - t267) * t218 + (t374 + t268) * t216) * t369;
t270 = t1 * qJD(1) + t5 * qJD(2);
t93 = -Ifges(7,2) * t379 - t155;
t94 = -Ifges(7,1) * t246 - t336;
t6 = t126 * t90 + t299 * t356 + t33 * t114 - t34 * t116 + (t94 / 0.2e1 - t79 / 0.2e1 - t34 * mrSges(7,3)) * t379 - (t81 / 0.2e1 + t93 / 0.2e1 - t33 * mrSges(7,3)) * t246;
t269 = t6 * qJD(1) + t7 * qJD(2);
t8 = t314 - t316 + (-t181 - t92 + t340) * t218 + (-t303 + t305 + t341) * t216 + m(7) * t251 + m(6) * (t216 * t268 - t311) + m(5) * (-t216 * t374 - t311);
t266 = -t8 * qJD(1) - t27 * qJD(2);
t54 = 0.2e1 * mrSges(7,1) * t360 + 0.2e1 * t363 * mrSges(7,2);
t265 = qJD(1) * t54 + qJD(3) * t188;
t261 = t114 * t353 + t116 * t351;
t241 = (t236 * t99 + t237 * t98) * t369 + (-t220 * t39 - t373 * t38) * t368 + t115 * t353 + t113 * t352 + t184 * t350 + t186 * t349 + t289 * t371;
t174 = -t206 * t347 - t207 * t345;
t221 = -t237 * pkin(5) + t232;
t248 = (-t216 * t224 * t295 - t218 * t232) * t369 + (-t216 * t318 + t218 * t319) * pkin(3) * t371 + (-t161 * t174 + t164 * t175 - t218 * t221) * t368;
t242 = (-t313 / 0.2e1 - t315 / 0.2e1) * mrSges(7,3) + (-t235 / 0.2e1 - t234 / 0.2e1) * t339 + t248;
t9 = t208 + (mrSges(5,1) + t278) * t218 - t241 + t242;
t259 = -t9 * qJD(1) + qJD(2) * t290;
t243 = -t330 * t364 + t123 + t261;
t13 = t243 + t262;
t256 = t13 * qJD(1);
t18 = m(7) * (-t246 * t34 - t33 * t379) - t246 * t114 - t379 * t116 + (t236 * t185 + t237 * t187 + m(6) * (t236 * t89 + t237 * t88)) * t218;
t255 = t18 * qJD(1);
t240 = (-t246 * t353 - t351 * t379) * mrSges(7,3) + t268 * t370 + (-t174 * t379 - t175 * t246 + t220 * t33 - t34 * t373) * t368 - t260 + t261;
t11 = (t324 / 0.2e1 + t328 / 0.2e1) * t216 + t240 + t378;
t45 = (t220 ^ 2 + t373 ^ 2) * mrSges(7,3) + t295 * mrSges(6,3) + m(7) * (t174 * t220 - t175 * t373) + t224 * t377;
t47 = t384 + (t377 / 0.4e1 + t293) * t372;
t250 = t11 * qJD(1) - t47 * qJD(2) + t45 * qJD(3);
t249 = t188 * t356 + t123 - t288 / 0.2e1;
t26 = t252 + t253;
t247 = -t5 * qJD(1) - t26 * qJD(2) - t380 / 0.2e1;
t21 = t249 + t262;
t190 = Ifges(7,2) * t220 - t212;
t192 = -Ifges(7,1) * t373 + t335;
t24 = t221 * t188 + (-t192 / 0.2e1 + t358) * t220 - (t357 + t190 / 0.2e1) * t373;
t239 = -(t81 / 0.4e1 + t93 / 0.4e1) * t373 + (-t94 / 0.4e1 + t79 / 0.4e1) * t220 - (-t174 * mrSges(7,3) / 0.2e1 + t193 / 0.4e1 + t190 / 0.4e1) * t246 + (mrSges(7,3) * t359 + t192 / 0.4e1 - t191 / 0.4e1) * t379 + t126 * t188 / 0.2e1 + t174 * t114 / 0.2e1 + t116 * t359 + t216 * t296 / 0.4e1 + t221 * t90 / 0.2e1;
t244 = Ifges(7,3) * t354 + t38 * t381 + t39 * mrSges(7,2) / 0.2e1 - t263;
t4 = t239 + t244;
t245 = -t4 * qJD(1) - t21 * qJD(2) - t24 * qJD(3);
t55 = t331 / 0.2e1 + t379 * t381;
t50 = t274 * t369 + t272 + t384;
t49 = -t377 * t354 + t272 + t385;
t22 = t249 - t262;
t14 = t243 - t262;
t12 = t218 * t278 + t241 + t242;
t10 = mrSges(6,2) * t279 - mrSges(6,1) * t310 / 0.2e1 + t240 - t378;
t3 = t239 - t244;
t2 = t5 * qJD(3) + t27 * qJD(4) + t7 * qJD(6);
t15 = [qJD(3) * t1 + qJD(4) * t8 + qJD(5) * t18 + qJD(6) * t6, t2 (mrSges(4,2) * t282 - (Ifges(6,5) * t236 - Ifges(7,5) * t220 + Ifges(6,6) * t237 - Ifges(7,6) * t373) * t218 / 0.2e1 - mrSges(4,1) * t285 + (m(6) * t267 + t304 - t306) * t224 + m(7) * (t125 * t221 + t174 * t38 + t175 * t39) + t267 * mrSges(6,3) + t232 * t180 + t221 * t91 - Ifges(5,5) * t216 + Ifges(5,6) * t218 + t38 * t329 - t39 * t330 + t280 * t340 + t283 * t341 + t136 * t349 + t137 * t350 + t80 * t352 + t78 * t353 + t164 * t357 - t161 * t358 + (-m(5) * t283 + m(6) * t232 - t376) * t374 - (m(5) * t280 - mrSges(5,2)) * t177 + (Ifges(6,1) * t236 + t337) * t279 + (Ifges(6,2) * t237 + t338) * t310 / 0.2e1 - Ifges(4,6) * t346 + Ifges(4,5) * t348 + t174 * t115 + t175 * t113 + t125 * t189) * qJD(3) + t12 * qJD(4) + t10 * qJD(5) + t3 * qJD(6) + t270, t12 * qJD(3) + t50 * qJD(5) + t14 * qJD(6) - t266 + t380, t10 * qJD(3) + t50 * qJD(4) + t55 * qJD(6) + t255, t3 * qJD(3) + t14 * qJD(4) + t55 * qJD(5) + (-mrSges(7,1) * t34 - mrSges(7,2) * t33 + t299) * qJD(6) + t269; t2, t26 * qJD(3), t49 * qJD(5) + t22 * qJD(6) - t247 + (t208 - t254 + (-t189 + t376) * t218 + 0.2e1 * t248 - t295 * t339 + (-t313 - t315) * mrSges(7,3)) * qJD(3), t317 + t343, t49 * qJD(3), t22 * qJD(3) - t90 * qJD(6) + t321; qJD(4) * t9 + qJD(5) * t11 + qJD(6) * t4 - t270, -t47 * qJD(5) + t21 * qJD(6) + t247, qJD(5) * t45 + qJD(6) * t24, -t259, t250 (-mrSges(7,1) * t175 - mrSges(7,2) * t174 + t296) * qJD(6) - t245; -qJD(3) * t9 - qJD(5) * t48 + qJD(6) * t13 + t266, -t317 + t343, t259, 0, -t302, t256 - t294; -t11 * qJD(3) + t48 * qJD(4) + t54 * qJD(6) - t255, t47 * qJD(3), -t250 + t294, t302, 0, t265; -qJD(3) * t4 - qJD(4) * t13 - qJD(5) * t54 - t269, -t21 * qJD(3) - t321, -t188 * qJD(5) + t245, -t256, -t265, 0;];
Cq  = t15;
