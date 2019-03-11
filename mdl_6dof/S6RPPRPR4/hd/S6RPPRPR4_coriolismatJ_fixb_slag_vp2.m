% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 01:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPRPR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR4_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR4_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:46:19
% EndTime: 2019-03-09 01:46:25
% DurationCPUTime: 3.51s
% Computational Cost: add. (8725->329), mult. (15883->470), div. (0->0), fcn. (16630->8), ass. (0->208)
t206 = sin(qJ(6));
t198 = t206 ^ 2;
t208 = cos(qJ(6));
t200 = t208 ^ 2;
t359 = t198 + t200;
t371 = t359 * mrSges(7,3);
t209 = cos(qJ(4));
t203 = sin(pkin(9));
t193 = t203 * qJ(2);
t204 = cos(pkin(9));
t210 = -pkin(1) - pkin(2);
t253 = -t204 * t210 + pkin(3) + t193;
t168 = t209 * pkin(4) + t253;
t370 = m(6) * t168;
t369 = t204 * qJ(2);
t202 = sin(pkin(10));
t331 = pkin(4) * t202;
t190 = pkin(8) + t331;
t257 = t190 * t359;
t207 = sin(qJ(4));
t297 = cos(pkin(10));
t177 = t202 * t209 + t207 * t297;
t259 = t297 * t209;
t228 = -t202 * t207 + t259;
t330 = t207 * pkin(4);
t125 = -pkin(5) * t177 + pkin(8) * t228 - t330;
t236 = t203 * t210 - pkin(7) + t369;
t221 = t207 * (qJ(5) - t236);
t172 = t209 * t236;
t256 = qJ(5) * t209 - t172;
t99 = t202 * t256 + t297 * t221;
t54 = t125 * t208 - t206 * t99;
t55 = t125 * t206 + t208 * t99;
t245 = -t54 * t206 + t55 * t208;
t348 = m(6) * pkin(4);
t368 = t202 * t348 - mrSges(6,2);
t179 = -t208 * mrSges(7,1) + t206 * mrSges(7,2);
t269 = t297 * pkin(4);
t191 = -t269 - pkin(5);
t367 = m(7) * t191 - t297 * t348 - mrSges(6,1) + t179;
t115 = t177 * t179;
t366 = -t115 / 0.2e1;
t365 = m(7) * t228;
t317 = Ifges(7,6) * t206;
t320 = Ifges(7,5) * t208;
t235 = -t317 / 0.2e1 + t320 / 0.2e1;
t364 = Ifges(6,4) - t235;
t306 = t208 * mrSges(7,2);
t310 = t206 * mrSges(7,1);
t180 = t306 + t310;
t117 = t180 * t177;
t315 = t177 * mrSges(6,3);
t363 = t117 + t315;
t136 = -t177 * mrSges(6,1) - mrSges(6,2) * t228;
t302 = t209 * mrSges(5,2);
t181 = -t207 * mrSges(5,1) - t302;
t362 = t136 + t181;
t196 = Ifges(7,4) * t208;
t361 = -Ifges(7,2) * t206 + t196;
t360 = Ifges(7,1) * t206 + t196;
t199 = t207 ^ 2;
t358 = t209 ^ 2 + t199;
t234 = t306 / 0.2e1 + t310 / 0.2e1;
t357 = t234 + t180 / 0.2e1;
t322 = Ifges(7,4) * t206;
t182 = Ifges(7,2) * t208 + t322;
t332 = t208 / 0.2e1;
t335 = -t206 / 0.2e1;
t356 = t182 * t335 + t332 * t360;
t285 = t177 * t206;
t129 = -mrSges(7,2) * t228 + mrSges(7,3) * t285;
t277 = t208 * t129;
t284 = t177 * t208;
t131 = mrSges(7,1) * t228 + mrSges(7,3) * t284;
t280 = t206 * t131;
t231 = t280 / 0.2e1 - t277 / 0.2e1;
t333 = -t208 / 0.2e1;
t355 = t129 * t335 + t131 * t333;
t353 = t202 * t221 - t256 * t297;
t352 = -m(7) * t257 - t371;
t351 = -m(6) / 0.2e1;
t350 = m(6) / 0.2e1;
t349 = m(7) / 0.2e1;
t347 = mrSges(7,1) / 0.2e1;
t346 = -mrSges(7,2) / 0.2e1;
t345 = -Ifges(7,3) / 0.2e1;
t279 = t207 * t203;
t156 = t202 * t279 - t203 * t259;
t139 = t156 * t206 - t204 * t208;
t344 = t139 / 0.2e1;
t157 = t177 * t203;
t343 = t157 / 0.2e1;
t342 = -t228 / 0.2e1;
t341 = t177 / 0.2e1;
t340 = -t180 / 0.2e1;
t338 = t190 / 0.2e1;
t336 = -t204 / 0.2e1;
t334 = t206 / 0.2e1;
t140 = -t156 * t208 - t206 * t204;
t292 = t140 * t208;
t293 = t139 * t206;
t241 = t292 - t293;
t232 = -t156 - t241;
t272 = -0.1e1 + t359;
t254 = t272 * t177;
t20 = (-t157 * t254 - t228 * t232) * t349;
t328 = t20 * qJD(4);
t326 = mrSges(6,3) * t228;
t325 = mrSges(7,3) * t177;
t324 = mrSges(7,3) * t190;
t321 = Ifges(7,5) * t228;
t318 = Ifges(7,6) * t228;
t314 = t177 * t99;
t309 = t206 * mrSges(7,3);
t94 = -t177 * t361 + t318;
t308 = t206 * t94;
t307 = t207 * mrSges(5,2);
t305 = t208 * mrSges(7,3);
t252 = Ifges(7,1) * t208 - t322;
t96 = -t177 * t252 + t321;
t304 = t208 * t96;
t303 = t209 * mrSges(5,1);
t158 = t177 * t204;
t299 = t99 * t158;
t290 = t157 * t177;
t214 = (-t228 * t241 - t290) * t349 + (t156 * t228 - t290) * t350;
t219 = -(m(7) * t359 + m(6)) * t203 / 0.2e1;
t26 = t214 + t219;
t296 = qJD(1) * t26;
t286 = t228 * t158;
t159 = t228 * t204;
t289 = t159 * t177;
t35 = (t289 * t359 - t286) * t349 + (-t286 + t289) * t350;
t295 = qJD(1) * t35;
t215 = (t292 / 0.2e1 - t293 / 0.2e1) * t325 + t129 * t344 - t140 * t131 / 0.2e1 + t115 * t343;
t239 = -t159 * t206 + t203 * t208;
t240 = t159 * t208 + t203 * t206;
t217 = t239 * t347 + t240 * t346;
t10 = t215 - t217;
t294 = t10 * qJD(1);
t291 = t157 * t158;
t229 = t177 * t297 - t202 * t228;
t273 = t348 / 0.2e1;
t283 = t191 * t177;
t213 = (-t228 * t257 - t283) * t349 + t366 + t229 * t273 + t342 * t371;
t128 = mrSges(7,2) * t177 + t228 * t309;
t130 = -mrSges(7,1) * t177 + t228 * t305;
t216 = (t206 * t55 + t208 * t54) * t349 + t128 * t334 + t130 * t332 + t330 * t351;
t16 = -t136 + t213 - t216;
t288 = t16 * qJD(1);
t261 = t200 / 0.2e1 + t198 / 0.2e1;
t17 = t115 * t342 + (t261 * t325 + t355) * t177;
t287 = t17 * qJD(1);
t282 = t204 * t207;
t281 = t206 * t130;
t278 = t208 * t128;
t226 = t234 * t228;
t22 = t226 + t231;
t276 = t22 * qJD(1);
t38 = t272 * t341 * t365;
t275 = t38 * qJD(1);
t271 = -Ifges(7,2) / 0.4e1 + Ifges(7,1) / 0.4e1;
t249 = -t317 + t320;
t248 = Ifges(7,5) * t206 + Ifges(7,6) * t208;
t116 = t228 * t180;
t103 = pkin(5) * t228 + t177 * pkin(8) + t168;
t47 = t103 * t208 - t206 * t353;
t48 = t103 * t206 + t208 * t353;
t93 = -Ifges(7,6) * t177 - t228 * t361;
t95 = -Ifges(7,5) * t177 - t228 * t252;
t1 = t99 * t116 - t353 * t117 + t48 * t128 + t55 * t129 + t47 * t130 + t54 * t131 + t168 * t136 - Ifges(5,4) * t199 + t253 * t181 - t330 * t370 + m(7) * (-t353 * t99 + t47 * t54 + t48 * t55) - (t304 / 0.2e1 - t308 / 0.2e1 + mrSges(6,1) * t330 - t364 * t228) * t228 + (Ifges(5,4) * t209 + (Ifges(5,1) - Ifges(5,2)) * t207) * t209 + (t95 * t333 + t93 * t334 + mrSges(6,2) * t330 - t364 * t177 - (Ifges(6,2) - Ifges(6,1) + Ifges(7,3)) * t228) * t177;
t222 = -t116 / 0.2e1 + t231;
t246 = -t206 * t47 + t208 * t48;
t237 = -t246 + t353;
t8 = (-t237 * t349 - t222) * t228 + (-t117 / 0.2e1 + t278 / 0.2e1 - t281 / 0.2e1 + (t245 - t99) * t349) * t177;
t247 = t1 * qJD(1) + t8 * qJD(3);
t5 = -t99 * t115 + t47 * t129 - t48 * t131 + (t228 * t248 / 0.2e1 + t94 * t332 + t96 * t334 + (t182 * t334 + t333 * t360) * t177 + t246 * mrSges(7,3)) * t177;
t244 = t5 * qJD(1) + t17 * qJD(3);
t9 = -m(7) * (t239 * t47 + t240 * t48 - t299) - t240 * t129 - t239 * t131 - m(6) * (t159 * t353 - t299) - mrSges(3,3) + t159 * t326 + t363 * t158 - m(3) * qJ(2) + (-m(4) * t193 - m(5) * t253 - mrSges(6,1) * t228 + mrSges(6,2) * t177 - mrSges(4,1) - t303 + t307 - t370) * t203 + (-mrSges(4,2) + (-m(5) * t236 + mrSges(5,3)) * t358 - m(4) * t369) * t204;
t243 = -t9 * qJD(1) + t35 * qJD(3);
t13 = t363 * t177 - (t277 - t280 - t326) * t228 + m(7) * (-t228 * t246 + t314) + m(6) * (-t228 * t353 + t314);
t242 = -qJD(1) * t13 + qJD(3) * t38;
t238 = t346 * t55 + t347 * t54;
t233 = t191 * t366 - t340 * t99;
t230 = t177 * t248;
t225 = t340 + t234;
t31 = m(7) * t232 * t157;
t211 = -mrSges(5,1) * t282 / 0.2e1 + t302 * t336 + t240 * t305 / 0.2e1 - t239 * t309 / 0.2e1 + (t349 * t257 - mrSges(6,2) / 0.2e1 + t202 * t273) * t159 + (t191 * t349 - mrSges(6,1) / 0.2e1 + t179 / 0.2e1 - t297 * t273) * t158;
t212 = pkin(4) * t282 * t350 + (t139 * t54 + t140 * t55 + t156 * t99 + t157 * t237) * t349 + t130 * t344 + t140 * t128 / 0.2e1 + t156 * t117 / 0.2e1;
t4 = -t116 * t343 + t157 * t231 + t336 * t362 - t211 + t212;
t224 = t4 * qJD(1) + t31 * qJD(2) + t20 * qJD(3);
t39 = t254 * t365;
t223 = t8 * qJD(1) + t20 * qJD(2) + t39 * qJD(3);
t32 = t225 * t157;
t6 = (0.3e1 / 0.4e1 * t318 + t94 / 0.4e1 + t129 * t338) * t206 + (t345 - t261 * t324 + (-t196 / 0.4e1 - t360 / 0.4e1 - t271 * t206) * t206) * t177 + (-0.3e1 / 0.4e1 * t321 - t96 / 0.4e1 + t131 * t338 + (-0.3e1 / 0.4e1 * t322 - t182 / 0.4e1 + t271 * t208) * t177) * t208 + t233 + t238;
t62 = t191 * t180 + t252 * t334 + t332 * t361 + t356;
t71 = t225 * t228;
t218 = t6 * qJD(1) + t32 * qJD(2) - t71 * qJD(3) - t62 * qJD(4);
t72 = t357 * t228;
t33 = t357 * t157;
t25 = t214 - t219;
t23 = t226 - t231;
t18 = t213 + t216;
t11 = t215 + t217;
t7 = t304 / 0.4e1 - t308 / 0.4e1 + t177 * t345 - t233 + t238 + (t182 / 0.2e1 - t252 / 0.4e1) * t284 + t359 * t324 * t341 + t355 * t190 - (-t249 / 0.4e1 + t235) * t228 + (0.2e1 * t360 + t361) * t285 / 0.4e1;
t3 = t211 + t222 * t157 + (-t136 / 0.2e1 - t181 / 0.2e1) * t204 + t212;
t2 = qJD(2) * t35 + qJD(4) * t8 - qJD(5) * t38 + qJD(6) * t17;
t12 = [-qJD(2) * t9 + qJD(4) * t1 + qJD(5) * t13 + qJD(6) * t5, 0.2e1 * ((t139 * t239 + t140 * t240 + t291) * t349 + (-t156 * t159 + t291) * t350 + (t351 + m(5) * (-0.1e1 + t358) / 0.2e1) * t203 * t204) * qJD(2) + t3 * qJD(4) + t25 * qJD(5) + t11 * qJD(6) + t243, t2, t3 * qJD(2) + t18 * qJD(5) + t7 * qJD(6) + t247 + (t315 * t331 + t269 * t326 - Ifges(5,5) * t209 + Ifges(5,6) * t207 + Ifges(6,6) * t177 - t191 * t116 - t230 / 0.2e1 + t95 * t334 + t93 * t332 + t236 * t307 - mrSges(5,1) * t172 - (Ifges(6,5) + t356) * t228 + t368 * t99 + t367 * t353 + (m(7) * t245 + t278 - t281) * t190 + t245 * mrSges(7,3)) * qJD(4), qJD(2) * t25 + qJD(4) * t18 + qJD(6) * t23 - t242, t11 * qJD(2) + t7 * qJD(4) + t23 * qJD(5) + (-mrSges(7,1) * t48 - mrSges(7,2) * t47 + t230) * qJD(6) + t244; qJD(4) * t4 + qJD(5) * t26 + qJD(6) * t10 - t243, t31 * qJD(4), -t295 + t328, t33 * qJD(6) + t224 + (mrSges(5,2) * t279 - t203 * t303 - t367 * t156 + (t352 - t368) * t157) * qJD(4), t296, t294 + t33 * qJD(4) + (-mrSges(7,1) * t140 - mrSges(7,2) * t139) * qJD(6); t2, t295 + t328, t39 * qJD(4) (m(7) * t283 - t228 * t352 - t229 * t348 + t115 + t362) * qJD(4) - t72 * qJD(6) + t223, -t275, -t72 * qJD(4) + qJD(6) * t115 + t287; -qJD(2) * t4 + qJD(5) * t16 - qJD(6) * t6 - t247, -qJD(6) * t32 - t224, qJD(6) * t71 - t223, t62 * qJD(6), t288 (t179 * t190 + t249) * qJD(6) - t218; -qJD(2) * t26 - qJD(4) * t16 - qJD(6) * t22 + t242, -t296, t275, -t288, 0, -t180 * qJD(6) - t276; -qJD(2) * t10 + qJD(4) * t6 + qJD(5) * t22 - t244, t32 * qJD(4) - t294, -t71 * qJD(4) - t287, t218, t276, 0;];
Cq  = t12;
