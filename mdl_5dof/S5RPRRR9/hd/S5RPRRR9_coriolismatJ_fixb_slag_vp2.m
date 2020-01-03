% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRR9
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
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRR9_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR9_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR9_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR9_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR9_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:07:23
% EndTime: 2019-12-31 19:07:31
% DurationCPUTime: 4.07s
% Computational Cost: add. (13271->344), mult. (26265->480), div. (0->0), fcn. (30724->8), ass. (0->216)
t202 = sin(qJ(5));
t204 = cos(qJ(5));
t185 = Ifges(6,5) * t202 + Ifges(6,6) * t204;
t200 = sin(pkin(9));
t201 = cos(pkin(9));
t330 = sin(qJ(3));
t332 = cos(qJ(3));
t180 = -t330 * t200 + t332 * t201;
t181 = t332 * t200 + t330 * t201;
t203 = sin(qJ(4));
t331 = cos(qJ(4));
t221 = t203 * t180 + t331 * t181;
t363 = t221 * t185;
t365 = Ifges(5,6) * t221;
t232 = mrSges(6,1) * t204 - t202 * mrSges(6,2);
t320 = pkin(6) + qJ(2);
t242 = t320 * t201;
t243 = t320 * t200;
t166 = -t330 * t242 - t332 * t243;
t137 = -t181 * pkin(7) + t166;
t167 = t332 * t242 - t330 * t243;
t216 = t180 * pkin(7) + t167;
t355 = t137 * t203 + t331 * t216;
t375 = t355 * t232;
t378 = t355 * mrSges(5,1);
t163 = -t331 * t180 + t181 * t203;
t195 = Ifges(6,4) * t204;
t360 = -Ifges(6,2) * t202 + t195;
t371 = Ifges(6,6) * t221 - t163 * t360;
t387 = t204 * t371;
t316 = Ifges(6,4) * t202;
t189 = t204 * Ifges(6,1) - t316;
t372 = Ifges(6,5) * t221 - t189 * t163;
t388 = t202 * t372;
t76 = -t331 * t137 + t203 * t216;
t392 = t76 * mrSges(5,2);
t394 = t363 / 0.2e1 - t365 - t375 - t378 + t392 + t387 / 0.2e1 + t388 / 0.2e1;
t393 = -t375 / 0.2e1 - t378 / 0.2e1 + t387 / 0.4e1 + t388 / 0.4e1 + t392 / 0.2e1;
t391 = t202 * t76;
t390 = t203 * t76;
t389 = t204 * t76;
t324 = t355 * t76;
t301 = t204 * mrSges(6,2);
t305 = t202 * mrSges(6,1);
t184 = t301 + t305;
t377 = t184 * t163;
t386 = -m(6) * t355 + t377;
t109 = t184 * t221;
t383 = t355 * t109 - t76 * t377;
t382 = -t163 / 0.2e1;
t380 = mrSges(6,3) * t163;
t161 = Ifges(5,4) * t163;
t379 = Ifges(6,5) * t163;
t328 = pkin(4) * t221;
t118 = pkin(8) * t163 + t328;
t326 = t181 * pkin(3);
t80 = t118 + t326;
t38 = t202 * t80 - t389;
t295 = t38 * t204;
t37 = t204 * t80 + t391;
t296 = t37 * t202;
t229 = t295 - t296;
t373 = t363 / 0.4e1 - t365 / 0.2e1;
t346 = -t221 / 0.2e1;
t370 = t221 / 0.2e1;
t369 = mrSges(6,1) * t221;
t368 = mrSges(6,2) * t221;
t367 = mrSges(5,3) * t221;
t317 = Ifges(5,4) * t221;
t313 = Ifges(6,2) * t204;
t186 = t313 + t316;
t273 = t202 * t186;
t176 = -t273 / 0.2e1;
t318 = Ifges(6,1) * t202;
t188 = t195 + t318;
t268 = t204 * t188;
t362 = t176 + t268 / 0.2e1;
t194 = Ifges(6,5) * t204;
t310 = Ifges(6,6) * t202;
t361 = t194 - t310;
t198 = t202 ^ 2;
t199 = t204 ^ 2;
t359 = t198 + t199;
t358 = t268 / 0.4e1 - t273 / 0.4e1;
t65 = t189 * t221 + t379;
t297 = t204 * t65;
t62 = Ifges(6,6) * t163 + t221 * t360;
t304 = t202 * t62;
t357 = t297 / 0.4e1 - t304 / 0.4e1;
t356 = (-t194 / 0.2e1 + t310 / 0.2e1) * t163;
t274 = t163 * t202;
t111 = mrSges(6,3) * t274 - t368;
t271 = t204 * t111;
t269 = t204 * t163;
t114 = mrSges(6,3) * t269 + t369;
t276 = t202 * t114;
t351 = -pkin(4) / 0.2e1;
t352 = m(6) / 0.2e1;
t354 = -pkin(4) * t355 * t352 - t377 * t351 + (t229 * t352 - t276 / 0.2e1 + t271 / 0.2e1) * pkin(8) + (t295 / 0.2e1 - t296 / 0.2e1) * mrSges(6,3) + t393;
t353 = t181 ^ 2;
t350 = m(5) * pkin(3);
t349 = m(6) * pkin(3);
t348 = -mrSges(6,1) / 0.2e1;
t347 = mrSges(6,2) / 0.2e1;
t343 = -t163 / 0.4e1;
t342 = t163 / 0.2e1;
t341 = t163 / 0.4e1;
t329 = pkin(3) * t203;
t192 = pkin(8) + t329;
t338 = -t192 / 0.2e1;
t263 = t331 * pkin(3);
t193 = -t263 - pkin(4);
t337 = t193 / 0.2e1;
t336 = -t202 / 0.2e1;
t335 = t202 / 0.2e1;
t334 = -t204 / 0.2e1;
t333 = t204 / 0.2e1;
t327 = pkin(4) * t184;
t315 = Ifges(5,5) * t163;
t309 = Ifges(6,3) * t221;
t281 = t221 * t202;
t262 = mrSges(6,3) * t281;
t112 = -mrSges(6,2) * t163 - t262;
t280 = t221 * t204;
t115 = t163 * mrSges(6,1) - mrSges(6,3) * t280;
t116 = -t163 * Ifges(5,2) + t317;
t117 = Ifges(5,1) * t221 - t161;
t253 = -pkin(2) * t201 - pkin(1);
t168 = -pkin(3) * t180 + t253;
t160 = t163 * mrSges(5,2);
t307 = t221 * mrSges(5,1);
t240 = -t160 + t307;
t241 = t181 * mrSges(4,1) + t180 * mrSges(4,2);
t246 = -t269 / 0.2e1;
t248 = t274 / 0.2e1;
t79 = pkin(4) * t163 - pkin(8) * t221 + t168;
t33 = -t202 * t355 + t204 * t79;
t34 = t202 * t79 + t204 * t355;
t59 = Ifges(6,3) * t163 + t221 * t361;
t1 = t383 - t371 * t281 / 0.2e1 + t372 * t280 / 0.2e1 + (mrSges(5,1) * t326 - Ifges(5,1) * t370 - t361 * t342) * t163 + (t59 - t317) * t370 + mrSges(5,2) * t221 * t326 + t309 * t342 + t116 * t346 - t353 * Ifges(4,4) + (Ifges(4,4) * t180 + (Ifges(4,1) - Ifges(4,2)) * t181) * t180 + m(6) * (t33 * t37 + t34 * t38 + t324) + t253 * t241 + t65 * t246 + t62 * t248 + t34 * t111 + t38 * t112 + t33 * t114 + t37 * t115 + (-Ifges(5,2) * t221 + t117 - t161) * t382 + (m(5) * t326 + t240) * t168;
t308 = t1 * qJD(1);
t300 = t204 * t34;
t110 = t202 * t380 - t368;
t113 = t204 * t380 + t369;
t43 = t118 * t204 + t391;
t44 = t202 * t118 - t389;
t4 = m(6) * (t33 * t43 + t34 * t44 + t324) + t44 * t112 + t34 * t110 + t43 * t115 + t33 * t113 + (t372 * t333 + t371 * t336 + t59 / 0.2e1 - t116 / 0.2e1 + t168 * mrSges(5,1) - t317 / 0.2e1) * t221 + (-t297 / 0.2e1 + t304 / 0.2e1 - t117 / 0.2e1 + t161 / 0.2e1 - t168 * mrSges(5,2) + t356 + (Ifges(6,3) / 0.2e1 - Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t221) * t163 + t383;
t294 = t4 * qJD(1);
t293 = t43 * t202;
t292 = t44 * t204;
t106 = t232 * t221;
t5 = t34 * t115 - t76 * t106 + (t65 * t335 + t62 * t333 + mrSges(6,3) * t300 + t185 * t342 + (-t188 * t334 + t176) * t221) * t221 + (-t112 - t262) * t33;
t291 = t5 * qJD(1);
t289 = t76 * t221;
t270 = t204 * t112;
t275 = t202 * t115;
t10 = (t109 + t367) * t221 + (t180 ^ 2 + t353) * mrSges(4,3) - (-t163 * mrSges(5,3) + t270 - t275) * t163 + m(6) * (t289 - (-t202 * t33 + t300) * t163) + m(5) * (-t163 * t355 + t289) + m(4) * (-t166 * t181 + t167 * t180) + (m(3) * qJ(2) + mrSges(3,3)) * (t200 ^ 2 + t201 ^ 2);
t287 = qJD(1) * t10;
t238 = t359 * t163;
t249 = t232 * t370;
t266 = t350 / 0.2e1;
t208 = (-t192 * t238 + t193 * t221) * t352 - t249 + (-t163 * t203 - t221 * t331) * t266 - t359 * t380 / 0.2e1;
t214 = (t202 * t38 + t204 * t37) * t352 + t111 * t335 + t114 * t333 + t181 * t266;
t11 = t208 - t214 - t240 - t241;
t286 = t11 * qJD(1);
t244 = t199 / 0.2e1 + t198 / 0.2e1;
t211 = -t244 * t380 + t160 / 0.2e1 + (-pkin(8) * t238 - t328) * t352 - t249;
t215 = -m(6) * (t202 * t44 + t204 * t43) / 0.2e1 + mrSges(5,2) * t342 + t110 * t336 + t113 * t334;
t13 = 0.2e1 * t346 * mrSges(5,1) + t211 + t215;
t285 = t13 * qJD(1);
t219 = (-t301 / 0.2e1 - t305 / 0.2e1) * t163;
t222 = t275 / 0.2e1 - t270 / 0.2e1;
t17 = -t219 + t222;
t279 = t17 * qJD(1);
t278 = t193 * t184;
t277 = t202 * t113;
t272 = t204 * t110;
t265 = mrSges(6,3) * t293;
t264 = mrSges(6,3) * t292;
t254 = t76 * t184 / 0.2e1;
t252 = t202 * t331;
t251 = t204 * t331;
t236 = t189 * t335 + t333 * t360 + t362;
t235 = -mrSges(6,2) * t331 / 0.2e1;
t234 = -t252 / 0.2e1;
t233 = mrSges(6,3) * t244;
t228 = t292 - t293;
t205 = (t228 * t192 + t193 * t355) * t352 + Ifges(5,5) * t382 - t377 * t337 + t273 * t341 + t277 * t338 + t109 * t329 / 0.2e1 + t268 * t343 + t192 * t272 / 0.2e1 - t265 / 0.2e1 + t264 / 0.2e1 + ((t34 * t251 - t33 * t252 + t390) * t352 + t115 * t234 + t112 * t251 / 0.2e1) * pkin(3) + t373 + t393;
t3 = -t205 - t315 / 0.2e1 - t358 * t163 + t354 + t373;
t210 = (-mrSges(5,1) - t232) * t329 + (t359 * mrSges(6,3) - mrSges(5,2)) * t263;
t220 = t359 * t331;
t93 = (t220 * t192 + t193 * t203) * t349 + t210;
t227 = -t3 * qJD(1) + t93 * qJD(3);
t119 = t236 + t278;
t213 = (t189 / 0.4e1 - t186 / 0.4e1 - t313 / 0.4e1) * t204 + (-t188 / 0.4e1 - t360 / 0.4e1 - t195 / 0.2e1 - t318 / 0.4e1) * t202;
t207 = (-t192 * t233 + t213) * t221 + t106 * t337 + t254;
t217 = -t309 / 0.2e1 + t37 * t348 + t38 * t347;
t7 = t194 * t341 + (t112 * t338 - t62 / 0.4e1 + (t343 + t382) * Ifges(6,6)) * t202 + (t115 * t338 + t65 / 0.4e1 + t379 / 0.2e1) * t204 + t207 + t217;
t226 = t7 * qJD(1) + t119 * qJD(3);
t225 = t44 * t347 + t43 * t348;
t224 = t361 * t341;
t223 = t112 * t336 + t115 * t334;
t120 = t236 - t327;
t69 = (-t193 / 0.2e1 + pkin(4) / 0.2e1) * t184 + (pkin(3) * t235 - t188 / 0.2e1 - t360 / 0.2e1) * t204 + (t263 * t348 - t189 / 0.2e1 + t186 / 0.2e1) * t202;
t209 = -pkin(8) * t233 + t213;
t212 = t223 * pkin(8) + t106 * t351 + t254 + t357;
t9 = (-0.3e1 / 0.4e1 * t310 + 0.3e1 / 0.4e1 * t194) * t163 + (-Ifges(6,3) / 0.2e1 + t209) * t221 + t212 + t225;
t218 = t9 * qJD(1) - t69 * qJD(3) + t120 * qJD(4);
t70 = t278 / 0.2e1 - t327 / 0.2e1 + (mrSges(6,1) * t234 + t204 * t235) * pkin(3) + t236;
t18 = -t219 - t222;
t15 = t208 + t214;
t14 = mrSges(5,1) * t370 - t307 / 0.2e1 + t211 - t215;
t8 = Ifges(6,3) * t370 + t209 * t221 + t212 + t224 - t225 + t356;
t6 = Ifges(6,5) * t246 + Ifges(6,6) * t248 + t223 * t192 + t207 - t217 + t224 + t357;
t2 = t205 - (Ifges(5,5) / 0.2e1 + t358) * t163 + (-Ifges(5,6) / 0.2e1 + t185 / 0.4e1) * t221 + t354;
t12 = [qJD(2) * t10 + qJD(3) * t1 + qJD(4) * t4 - qJD(5) * t5, qJD(3) * t15 + qJD(4) * t14 + qJD(5) * t18 + t287, t15 * qJD(2) + t2 * qJD(4) + t6 * qJD(5) + t308 + (-t329 * t367 - t315 + Ifges(4,5) * t180 - Ifges(4,6) * t181 + (-t331 * t355 - t390) * t350 - t167 * mrSges(4,1) - t166 * mrSges(4,2) - (-mrSges(5,3) * t263 + t362) * t163 - t386 * t193 + (m(6) * t229 + t271 - t276) * t192 + t229 * mrSges(6,3) + t394) * qJD(3), t14 * qJD(2) + t2 * qJD(3) + t8 * qJD(5) + t294 + (t264 - t265 + (-t268 / 0.2e1 + t273 / 0.2e1 - Ifges(5,5)) * t163 + t386 * pkin(4) + (m(6) * t228 + t272 - t277) * pkin(8) + t394) * qJD(4), -t291 + t18 * qJD(2) + t6 * qJD(3) + t8 * qJD(4) + (-t34 * mrSges(6,1) - t33 * mrSges(6,2) - t363) * qJD(5); -qJD(3) * t11 - qJD(4) * t13 - qJD(5) * t17 - t287, 0, -t286, -t285, -qJD(5) * t184 - t279; qJD(2) * t11 - qJD(4) * t3 + qJD(5) * t7 - t308, t286, qJD(4) * t93 + qJD(5) * t119, ((-pkin(4) * t203 + pkin(8) * t220) * t349 + t210) * qJD(4) + t70 * qJD(5) + t227, t70 * qJD(4) + (-t192 * t232 + t361) * qJD(5) + t226; qJD(2) * t13 + qJD(3) * t3 + qJD(5) * t9 - t294, t285, -qJD(5) * t69 - t227, t120 * qJD(5), (-pkin(8) * t232 + t361) * qJD(5) + t218; qJD(2) * t17 - qJD(3) * t7 - qJD(4) * t9 + t291, t279, qJD(4) * t69 - t226, -t218, 0;];
Cq = t12;
