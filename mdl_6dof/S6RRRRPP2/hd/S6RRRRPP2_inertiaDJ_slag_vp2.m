% Calculate time derivative of joint inertia matrix for
% S6RRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:48:59
% EndTime: 2019-03-09 20:49:13
% DurationCPUTime: 6.08s
% Computational Cost: add. (4652->442), mult. (10618->588), div. (0->0), fcn. (9031->6), ass. (0->188)
t349 = Ifges(7,4) + Ifges(6,5);
t347 = Ifges(7,2) + Ifges(6,3);
t202 = sin(qJ(3));
t203 = sin(qJ(2));
t205 = cos(qJ(3));
t206 = cos(qJ(2));
t152 = t202 * t203 - t205 * t206;
t335 = Ifges(6,4) + Ifges(5,5);
t346 = -Ifges(7,5) + t335;
t348 = t152 * t346;
t339 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t204 = cos(qJ(4));
t345 = t349 * t204;
t201 = sin(qJ(4));
t344 = t349 * t201;
t342 = Ifges(6,6) - Ifges(7,6);
t341 = t201 ^ 2 + t204 ^ 2;
t340 = t201 * t347 + t345;
t292 = Ifges(5,4) * t201;
t337 = t204 * t339 - t292 + t344;
t336 = -t204 * t347 + t344;
t334 = Ifges(6,2) + Ifges(5,3);
t153 = t202 * t206 + t203 * t205;
t318 = qJD(2) + qJD(3);
t105 = t318 * t153;
t266 = qJD(4) * t201;
t104 = t318 * t152;
t279 = t104 * t204;
t214 = t153 * t266 + t279;
t265 = qJD(4) * t204;
t247 = t153 * t265;
t280 = t104 * t201;
t215 = t247 - t280;
t333 = t105 * t342 - t214 * t349 + t215 * t347;
t332 = t152 * t342 + t153 * t340;
t169 = -t204 * mrSges(5,1) + t201 * mrSges(5,2);
t331 = -mrSges(4,1) + t169;
t330 = t340 * qJD(4);
t329 = Ifges(6,6) * t266 + t265 * t335;
t285 = pkin(2) * qJD(3);
t252 = t205 * t285;
t328 = t341 * t252;
t305 = -pkin(8) - pkin(7);
t177 = t305 * t203;
t178 = t305 * t206;
t125 = t177 * t202 - t178 * t205;
t254 = pkin(2) * qJD(2) * t203;
t45 = pkin(3) * t105 + pkin(9) * t104 + t254;
t249 = qJD(2) * t305;
t164 = t203 * t249;
t234 = t206 * t249;
t319 = t177 * t205 + t178 * t202;
t53 = qJD(3) * t319 + t164 * t205 + t202 * t234;
t187 = -pkin(2) * t206 - pkin(1);
t88 = pkin(3) * t152 - pkin(9) * t153 + t187;
t262 = t201 * t45 + t204 * t53 + t265 * t88;
t10 = -t125 * t266 + t262;
t256 = t125 * t265 + t201 * t53 + t266 * t88;
t11 = t204 * t45 - t256;
t326 = t10 * t204 - t11 * t201;
t321 = qJ(5) * t105 + qJD(5) * t152;
t4 = t10 + t321;
t112 = t201 * t125;
t47 = t204 * t88 - t112;
t41 = -pkin(4) * t152 - t47;
t325 = t204 * t4 + t265 * t41;
t324 = (-Ifges(5,4) + t349) * t215 - t339 * t214 + t346 * t105;
t260 = t153 * t337 + t348;
t323 = t337 * qJD(4);
t291 = Ifges(5,4) * t204;
t322 = t201 * t339 + t291 - t345;
t29 = mrSges(5,1) * t215 - mrSges(5,2) * t214;
t54 = qJD(3) * t125 + t164 * t202 - t205 * t234;
t320 = m(5) * t54 + t29;
t196 = t201 * qJ(5);
t222 = pkin(4) * t204 + t196;
t317 = (-mrSges(4,2) + (mrSges(6,2) + mrSges(5,3)) * t341) * t252;
t316 = 0.2e1 * m(5);
t315 = 2 * m(6);
t314 = 2 * m(7);
t313 = -2 * Ifges(4,4);
t312 = 0.2e1 * t54;
t311 = -0.2e1 * t319;
t229 = t201 * mrSges(6,1) - t204 * mrSges(6,3);
t155 = t229 * qJD(4);
t310 = 0.2e1 * t155;
t156 = -mrSges(7,1) * t266 + mrSges(7,2) * t265;
t309 = 0.2e1 * t156;
t168 = mrSges(7,1) * t204 + mrSges(7,2) * t201;
t308 = 0.2e1 * t168;
t307 = -0.2e1 * t201;
t306 = -0.2e1 * t204;
t207 = -pkin(4) - pkin(5);
t297 = -Ifges(5,6) - Ifges(7,6);
t296 = pkin(9) - qJ(6);
t277 = t153 * t201;
t90 = -mrSges(5,2) * t152 - mrSges(5,3) * t277;
t94 = -mrSges(6,2) * t277 + mrSges(6,3) * t152;
t295 = t90 + t94;
t276 = t153 * t204;
t92 = mrSges(5,1) * t152 - mrSges(5,3) * t276;
t93 = -mrSges(6,1) * t152 + mrSges(6,2) * t276;
t294 = -t92 + t93;
t293 = mrSges(6,2) * t201;
t286 = Ifges(7,5) * t204;
t282 = t319 * t54;
t48 = t125 * t204 + t201 * t88;
t281 = qJ(5) * t204;
t278 = t319 * t202;
t274 = t201 * t205;
t272 = t204 * t205;
t185 = pkin(2) * t202 + pkin(9);
t271 = -qJ(6) + t185;
t270 = t328 * t185;
t269 = t328 * pkin(9);
t264 = qJD(5) * t204;
t225 = -Ifges(5,2) * t201 + t291;
t63 = Ifges(5,6) * t152 + t153 * t225;
t261 = -t63 + t332;
t258 = 0.2e1 * qJD(4);
t257 = pkin(3) + t222;
t255 = m(6) * t264;
t253 = t202 * t285;
t40 = qJ(5) * t152 + t48;
t251 = mrSges(7,3) * t266;
t132 = pkin(4) * t266 - qJ(5) * t265 - qJD(5) * t201;
t186 = -pkin(2) * t205 - pkin(3);
t182 = qJ(6) * t266;
t170 = t296 * t204;
t240 = t265 / 0.2e1;
t239 = 0.2e1 * t254;
t144 = t271 * t204;
t231 = -qJD(6) + t252;
t230 = t201 * mrSges(5,1) + t204 * mrSges(5,2);
t167 = -t204 * mrSges(6,1) - t201 * mrSges(6,3);
t221 = pkin(4) * t201 - t281;
t140 = t186 - t222;
t220 = t331 * t253;
t219 = Ifges(6,6) * t215 + t105 * t334 - t279 * t335;
t218 = qJ(6) * t104 - qJD(6) * t153;
t34 = -t105 * mrSges(7,1) + mrSges(7,3) * t214;
t217 = t201 * t207 + t281;
t216 = t201 * t297 - t286;
t126 = -pkin(5) * t266 - t132;
t36 = -t105 * mrSges(6,1) - mrSges(6,2) * t214;
t213 = qJ(5) * t251 + (mrSges(6,2) - mrSges(7,3)) * t264 + t329;
t28 = -mrSges(7,1) * t215 - mrSges(7,2) * t214;
t212 = (-mrSges(6,2) * qJ(5) + t297) * t201 + (-mrSges(6,2) * pkin(4) - mrSges(7,3) * t207 - Ifges(7,5)) * t204;
t211 = -m(6) * t222 + t167 + t169;
t160 = t225 * qJD(4);
t173 = Ifges(5,2) * t204 + t292;
t210 = (t160 - t330) * t204 + (-t173 + t336) * t266 + t322 * t265 + t323 * t201;
t35 = mrSges(5,1) * t105 + mrSges(5,3) * t214;
t38 = -mrSges(5,2) * t105 - mrSges(5,3) * t215;
t39 = -mrSges(6,2) * t215 + mrSges(6,3) * t105;
t7 = -pkin(4) * t105 - t11;
t209 = (t38 + t39) * t204 + (-t35 + t36) * t201 + (-t201 * t295 + t204 * t294) * qJD(4) + m(5) * (-t265 * t47 - t266 * t48 + t326) + m(6) * (t201 * t7 - t266 * t40 + t325);
t1 = t153 * t182 + t207 * t105 + (t218 - t45) * t204 + t256;
t13 = -t221 * t104 + (qJD(4) * t222 - t264) * t153 + t54;
t157 = t230 * qJD(4);
t2 = qJ(6) * t247 + (-qJD(4) * t125 - t218) * t201 + t262 + t321;
t22 = -Ifges(5,4) * t214 - Ifges(5,2) * t215 + Ifges(5,6) * t105;
t26 = t112 + (-qJ(6) * t153 - t88) * t204 + t207 * t152;
t31 = qJ(6) * t277 + t40;
t42 = t153 * t217 + t319;
t52 = t153 * t221 - t319;
t8 = -t217 * t104 + (t264 + (t204 * t207 - t196) * qJD(4)) * t153 - t54;
t208 = t7 * t293 + t331 * t54 + t325 * mrSges(6,2) + ((-t201 * t48 - t204 * t47) * qJD(4) + t326) * mrSges(5,3) - t215 * t173 / 0.2e1 + (-Ifges(5,6) * t266 + t329) * t152 / 0.2e1 + t260 * t240 + t324 * t201 / 0.2e1 - t322 * t279 / 0.2e1 + t323 * t276 / 0.2e1 + t332 * t266 / 0.2e1 - t333 * t204 / 0.2e1 + ((Ifges(5,6) - Ifges(6,6)) * t204 + t335 * t201) * t105 / 0.2e1 + (-t26 * t204 * mrSges(7,3) - t40 * t293 - t152 * (Ifges(7,6) * t201 + t286) / 0.2e1) * qJD(4) - t105 * (Ifges(7,5) * t201 - Ifges(7,6) * t204) / 0.2e1 + t204 * t22 / 0.2e1 + t13 * t167 + t8 * t168 + t52 * t155 + t42 * t156 - Ifges(4,5) * t104 - Ifges(4,6) * t105 - t53 * mrSges(4,2) + t336 * (t153 * t240 - t280 / 0.2e1) - (t322 * t153 + t63) * t266 / 0.2e1 - t319 * t157 + t31 * t251 + (-t1 * t201 - t2 * t204) * mrSges(7,3) + (t330 / 0.2e1 - t160 / 0.2e1) * t277;
t197 = t204 * pkin(5);
t190 = mrSges(6,2) * t265;
t166 = t296 * t201;
t143 = t271 * t201;
t141 = t197 + t257;
t133 = qJD(4) * t170 - qJD(6) * t201;
t131 = -pkin(9) * t266 - qJD(6) * t204 + t182;
t128 = -t140 + t197;
t127 = t132 + t253;
t114 = t126 - t253;
t103 = qJD(4) * t144 + t201 * t231;
t102 = -t185 * t266 + t204 * t231 + t182;
t91 = -mrSges(7,1) * t152 - mrSges(7,3) * t276;
t89 = mrSges(7,2) * t152 + mrSges(7,3) * t277;
t82 = t230 * t153;
t81 = (-mrSges(7,1) * t201 + mrSges(7,2) * t204) * t153;
t80 = t229 * t153;
t37 = mrSges(7,2) * t105 + mrSges(7,3) * t215;
t27 = mrSges(6,1) * t215 + mrSges(6,3) * t214;
t3 = [(mrSges(4,2) * t239 + mrSges(4,3) * t312 - 0.2e1 * Ifges(4,1) * t104 + t324 * t204 + (-t22 + t333) * t201 + (t313 + t346 * t204 + (Ifges(6,6) + t297) * t201) * t105 + ((t152 * t297 + t261) * t204 + (-t260 - t348) * t201) * qJD(4)) * t153 + (mrSges(4,1) * t239 - 0.2e1 * t53 * mrSges(4,3) + ((2 * Ifges(4,2)) + (2 * Ifges(7,3)) + t334) * t105 - (t313 + t216) * t104 + t219) * t152 + 0.2e1 * t2 * t89 + 0.2e1 * t10 * t90 + 0.2e1 * t1 * t91 + 0.2e1 * t11 * t92 + 0.2e1 * t7 * t93 + 0.2e1 * t4 * t94 + 0.2e1 * t13 * t80 + 0.2e1 * t8 * t81 + 0.2e1 * t47 * t35 + 0.2e1 * t48 * t38 + 0.2e1 * t52 * t27 + 0.2e1 * t26 * t34 + 0.2e1 * t31 * t37 + 0.2e1 * t40 * t39 + 0.2e1 * t41 * t36 + 0.2e1 * t42 * t28 - (0.2e1 * t187 * mrSges(4,2) + mrSges(4,3) * t311 + t261 * t201 + t260 * t204) * t104 + 0.2e1 * m(4) * (t125 * t53 + t187 * t254 - t282) + (t10 * t48 + t11 * t47 - t282) * t316 + 0.2e1 * ((-t203 ^ 2 + t206 ^ 2) * Ifges(3,4) - pkin(1) * (mrSges(3,1) * t203 + mrSges(3,2) * t206) + (-Ifges(3,2) + Ifges(3,1)) * t203 * t206) * qJD(2) + t29 * t311 + t82 * t312 + (t1 * t26 + t2 * t31 + t42 * t8) * t314 + (t13 * t52 + t4 * t40 + t41 * t7) * t315 + 0.2e1 * (mrSges(4,1) * t187 - mrSges(4,3) * t125) * t105; (Ifges(3,5) * t206 - Ifges(3,6) * t203 + (-mrSges(3,1) * t206 + mrSges(3,2) * t203) * pkin(7)) * qJD(2) + m(6) * (t127 * t52 + t13 * t140) + t140 * t27 + t143 * t34 + t144 * t37 + t127 * t80 + t128 * t28 + t114 * t81 + t102 * t89 + t103 * t91 + t208 + m(7) * (t1 * t143 + t102 * t31 + t103 * t26 + t114 * t42 + t128 * t8 + t144 * t2) + t209 * t185 + (m(4) * (t202 * t53 - t205 * t54) + (t104 * t205 - t105 * t202) * mrSges(4,3) + ((mrSges(4,3) * t153 + t82) * t202 + (-t152 * mrSges(4,3) + t201 * t294 + t204 * t295) * t205 + m(6) * (t272 * t40 + t274 * t41) + m(5) * (t272 * t48 - t274 * t47 - t278) + m(4) * (t125 * t205 - t278)) * qJD(3)) * pkin(2) + t320 * t186; 0.2e1 * t220 + 0.2e1 * t127 * t167 + t114 * t308 + 0.2e1 * t186 * t157 + t140 * t310 + t128 * t309 + 0.2e1 * t317 + (t186 * t253 + t270) * t316 + (t102 * t144 + t103 * t143 + t114 * t128) * t314 + (t127 * t140 + t270) * t315 + t210 + (t102 * t306 + t103 * t307 + (-t143 * t204 + t144 * t201) * t258) * mrSges(7,3); m(7) * (t1 * t166 + t126 * t42 + t131 * t31 + t133 * t26 + t141 * t8 + t170 * t2) + t170 * t37 - t257 * t27 + t166 * t34 + t132 * t80 + t133 * t91 + t141 * t28 + t126 * t81 + t131 * t89 + t208 + m(6) * (-t13 * t257 + t132 * t52) + t209 * pkin(9) - t320 * pkin(3); ((-t102 - t131) * t204 + (-t103 - t133) * t201 + ((-t143 - t166) * t204 + (t144 + t170) * t201) * qJD(4)) * mrSges(7,3) + m(5) * (-pkin(3) * t253 + t269) + m(7) * (t102 * t170 + t103 * t166 + t114 * t141 + t126 * t128 + t131 * t144 + t133 * t143) + m(6) * (-t127 * t257 + t132 * t140 + t269) + (t186 - pkin(3)) * t157 + (t141 + t128) * t156 + (t140 - t257) * t155 + t317 + t220 + (t126 + t114) * t168 + (t132 + t127) * t167 + t210; t126 * t308 + t141 * t309 - 0.2e1 * pkin(3) * t157 - t257 * t310 + 0.2e1 * (-m(6) * t257 + t167) * t132 + t210 + (t126 * t141 + t131 * t170 + t133 * t166) * t314 + (t131 * t306 + t133 * t307 + (-t166 * t204 + t170 * t201) * t258) * mrSges(7,3); t219 + t207 * t34 + Ifges(7,3) * t105 - pkin(4) * t36 - t10 * mrSges(5,2) + t11 * mrSges(5,1) - t7 * mrSges(6,1) + t2 * mrSges(7,2) + t4 * mrSges(6,3) - t1 * mrSges(7,1) + (-t201 * t346 + t204 * t297) * t153 * qJD(4) + (t89 + t94) * qJD(5) + (t37 + t39) * qJ(5) + m(7) * (qJ(5) * t2 + qJD(5) * t31 + t1 * t207) + m(6) * (-pkin(4) * t7 + qJ(5) * t4 + qJD(5) * t40) - t216 * t104; -t103 * mrSges(7,1) + t102 * mrSges(7,2) + t185 * t255 + m(7) * (qJ(5) * t102 + qJD(5) * t144 + t103 * t207) + (-m(6) * t221 - t229 - t230) * t252 + (t185 * t211 + t212) * qJD(4) + t213; -t133 * mrSges(7,1) + t131 * mrSges(7,2) + pkin(9) * t255 + m(7) * (qJ(5) * t131 + qJD(5) * t170 + t133 * t207) + (pkin(9) * t211 + t212) * qJD(4) + t213; 0.2e1 * (mrSges(7,2) + mrSges(6,3) + (m(6) + m(7)) * qJ(5)) * qJD(5); m(6) * t7 + m(7) * t1 + t34 + t36; m(7) * t103 - mrSges(7,3) * t265 + t190 + (t185 * t265 + t201 * t252) * m(6); m(7) * t133 + t190 + (m(6) * pkin(9) - mrSges(7,3)) * t265; 0; 0; m(7) * t8 + t28; m(7) * t114 + t156; m(7) * t126 + t156; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t3(1) t3(2) t3(4) t3(7) t3(11) t3(16); t3(2) t3(3) t3(5) t3(8) t3(12) t3(17); t3(4) t3(5) t3(6) t3(9) t3(13) t3(18); t3(7) t3(8) t3(9) t3(10) t3(14) t3(19); t3(11) t3(12) t3(13) t3(14) t3(15) t3(20); t3(16) t3(17) t3(18) t3(19) t3(20) t3(21);];
Mq  = res;
