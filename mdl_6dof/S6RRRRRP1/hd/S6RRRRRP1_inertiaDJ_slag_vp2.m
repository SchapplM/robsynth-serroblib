% Calculate time derivative of joint inertia matrix for
% S6RRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 00:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP1_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:55:38
% EndTime: 2019-03-10 00:55:53
% DurationCPUTime: 7.30s
% Computational Cost: add. (10507->483), mult. (22995->690), div. (0->0), fcn. (22220->8), ass. (0->212)
t336 = Ifges(6,4) + Ifges(7,4);
t335 = Ifges(6,1) + Ifges(7,1);
t334 = Ifges(6,2) + Ifges(7,2);
t196 = cos(qJ(5));
t333 = t336 * t196;
t192 = sin(qJ(5));
t332 = t336 * t192;
t286 = Ifges(6,5) + Ifges(7,5);
t323 = Ifges(6,6) + Ifges(7,6);
t330 = -t192 * t334 + t333;
t329 = t196 * t335 - t332;
t328 = t192 * t286 + t323 * t196;
t326 = t334 * t196 + t332;
t325 = t335 * t192 + t333;
t190 = t192 ^ 2;
t191 = t196 ^ 2;
t324 = t190 + t191;
t322 = Ifges(6,3) + Ifges(7,3);
t194 = sin(qJ(3));
t195 = sin(qJ(2));
t198 = cos(qJ(3));
t199 = cos(qJ(2));
t258 = t198 * t199;
t149 = -t194 * t195 + t258;
t150 = t194 * t199 + t198 * t195;
t193 = sin(qJ(4));
t197 = cos(qJ(4));
t110 = t149 * t193 + t150 * t197;
t248 = qJD(5) * t192;
t311 = qJD(2) + qJD(3);
t116 = t311 * t149;
t117 = t311 * t150;
t210 = t197 * t149 - t150 * t193;
t59 = qJD(4) * t210 + t116 * t197 - t117 * t193;
t271 = t196 * t59;
t207 = t110 * t248 - t271;
t247 = qJD(5) * t196;
t235 = t110 * t247;
t273 = t192 * t59;
t208 = t235 + t273;
t60 = qJD(4) * t110 + t116 * t193 + t197 * t117;
t321 = -t207 * t336 - t334 * t208 + t323 * t60;
t320 = -t335 * t207 - t208 * t336 + t286 * t60;
t319 = t330 * t110 - t323 * t210;
t282 = t329 * t110 - t286 * t210;
t318 = t330 * qJD(5);
t317 = t329 * qJD(5);
t254 = t286 * t247;
t209 = t317 * t192 + t318 * t196 + t325 * t247;
t203 = -t248 * t326 + t209;
t298 = -pkin(8) - pkin(7);
t172 = t298 * t195;
t151 = t194 * t172;
t173 = t298 * t199;
t121 = -t198 * t173 + t151;
t100 = pkin(9) * t149 + t121;
t120 = t198 * t172 + t173 * t194;
t99 = -pkin(9) * t150 + t120;
t69 = t100 * t193 - t197 * t99;
t70 = t100 * t197 + t193 * t99;
t67 = t196 * t70;
t181 = -pkin(2) * t199 - pkin(1);
t126 = -t149 * pkin(3) + t181;
t68 = -pkin(4) * t210 - t110 * pkin(10) + t126;
t34 = t192 * t68 + t67;
t314 = qJD(5) * t34;
t312 = t324 * t197;
t188 = t196 * qJD(6);
t284 = -qJ(6) - pkin(10);
t228 = qJD(5) * t284;
t131 = t192 * t228 + t188;
t281 = mrSges(7,3) * t196;
t125 = t131 * t281;
t154 = mrSges(7,1) * t248 + mrSges(7,2) * t247;
t291 = pkin(5) * t196;
t180 = -pkin(4) - t291;
t128 = t180 * t154;
t163 = -mrSges(7,1) * t196 + mrSges(7,2) * t192;
t186 = pkin(5) * t248;
t142 = t163 * t186;
t218 = mrSges(6,1) * t192 + mrSges(6,2) * t196;
t155 = t218 * qJD(5);
t292 = pkin(4) * t155;
t310 = t125 + t128 + t142 - t292;
t177 = pkin(3) * t193 + pkin(10);
t256 = -qJ(6) - t177;
t224 = qJD(5) * t256;
t250 = qJD(4) * t197;
t241 = pkin(3) * t250;
t114 = t192 * t224 + t196 * t241 + t188;
t107 = t114 * t281;
t293 = pkin(3) * t197;
t161 = t180 - t293;
t122 = t161 * t154;
t251 = qJD(4) * t193;
t242 = pkin(3) * t251;
t160 = t186 + t242;
t123 = t160 * t163;
t178 = -pkin(4) - t293;
t127 = t178 * t155;
t164 = -mrSges(6,1) * t196 + mrSges(6,2) * t192;
t143 = t164 * t242;
t223 = mrSges(6,3) * t241;
t170 = t190 * t223;
t171 = t191 * t223;
t309 = t107 + t122 + t123 + t127 + t143 + t170 + t171;
t308 = 2 * m(5);
t307 = 2 * m(6);
t306 = 2 * m(7);
t305 = -0.2e1 * mrSges(7,3);
t290 = t116 * pkin(9);
t252 = qJD(3) * t198;
t253 = qJD(3) * t194;
t81 = t150 * qJD(2) * t298 + t172 * t252 + t173 * t253;
t66 = -pkin(9) * t117 + t81;
t204 = (t258 * t298 - t151) * qJD(2);
t82 = -t121 * qJD(3) + t204;
t23 = qJD(4) * t70 + t193 * t66 - t197 * (t82 - t290);
t304 = 0.2e1 * t23;
t303 = 0.2e1 * t69;
t102 = qJD(2) * t195 * pkin(2) + pkin(3) * t117;
t302 = 0.2e1 * t102;
t301 = 0.2e1 * t181;
t300 = -0.2e1 * t192;
t297 = mrSges(7,3) * pkin(5);
t22 = t197 * t66 - t69 * qJD(4) + (-t172 * t253 + t173 * t252 + t204 - t290) * t193;
t31 = pkin(4) * t60 - pkin(10) * t59 + t102;
t243 = t192 * t31 + t196 * t22 + t68 * t247;
t5 = -t248 * t70 + t243;
t289 = t196 * t5;
t288 = t23 * t69;
t229 = -t192 * t22 + t196 * t31;
t6 = t229 - t314;
t287 = t6 * t192;
t276 = pkin(3) * qJD(4);
t179 = pkin(2) * t198 + pkin(3);
t262 = t193 * t194;
t104 = t179 * t250 + (-t194 * t251 + (t197 * t198 - t262) * qJD(3)) * pkin(2);
t275 = t104 * mrSges(5,2);
t261 = t194 * t197;
t105 = t179 * t251 + (t194 * t250 + (t193 * t198 + t261) * qJD(3)) * pkin(2);
t274 = t105 * t69;
t269 = t104 * t190;
t268 = t104 * t191;
t267 = t110 * t192;
t266 = t110 * t196;
t265 = t177 * t196;
t189 = t196 * qJ(6);
t134 = pkin(2) * t261 + t193 * t179;
t130 = pkin(10) + t134;
t257 = -qJ(6) - t130;
t249 = qJD(5) * t130;
t246 = -0.2e1 * t281;
t245 = mrSges(7,3) * t300;
t244 = 0.2e1 * t199;
t240 = m(7) * pkin(5) + mrSges(7,1);
t232 = t323 * t192;
t231 = -t248 / 0.2e1;
t33 = -t192 * t70 + t196 * t68;
t227 = t286 * t271 + t322 * t60;
t226 = t104 * t324;
t133 = -pkin(2) * t262 + t179 * t197;
t225 = qJD(5) * t257;
t129 = -pkin(4) - t133;
t220 = -(2 * Ifges(5,4)) - t232;
t219 = -t193 * mrSges(5,1) - t197 * mrSges(5,2);
t24 = -pkin(5) * t210 - t110 * t189 + t33;
t32 = -qJ(6) * t267 + t34;
t213 = -t192 * t32 - t196 * t24;
t212 = -t192 * t34 - t196 * t33;
t211 = -qJ(6) * t59 - qJD(6) * t110;
t206 = (-mrSges(4,1) * t194 - mrSges(4,2) * t198) * qJD(3) * pkin(2);
t17 = t208 * mrSges(7,1) - t207 * mrSges(7,2);
t103 = t105 * mrSges(5,1);
t124 = t129 - t291;
t106 = t124 * t154;
t113 = t129 * t155;
t74 = t104 * t196 + t192 * t225 + t188;
t71 = t74 * t281;
t96 = t186 + t105;
t83 = t96 * t163;
t94 = t105 * t164;
t97 = mrSges(6,3) * t269;
t98 = mrSges(6,3) * t268;
t202 = -t103 + t106 + t113 - t275 + t71 + t83 + t94 + t97 + t98 + t203;
t3 = -qJ(6) * t235 + (-qJD(5) * t70 + t211) * t192 + t243;
t36 = pkin(5) * t267 + t69;
t8 = pkin(5) * t208 + t23;
t201 = -t22 * mrSges(5,2) + mrSges(6,3) * t289 + Ifges(5,5) * t59 + t36 * t154 + t69 * t155 + t8 * t163 + t3 * t281 - (-t323 * t248 + t254) * t210 / 0.2e1 + t320 * t192 / 0.2e1 + t321 * t196 / 0.2e1 - t318 * t267 / 0.2e1 + t317 * t266 / 0.2e1 + t319 * t231 + t282 * t247 / 0.2e1 + (t164 - mrSges(5,1)) * t23 + t326 * (-t235 / 0.2e1 - t273 / 0.2e1) + t325 * (t110 * t231 + t271 / 0.2e1) + (-Ifges(5,6) + t328 / 0.2e1) * t60;
t200 = t82 * mrSges(4,1) - t81 * mrSges(4,2) + Ifges(4,5) * t116 - Ifges(4,6) * t117 + t201;
t165 = pkin(10) * t196 + t189;
t162 = t284 * t192;
t145 = t189 + t265;
t144 = t256 * t192;
t132 = -qJD(6) * t192 + t196 * t228;
t119 = t130 * t196 + t189;
t118 = t257 * t192;
t115 = (-qJD(6) - t241) * t192 + t196 * t224;
t79 = -mrSges(6,1) * t210 - mrSges(6,3) * t266;
t78 = -mrSges(7,1) * t210 - mrSges(7,3) * t266;
t77 = mrSges(6,2) * t210 - mrSges(6,3) * t267;
t76 = mrSges(7,2) * t210 - mrSges(7,3) * t267;
t75 = (-qJD(6) - t104) * t192 + t196 * t225;
t73 = t218 * t110;
t72 = (mrSges(7,1) * t192 + mrSges(7,2) * t196) * t110;
t28 = -mrSges(6,2) * t60 - mrSges(6,3) * t208;
t27 = -mrSges(7,2) * t60 - mrSges(7,3) * t208;
t26 = mrSges(6,1) * t60 + mrSges(6,3) * t207;
t25 = mrSges(7,1) * t60 + mrSges(7,3) * t207;
t18 = mrSges(6,1) * t208 - mrSges(6,2) * t207;
t1 = pkin(5) * t60 + t211 * t196 + (-t67 + (qJ(6) * t110 - t68) * t192) * qJD(5) + t229;
t2 = [(mrSges(4,1) * t117 + mrSges(4,2) * t116) * t301 + 0.2e1 * (t116 * t149 - t117 * t150) * Ifges(4,4) + 0.2e1 * (-t116 * t120 - t117 * t121 + t149 * t81 - t150 * t82) * mrSges(4,3) - (mrSges(5,1) * t302 - 0.2e1 * mrSges(5,3) * t22 + ((2 * Ifges(5,2)) + t322) * t60 + t220 * t59 + t227) * t210 + t18 * t303 + t73 * t304 + (t1 * t24 + t3 * t32 + t36 * t8) * t306 + (t33 * t6 + t34 * t5 + t288) * t307 + (t102 * t126 + t22 * t70 + t288) * t308 + 0.2e1 * t3 * t76 + 0.2e1 * t5 * t77 + 0.2e1 * t1 * t78 + 0.2e1 * t6 * t79 + 0.2e1 * t8 * t72 + 0.2e1 * m(4) * (t120 * t82 + t121 * t81) + 0.2e1 * t32 * t27 + 0.2e1 * t33 * t26 + 0.2e1 * t34 * t28 + 0.2e1 * t36 * t17 + 0.2e1 * t24 * t25 + (0.2e1 * t126 * mrSges(5,2) + mrSges(5,3) * t303 - t192 * t319 + t282 * t196) * t59 + (mrSges(5,2) * t302 + mrSges(5,3) * t304 + 0.2e1 * Ifges(5,1) * t59 + t320 * t196 - t321 * t192 + (t196 * t286 + t220) * t60 + (-t192 * t282 - t196 * t319 + t210 * t328) * qJD(5)) * t110 - 0.2e1 * t149 * Ifges(4,2) * t117 + 0.2e1 * t116 * t150 * Ifges(4,1) + ((-mrSges(3,2) * pkin(1) + Ifges(3,4) * t199) * t244 + (m(4) * pkin(2) * t301 + 0.2e1 * pkin(2) * (-mrSges(4,1) * t149 + mrSges(4,2) * t150) - 0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t195 + (Ifges(3,1) - Ifges(3,2)) * t244) * t195) * qJD(2) + 0.2e1 * (mrSges(5,1) * t126 - mrSges(5,3) * t70) * t60; (t104 * t210 + t105 * t110 - t133 * t59 - t134 * t60) * mrSges(5,3) + t200 + (Ifges(3,5) * t199 - Ifges(3,6) * t195 + (-mrSges(3,1) * t199 + mrSges(3,2) * t195) * pkin(7)) * qJD(2) + t124 * t17 + t129 * t18 + t118 * t25 + t119 * t27 + t105 * t73 + t96 * t72 + t74 * t76 + t75 * t78 + m(6) * (t129 * t23 + t274) + (m(4) * (t194 * t81 + t198 * t82 + (-t120 * t194 + t121 * t198) * qJD(3)) + (-t198 * t116 - t194 * t117 + (t149 * t198 + t150 * t194) * qJD(3)) * mrSges(4,3)) * pkin(2) + (m(6) * (-t104 * t33 - t130 * t6 - t249 * t34) - t130 * t26 - t104 * t79 - t77 * t249 + (-qJD(5) * t32 - t1) * mrSges(7,3) + (-t6 - t314) * mrSges(6,3)) * t192 + (m(6) * (t104 * t34 + t130 * t5 - t249 * t33) + t130 * t28 + t104 * t77 - t79 * t249 + (-mrSges(6,3) * t33 - mrSges(7,3) * t24) * qJD(5)) * t196 + m(7) * (t1 * t118 + t119 * t3 + t124 * t8 + t24 * t75 + t32 * t74 + t36 * t96) + m(5) * (t104 * t70 - t133 * t23 + t134 * t22 + t274); (t104 * t134 - t105 * t133) * t308 + (t118 * t75 + t119 * t74 + t124 * t96) * t306 + (t105 * t129 + t130 * t226) * t307 + 0.2e1 * t113 + 0.2e1 * t106 + (t118 * t246 + (t119 * t305 - t326) * t192) * qJD(5) - 0.2e1 * t103 + 0.2e1 * t98 + 0.2e1 * t97 + 0.2e1 * t94 + t75 * t245 - 0.2e1 * t275 + 0.2e1 * t206 + t209 + 0.2e1 * t83 + 0.2e1 * t71; t200 + t28 * t265 + t178 * t18 + t160 * t72 + t161 * t17 + t144 * t25 + t145 * t27 + t114 * t76 + t115 * t78 + m(7) * (t1 * t144 + t114 * t32 + t115 * t24 + t145 * t3 + t160 * t36 + t161 * t8) + (-t6 * mrSges(6,3) - t1 * mrSges(7,3) - t177 * t26) * t192 + m(6) * (-t177 * t287 + t178 * t23 + t5 * t265) + (m(5) * (t193 * t22 - t197 * t23) + (-t193 * t60 - t197 * t59) * mrSges(5,3) + ((t210 * mrSges(5,3) - t192 * t79 + t196 * t77 + m(6) * (-t192 * t33 + t196 * t34) + m(5) * t70) * t197 + (t110 * mrSges(5,3) + t73 + (m(6) + m(5)) * t69) * t193) * qJD(4)) * pkin(3) + (t213 * mrSges(7,3) + t212 * mrSges(6,3) + (m(6) * t212 - t192 * t77 - t196 * t79) * t177) * qJD(5); t202 + ((-t115 - t75) * t192 + ((-t118 - t144) * t196 + (-t119 - t145) * t192) * qJD(5)) * mrSges(7,3) + m(6) * (t105 * t178 + (t268 + t269) * t177) + m(7) * (t114 * t119 + t115 * t118 + t124 * t160 + t144 * t75 + t145 * t74 + t161 * t96) + (m(5) * (t104 * t193 - t105 * t197) + (m(6) * (t129 * t193 + t130 * t312) + m(5) * (-t133 * t193 + t134 * t197) + t219) * qJD(4)) * pkin(3) + t206 + t309; 0.2e1 * t170 + 0.2e1 * t107 + t115 * t245 + 0.2e1 * t143 + 0.2e1 * t127 + 0.2e1 * t123 + 0.2e1 * t122 + (t114 * t145 + t115 * t144 + t160 * t161) * t306 + 0.2e1 * t171 + (t144 * t246 + (t145 * t305 - t326) * t192) * qJD(5) + 0.2e1 * (m(6) * (t177 * t312 + t178 * t193) + t219) * t276 + t209; t201 + t72 * t186 + t180 * t17 + t162 * t25 + t165 * t27 + t131 * t76 + t132 * t78 + m(7) * (t1 * t162 + t131 * t32 + t132 * t24 + t165 * t3 + t180 * t8 + t186 * t36) + (m(6) * (-t247 * t33 - t248 * t34 - t287 + t289) + t196 * t28 - t192 * t26 - t79 * t247 - t77 * t248) * pkin(10) + (qJD(5) * t213 - t1 * t192) * mrSges(7,3) + (qJD(5) * t212 - t287) * mrSges(6,3) + (-m(6) * t23 - t18) * pkin(4); t202 + m(6) * (-pkin(4) * t105 + pkin(10) * t226) + m(7) * (t118 * t132 + t119 * t131 + t124 * t186 + t162 * t75 + t165 * t74 + t180 * t96) + ((-t132 - t75) * t192 + ((-t118 - t162) * t196 + (-t119 - t165) * t192) * qJD(5)) * mrSges(7,3) + t310; (m(6) * (-pkin(4) * t193 + pkin(10) * t312) + t219) * t276 + t203 + m(7) * (t114 * t165 + t115 * t162 + t131 * t145 + t132 * t144 + t160 * t180 + t161 * t186) + ((-t115 - t132) * t192 + ((-t144 - t162) * t196 + (-t145 - t165) * t192) * qJD(5)) * mrSges(7,3) + t309 + t310; -0.2e1 * t292 + 0.2e1 * t142 + 0.2e1 * t128 + (t131 * t165 + t132 * t162 + t180 * t186) * t306 + 0.2e1 * t125 + (t132 * t300 + 0.2e1 * (-t162 * t196 - t165 * t192) * qJD(5)) * mrSges(7,3) + t203; mrSges(6,1) * t6 + mrSges(7,1) * t1 - mrSges(6,2) * t5 - mrSges(7,2) * t3 - t59 * t232 + (m(7) * t1 + t25) * pkin(5) - t328 * t110 * qJD(5) + t227; -mrSges(7,2) * t74 + t240 * t75 - t218 * t104 + ((-mrSges(6,1) * t130 - t297) * t196 + (mrSges(6,2) * t130 - t323) * t192) * qJD(5) + t254; -mrSges(7,2) * t114 + t240 * t115 - t218 * t241 + ((-mrSges(6,1) * t177 - t297) * t196 + (mrSges(6,2) * t177 - t323) * t192) * qJD(5) + t254; -mrSges(7,2) * t131 + t240 * t132 + ((-mrSges(6,1) * pkin(10) - t297) * t196 + (mrSges(6,2) * pkin(10) - t323) * t192) * qJD(5) + t254; 0; m(7) * t8 + t17; m(7) * t96 + t154; m(7) * t160 + t154; m(7) * t186 + t154; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t2(1) t2(2) t2(4) t2(7) t2(11) t2(16); t2(2) t2(3) t2(5) t2(8) t2(12) t2(17); t2(4) t2(5) t2(6) t2(9) t2(13) t2(18); t2(7) t2(8) t2(9) t2(10) t2(14) t2(19); t2(11) t2(12) t2(13) t2(14) t2(15) t2(20); t2(16) t2(17) t2(18) t2(19) t2(20) t2(21);];
Mq  = res;
