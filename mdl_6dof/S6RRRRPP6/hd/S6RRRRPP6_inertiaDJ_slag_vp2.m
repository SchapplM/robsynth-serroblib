% Calculate time derivative of joint inertia matrix for
% S6RRRRPP6
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
% Datum: 2019-03-09 21:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP6_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP6_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP6_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP6_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:11:21
% EndTime: 2019-03-09 21:11:36
% DurationCPUTime: 6.92s
% Computational Cost: add. (4552->487), mult. (11088->673), div. (0->0), fcn. (9022->6), ass. (0->216)
t325 = Ifges(5,1) + Ifges(7,3);
t324 = -Ifges(5,4) + Ifges(7,6);
t323 = Ifges(7,4) + Ifges(6,5);
t322 = Ifges(5,5) + Ifges(7,5);
t321 = Ifges(7,2) + Ifges(6,3);
t320 = Ifges(6,6) - Ifges(7,6);
t319 = -Ifges(6,1) - Ifges(7,1) - Ifges(5,3);
t297 = -Ifges(6,4) + t322;
t318 = Ifges(5,6) - t323;
t184 = cos(qJ(2));
t271 = cos(qJ(3));
t219 = qJD(2) * t271;
t207 = t184 * t219;
t182 = sin(qJ(3));
t183 = sin(qJ(2));
t240 = qJD(3) * t183;
t220 = t182 * t240;
t189 = t207 - t220;
t243 = qJD(2) * t183;
t218 = qJD(3) * t271;
t242 = qJD(2) * t184;
t222 = t182 * t242;
t302 = -t183 * t218 - t222;
t317 = -Ifges(4,5) * t189 - t302 * Ifges(4,6) - Ifges(4,3) * t243;
t316 = -mrSges(5,1) + mrSges(6,2);
t299 = qJD(3) + qJD(4);
t181 = sin(qJ(4));
t270 = cos(qJ(4));
t315 = (-mrSges(5,2) * t270 + t181 * t316) * pkin(3) * qJD(4);
t314 = 2 * mrSges(6,1) + 2 * mrSges(5,3);
t313 = mrSges(5,2) - mrSges(6,3);
t312 = mrSges(7,2) + mrSges(6,3);
t311 = Ifges(3,1) - Ifges(3,2);
t223 = t270 * t182;
t137 = t181 * t271 + t223;
t99 = t299 * t137;
t60 = t181 * t222 + t183 * t99 - t207 * t270;
t208 = t270 * t271;
t186 = t299 * t208;
t246 = t182 * t183;
t226 = t181 * t246;
t61 = -qJD(4) * t226 + t137 * t242 - t181 * t220 + t183 * t186;
t310 = t323 * t243 + t320 * t60 + t321 * t61;
t309 = t322 * t243 + t324 * t61 - t325 * t60;
t247 = t181 * t182;
t289 = t299 * t247;
t98 = t289 - t186;
t308 = t320 * t98 + t321 * t99;
t307 = t324 * t99 - t325 * t98;
t125 = t137 * t183;
t126 = t183 * t208 - t226;
t306 = t321 * t125 - t320 * t126 - t323 * t184;
t305 = t324 * t125 + t325 * t126 - t322 * t184;
t136 = -t208 + t247;
t304 = t321 * t136 - t320 * t137;
t303 = t324 * t136 + t325 * t137;
t146 = -pkin(2) * t184 - t183 * pkin(8) - pkin(1);
t224 = t184 * t271;
t164 = pkin(7) * t224;
t122 = t182 * t146 + t164;
t144 = (pkin(2) * t183 - pkin(8) * t184) * qJD(2);
t241 = qJD(3) * t182;
t68 = t182 * t144 + t146 * t218 + (-t183 * t219 - t184 * t241) * pkin(7);
t231 = pkin(7) * t243;
t245 = t271 * t144 + t182 * t231;
t69 = -t122 * qJD(3) + t245;
t301 = -t182 * t69 + t271 * t68;
t217 = qJD(4) * t270;
t209 = pkin(3) * t217;
t159 = t209 + qJD(5);
t269 = pkin(3) * t181;
t165 = qJ(5) + t269;
t300 = -t136 * t159 - t165 * t99;
t298 = -0.2e1 * pkin(1) * qJD(2);
t34 = mrSges(6,1) * t61 - mrSges(6,3) * t243;
t109 = -pkin(9) * t246 + t122;
t239 = qJD(4) * t181;
t32 = (pkin(3) * t183 - pkin(9) * t224) * qJD(2) + (-t164 + (pkin(9) * t183 - t146) * t182) * qJD(3) + t245;
t52 = pkin(9) * t302 + t68;
t135 = t271 * t146;
t225 = t183 * t271;
t268 = pkin(7) * t182;
t89 = -pkin(9) * t225 + t135 + (-pkin(3) - t268) * t184;
t9 = -t109 * t239 + t181 * t32 + t89 * t217 + t270 * t52;
t5 = -qJ(5) * t243 + qJD(5) * t184 - t9;
t295 = -m(6) * t5 - t34;
t36 = -t60 * mrSges(6,1) + mrSges(6,2) * t243;
t195 = t109 * t217 + t181 * t52 + t89 * t239 - t270 * t32;
t7 = -pkin(4) * t243 + t195;
t294 = m(6) * t7 + t36;
t3 = -t61 * pkin(5) - t5;
t35 = -t61 * mrSges(7,1) + mrSges(7,2) * t243;
t293 = m(7) * t3 + t35;
t292 = t319 * t243 + t297 * t60 + t318 * t61;
t291 = -t297 * t98 - t318 * t99;
t113 = t125 * mrSges(6,1) + mrSges(6,3) * t184;
t40 = t270 * t109 + t181 * t89;
t29 = qJ(5) * t184 - t40;
t290 = m(6) * t29 + t113;
t115 = t126 * mrSges(6,1) - mrSges(6,2) * t184;
t117 = -mrSges(5,1) * t184 - t126 * mrSges(5,3);
t39 = -t181 * t109 + t270 * t89;
t30 = t184 * pkin(4) - t39;
t288 = m(6) * t30 + t115 - t117;
t287 = 2 * m(5);
t286 = 0.2e1 * m(6);
t285 = 0.2e1 * m(7);
t284 = -0.2e1 * mrSges(7,1);
t283 = m(5) * pkin(3);
t278 = -pkin(9) - pkin(8);
t178 = t183 * pkin(7);
t267 = t98 * mrSges(6,1);
t266 = t98 * mrSges(7,1);
t265 = t98 * mrSges(5,3);
t264 = t99 * mrSges(7,1);
t263 = t99 * mrSges(5,3);
t261 = pkin(4) + qJ(6);
t255 = Ifges(4,4) * t182;
t253 = t136 * mrSges(5,3);
t252 = t137 * mrSges(7,1);
t248 = t159 * t165;
t145 = pkin(3) * t246 + t178;
t244 = qJ(5) * qJD(5);
t234 = t271 * pkin(8);
t233 = t270 * pkin(3);
t232 = pkin(3) * t239;
t230 = pkin(8) * t241;
t177 = pkin(7) * t242;
t176 = pkin(3) * t241;
t228 = Ifges(4,4) * t271;
t120 = -pkin(3) * t302 + t177;
t221 = t183 * t242;
t214 = -t241 / 0.2e1;
t138 = t278 * t223;
t152 = pkin(9) * t271 + t234;
t110 = t152 * t181 - t138;
t211 = t137 * t232;
t210 = pkin(8) * t218;
t168 = -pkin(3) * t271 - pkin(2);
t167 = -t233 - pkin(4);
t33 = -t60 * mrSges(7,1) - mrSges(7,3) * t243;
t205 = Ifges(4,5) * t218 - Ifges(4,6) * t241;
t204 = -qJ(5) * t126 + t145;
t111 = t152 * t270 + t247 * t278;
t190 = qJD(3) * t152;
t64 = -t138 * t299 + t152 * t239 + t181 * t190;
t65 = t152 * t217 + t270 * t190 + t278 * t289;
t203 = t110 * t65 - t111 * t64;
t202 = -qJ(5) * t99 - qJD(5) * t136;
t201 = qJ(5) * t159 + qJD(5) * t165;
t200 = mrSges(4,1) * t182 + mrSges(4,2) * t271;
t199 = -Ifges(4,1) * t182 - t228;
t198 = -Ifges(4,2) * t182 + t228;
t197 = -Ifges(4,2) * t271 - t255;
t196 = Ifges(4,1) * t271 - t255;
t194 = qJ(5) * t98 - qJD(5) * t137 + t176;
t192 = -t137 * qJ(5) + t168;
t191 = qJ(5) * t60 - qJD(5) * t126 + t120;
t25 = -pkin(5) * t99 - t64;
t26 = -t98 * pkin(5) + t65;
t187 = t25 * mrSges(7,2) - t26 * mrSges(7,3) + t316 * t65 + t291;
t1 = -t60 * pkin(5) + t184 * qJD(6) - t243 * t261 + t195;
t185 = -mrSges(5,1) * t195 - t9 * mrSges(5,2) + t7 * mrSges(6,2) + t3 * mrSges(7,2) - t5 * mrSges(6,3) - t1 * mrSges(7,3) - t292;
t162 = -qJ(6) + t167;
t158 = -qJD(6) + t232;
t143 = -t184 * mrSges(4,1) - mrSges(4,3) * t225;
t142 = mrSges(4,2) * t184 - mrSges(4,3) * t246;
t141 = t196 * qJD(3);
t140 = t198 * qJD(3);
t139 = t200 * qJD(3);
t124 = -Ifges(4,5) * t184 + t183 * t196;
t123 = -Ifges(4,6) * t184 + t183 * t198;
t121 = -t184 * t268 + t135;
t119 = -mrSges(4,2) * t243 + mrSges(4,3) * t302;
t118 = mrSges(4,1) * t243 - mrSges(4,3) * t189;
t116 = mrSges(5,2) * t184 - t125 * mrSges(5,3);
t114 = -t125 * mrSges(7,1) - mrSges(7,2) * t184;
t112 = t126 * mrSges(7,1) + mrSges(7,3) * t184;
t107 = Ifges(5,4) * t137 - Ifges(5,2) * t136;
t106 = -Ifges(6,2) * t137 + Ifges(6,6) * t136;
t102 = -mrSges(6,2) * t136 - mrSges(6,3) * t137;
t101 = mrSges(5,1) * t136 + mrSges(5,2) * t137;
t100 = -mrSges(7,2) * t137 + mrSges(7,3) * t136;
t87 = -mrSges(4,1) * t302 + mrSges(4,2) * t189;
t86 = t136 * pkin(4) + t192;
t82 = -t136 * pkin(5) + t111;
t81 = pkin(5) * t137 + t110;
t80 = -mrSges(6,2) * t125 - mrSges(6,3) * t126;
t79 = mrSges(5,1) * t125 + mrSges(5,2) * t126;
t78 = -mrSges(7,2) * t126 + mrSges(7,3) * t125;
t77 = t199 * t240 + (Ifges(4,5) * t183 + t184 * t196) * qJD(2);
t76 = t197 * t240 + (Ifges(4,6) * t183 + t184 * t198) * qJD(2);
t74 = Ifges(5,4) * t126 - Ifges(5,2) * t125 - Ifges(5,6) * t184;
t73 = -Ifges(6,4) * t184 - Ifges(6,2) * t126 + Ifges(6,6) * t125;
t67 = t136 * t261 + t192;
t66 = pkin(4) * t125 + t204;
t50 = -Ifges(5,4) * t98 - Ifges(5,2) * t99;
t49 = Ifges(6,2) * t98 + Ifges(6,6) * t99;
t45 = mrSges(7,2) * t98 + mrSges(7,3) * t99;
t44 = mrSges(5,1) * t99 - mrSges(5,2) * t98;
t43 = -mrSges(6,2) * t99 + mrSges(6,3) * t98;
t38 = -mrSges(5,2) * t243 - mrSges(5,3) * t61;
t37 = mrSges(5,1) * t243 + mrSges(5,3) * t60;
t31 = t125 * t261 + t204;
t27 = pkin(4) * t99 + t194;
t24 = -t125 * pkin(5) - t29;
t22 = t126 * pkin(5) + t184 * qJ(6) + t30;
t21 = mrSges(7,2) * t60 + mrSges(7,3) * t61;
t20 = mrSges(5,1) * t61 - mrSges(5,2) * t60;
t19 = -mrSges(6,2) * t61 + mrSges(6,3) * t60;
t17 = -Ifges(5,4) * t60 - Ifges(5,2) * t61 + Ifges(5,6) * t243;
t16 = Ifges(6,4) * t243 + Ifges(6,2) * t60 + Ifges(6,6) * t61;
t12 = qJD(6) * t136 + t261 * t99 + t194;
t11 = pkin(4) * t61 + t191;
t4 = qJD(6) * t125 + t261 * t61 + t191;
t2 = [(t297 * t126 - t318 * t125 + (Ifges(4,5) * t271 - Ifges(4,6) * t182 - (2 * Ifges(3,4))) * t183) * t243 + t77 * t225 + 0.2e1 * t68 * t142 + 0.2e1 * t69 * t143 + 0.2e1 * t145 * t20 + (0.2e1 * Ifges(3,4) * t242 + (-Ifges(4,3) + t311 + t319) * t243 + mrSges(3,2) * t298 + t292 + t317) * t184 + 0.2e1 * t1 * t112 + 0.2e1 * t5 * t113 + 0.2e1 * t3 * t114 + 0.2e1 * t7 * t115 + 0.2e1 * t9 * t116 + 0.2e1 * t120 * t79 + 0.2e1 * t121 * t118 + 0.2e1 * t122 * t119 + 0.2e1 * t11 * t80 + 0.2e1 * t66 * t19 + 0.2e1 * t4 * t78 + 0.2e1 * t29 * t34 + 0.2e1 * t24 * t35 + 0.2e1 * t30 * t36 + 0.2e1 * t39 * t37 + 0.2e1 * t40 * t38 + 0.2e1 * t31 * t21 + 0.2e1 * t22 * t33 + (t120 * t145 - t195 * t39 + t40 * t9) * t287 - 0.2e1 * t195 * t117 + (-t16 + t309) * t126 + (-t17 + t310) * t125 + (mrSges(3,1) * t298 + t311 * t242) * t183 + (t73 - t305) * t60 + (-t74 + t306) * t61 + t302 * t123 + t189 * t124 + 0.2e1 * t87 * t178 - t76 * t246 + 0.2e1 * pkin(7) * t200 * t221 + (t1 * t22 + t24 * t3 + t31 * t4) * t285 + (t11 * t66 + t29 * t5 + t30 * t7) * t286 + 0.2e1 * m(4) * (pkin(7) ^ 2 * t221 + t121 * t69 + t122 * t68); (Ifges(4,5) * t182 + Ifges(4,6) * t271 - t136 * t318 + t297 * t137) * t243 / 0.2e1 + (-t49 / 0.2e1 + t307 / 0.2e1) * t126 + (-t50 / 0.2e1 + t308 / 0.2e1) * t125 + m(5) * (t120 * t168 + t145 * t176) + m(6) * (t11 * t86 + t27 * t66) + (t136 * t5 + t137 * t7) * mrSges(6,1) + (-t74 / 0.2e1 + t29 * mrSges(6,1) + t306 / 0.2e1) * t99 + (-t17 / 0.2e1 - t3 * mrSges(7,1)) * t136 + t39 * t265 + t139 * t178 + t1 * t252 + (-t107 / 0.2e1 + t304 / 0.2e1) * t61 + Ifges(3,5) * t242 + t119 * t234 + t79 * t176 + mrSges(3,2) * t231 + t182 * t77 / 0.2e1 + t168 * t20 + t145 * t44 + (t106 / 0.2e1 - t303 / 0.2e1) * t60 + t26 * t112 + t25 * t114 + t120 * t101 + (t73 / 0.2e1 - t305 / 0.2e1) * t98 + t4 * t100 + t11 * t102 + t81 * t33 + t82 * t35 + t86 * t19 - pkin(2) * t87 + t66 * t43 + t67 * t21 + t12 * t78 + t27 * t80 + t31 * t45 + t123 * t214 + (-t16 / 0.2e1 + t195 * mrSges(5,3)) * t137 + (m(5) * t195 + t294 - t37) * t110 - t302 * t197 / 0.2e1 - (t183 * t214 + t207 / 0.2e1) * t199 + t309 * t137 / 0.2e1 + t310 * t136 / 0.2e1 + (m(4) * ((-t121 * t271 - t122 * t182) * qJD(3) + t301) - t182 * t118) * pkin(8) + (-t121 * t218 - t122 * t241 + t301) * mrSges(4,3) + (m(5) * t9 + t295 + t38) * t111 + (-m(5) * t40 - t116 + t290) * t64 - (t205 + t291) * t184 / 0.2e1 + (-m(5) * t39 + t288) * t65 + m(7) * (t1 * t81 + t12 * t31 + t22 * t26 + t24 * t25 + t3 * t82 + t4 * t67) + t271 * t76 / 0.2e1 - t30 * t267 - t40 * t263 - t24 * t264 - t22 * t266 - t9 * t253 - t140 * t246 / 0.2e1 - Ifges(3,6) * t243 + (-m(4) * pkin(2) - mrSges(4,1) * t271 + t182 * mrSges(4,2) - mrSges(3,1)) * t177 - t142 * t230 - t143 * t210 + t124 * t218 / 0.2e1 + t141 * t225 / 0.2e1; t271 * t140 - 0.2e1 * pkin(2) * t139 + 0.2e1 * t12 * t100 + 0.2e1 * t27 * t102 + t182 * t141 + 0.2e1 * t168 * t44 + 0.2e1 * t86 * t43 + 0.2e1 * t67 * t45 + (-t271 * t199 + (0.2e1 * pkin(3) * t101 + t197) * t182) * qJD(3) + (t168 * t176 + t203) * t287 + (t27 * t86 + t203) * t286 + (t12 * t67 + t25 * t82 + t26 * t81) * t285 + (-t111 * t314 + t284 * t82 - t107 + t304) * t99 + (-t110 * t314 + t284 * t81 + t106 - t303) * t98 + (0.2e1 * t26 * mrSges(7,1) + t314 * t65 + t307 - t49) * t137 + (t25 * t284 + t314 * t64 + t308 - t50) * t136; t185 + (-t270 * t195 + t181 * t9 + (-t181 * t39 + t270 * t40) * qJD(4)) * t283 + t116 * t209 + t158 * t112 + t162 * t33 - t68 * mrSges(4,2) + t69 * mrSges(4,1) + m(7) * (t1 * t162 + t158 * t22) + t38 * t269 + t37 * t233 + t294 * t167 + t288 * t232 + (t293 + t295) * t165 + (m(7) * t24 + t114 - t290) * t159 - t317; t187 + (-t270 * t65 + (t110 * t181 + t111 * t270) * qJD(4)) * t283 + mrSges(4,2) * t230 - t263 * t269 + t233 * t265 + t205 - mrSges(4,1) * t210 + m(6) * (t110 * t232 + t111 * t159 + t167 * t65) + m(7) * (t158 * t81 + t159 * t82 + t162 * t26 + t165 * t25) - t167 * t267 + t158 * t252 - t162 * t266 - t209 * t253 + mrSges(5,3) * t211 + (-m(6) * t165 - t181 * t283 + t313) * t64 + t300 * mrSges(7,1) + (t211 + t300) * mrSges(6,1); -0.2e1 * t158 * mrSges(7,3) + (t158 * t162 + t248) * t285 + (t167 * t232 + t248) * t286 + 0.2e1 * t312 * t159 + 0.2e1 * t315; t185 - t261 * t33 - qJD(6) * t112 - pkin(4) * t36 + (-t113 + t114) * qJD(5) + (-t34 + t35) * qJ(5) + m(7) * (qJ(5) * t3 + qJD(5) * t24 - qJD(6) * t22 - t1 * t261) + m(6) * (-pkin(4) * t7 - qJ(5) * t5 - qJD(5) * t29); t187 + t313 * t64 + m(6) * (-pkin(4) * t65 - qJ(5) * t64 + qJD(5) * t111) + m(7) * (qJ(5) * t25 + qJD(5) * t82 - qJD(6) * t81 - t26 * t261) + (pkin(4) * t98 + t202) * mrSges(6,1) + (-qJD(6) * t137 + t261 * t98 + t202) * mrSges(7,1); (qJD(6) - t158) * mrSges(7,3) + t315 + m(7) * (-qJD(6) * t162 - t158 * t261 + t201) + m(6) * (-pkin(4) * t232 + t201) + t312 * (qJD(5) + t159); 0.2e1 * m(6) * t244 + 0.2e1 * qJD(6) * mrSges(7,3) + 0.2e1 * m(7) * (qJD(6) * t261 + t244) + 0.2e1 * t312 * qJD(5); m(7) * t1 + t294 + t33; (-mrSges(7,1) - mrSges(6,1)) * t98 + m(6) * t65 + m(7) * t26; m(6) * t232 + m(7) * t158; -m(7) * qJD(6); 0; t293; m(7) * t25 - t264; m(7) * t159; m(7) * qJD(5); 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t2(1) t2(2) t2(4) t2(7) t2(11) t2(16); t2(2) t2(3) t2(5) t2(8) t2(12) t2(17); t2(4) t2(5) t2(6) t2(9) t2(13) t2(18); t2(7) t2(8) t2(9) t2(10) t2(14) t2(19); t2(11) t2(12) t2(13) t2(14) t2(15) t2(20); t2(16) t2(17) t2(18) t2(19) t2(20) t2(21);];
Mq  = res;
