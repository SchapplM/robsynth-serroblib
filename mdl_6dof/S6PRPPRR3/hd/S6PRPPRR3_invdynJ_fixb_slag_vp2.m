% Calculate vector of inverse dynamics joint torques for
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPPRR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPPRR3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:21:29
% EndTime: 2019-03-08 19:21:43
% DurationCPUTime: 8.99s
% Computational Cost: add. (3903->493), mult. (7854->681), div. (0->0), fcn. (5596->12), ass. (0->237)
t155 = sin(pkin(11));
t158 = cos(pkin(11));
t167 = cos(qJ(2));
t157 = sin(pkin(6));
t260 = qJD(1) * t157;
t227 = t167 * t260;
t164 = sin(qJ(2));
t228 = t164 * t260;
t78 = t155 * t227 - t158 * t228;
t343 = qJD(3) * t155 - t78;
t342 = m(6) + m(7);
t163 = sin(qJ(5));
t166 = cos(qJ(5));
t206 = -pkin(5) * t163 + pkin(9) * t166;
t168 = -pkin(2) - pkin(3);
t128 = t158 * qJ(3) + t155 * t168;
t120 = -pkin(8) + t128;
t274 = t120 * t166;
t341 = qJD(5) * t206 - qJD(6) * t274 + t343;
t292 = mrSges(3,1) + mrSges(4,1);
t328 = mrSges(3,2) - mrSges(4,3);
t340 = t164 * t328 - t167 * t292;
t339 = m(4) + m(3);
t162 = sin(qJ(6));
t165 = cos(qJ(6));
t253 = qJD(5) * t165;
t258 = qJD(2) * t163;
t117 = t162 * t258 + t253;
t338 = t117 / 0.2e1;
t245 = t162 * qJD(5);
t118 = t165 * t258 - t245;
t301 = -t118 / 0.2e1;
t257 = qJD(2) * t166;
t147 = qJD(6) + t257;
t337 = t147 / 0.2e1;
t160 = cos(pkin(6));
t333 = qJDD(1) * t160;
t142 = qJDD(4) - t333;
t336 = m(5) * t142;
t294 = -qJD(2) / 0.2e1;
t235 = mrSges(6,3) * t258;
t278 = qJD(5) * mrSges(6,1) + mrSges(7,1) * t117 + mrSges(7,2) * t118 + t235;
t144 = -qJD(1) * t160 + qJD(4);
t191 = qJD(3) - t227;
t108 = qJD(2) * t168 + t191;
t124 = qJD(2) * qJ(3) + t228;
t63 = t155 * t108 + t158 * t124;
t54 = -qJD(2) * pkin(8) + t63;
t37 = t144 * t163 + t166 * t54;
t29 = qJD(5) * pkin(9) + t37;
t207 = pkin(5) * t166 + pkin(9) * t163;
t62 = t108 * t158 - t155 * t124;
t53 = qJD(2) * pkin(4) - t62;
t40 = qJD(2) * t207 + t53;
t10 = t162 * t40 + t165 * t29;
t242 = qJD(2) * qJD(5);
t125 = -qJDD(2) * t166 + t163 * t242;
t126 = -qJDD(2) * t163 - t166 * t242;
t219 = qJD(2) * t260;
t138 = t164 * t219;
t241 = qJDD(1) * t157;
t97 = t167 * t241 - t138;
t184 = qJDD(3) - t97;
t73 = qJDD(2) * t168 + t184;
t243 = qJD(2) * qJD(3);
t315 = qJDD(2) * qJ(3) + t243;
t139 = t167 * t219;
t98 = t164 * t241 + t139;
t74 = t98 + t315;
t22 = -t155 * t74 + t158 * t73;
t20 = qJDD(2) * pkin(4) - t22;
t15 = -pkin(5) * t125 - pkin(9) * t126 + t20;
t23 = t155 * t73 + t158 * t74;
t21 = -qJDD(2) * pkin(8) + t23;
t252 = qJD(5) * t166;
t254 = qJD(5) * t163;
t7 = t163 * t142 + t144 * t252 + t166 * t21 - t254 * t54;
t5 = qJDD(5) * pkin(9) + t7;
t9 = -t162 * t29 + t165 * t40;
t1 = qJD(6) * t9 + t15 * t162 + t165 * t5;
t2 = -qJD(6) * t10 + t15 * t165 - t162 * t5;
t205 = t1 * t165 - t162 * t2;
t249 = qJD(6) * t165;
t251 = qJD(6) * t162;
t334 = -t10 * t251 - t9 * t249 + t205;
t36 = t144 * t166 - t163 * t54;
t28 = -qJD(5) * pkin(5) - t36;
t331 = m(7) * t28;
t330 = -t125 / 0.2e1;
t329 = -t126 / 0.2e1;
t255 = qJD(3) * t158;
t319 = -t120 * t254 + t166 * t255;
t127 = -t155 * qJ(3) + t158 * t168;
t119 = pkin(4) - t127;
t82 = t119 + t207;
t172 = qJD(6) * t82 + t319;
t265 = t162 * t166;
t89 = (t155 * t164 + t158 * t167) * t157;
t81 = qJD(1) * t89;
t327 = -t162 * t172 + t341 * t165 + t265 * t81;
t263 = t165 * t166;
t326 = t341 * t162 + t165 * t172 - t263 * t81;
t107 = t155 * t263 - t158 * t162;
t222 = t155 * t254;
t325 = -qJD(6) * t107 + t162 * t222 - (t155 * t165 - t158 * t265) * qJD(2);
t106 = -t155 * t265 - t158 * t165;
t324 = qJD(6) * t106 - t165 * t222 - (t155 * t162 + t158 * t263) * qJD(2);
t202 = -mrSges(7,1) * t165 + mrSges(7,2) * t162;
t180 = m(7) * pkin(5) - t202;
t323 = mrSges(6,1) + t180;
t236 = m(7) * pkin(9) + mrSges(7,3);
t322 = mrSges(6,2) - t236;
t60 = qJD(6) * t117 + qJDD(5) * t162 + t126 * t165;
t61 = qJD(6) * t118 + qJDD(5) * t165 - t126 * t162;
t16 = -mrSges(7,1) * t61 + mrSges(7,2) * t60;
t279 = -qJDD(5) * mrSges(6,1) + mrSges(6,3) * t126 + t16;
t77 = -qJDD(2) * pkin(2) + t184;
t321 = qJD(2) * t124 - t77;
t289 = Ifges(6,4) * t166;
t290 = Ifges(6,4) * t163;
t320 = t163 * (-Ifges(6,1) * t166 + t290) + t166 * (Ifges(6,2) * t163 - t289);
t277 = qJD(5) * t37;
t8 = t142 * t166 - t163 * t21 - t277;
t318 = (t255 - t81) * t163;
t270 = t157 * t167;
t272 = t157 * t164;
t317 = t155 * t270 - t158 * t272;
t200 = -t163 * Ifges(6,1) - t289;
t114 = Ifges(7,4) * t117;
t57 = -Ifges(7,1) * t118 + Ifges(7,5) * t147 + t114;
t314 = Ifges(6,5) * qJD(5) + qJD(2) * t200 + t165 * t57;
t115 = qJDD(6) - t125;
t38 = mrSges(7,1) * t115 - mrSges(7,3) * t60;
t39 = -mrSges(7,2) * t115 + mrSges(7,3) * t61;
t313 = -t162 * t38 + t165 * t39;
t312 = -t163 * t8 + t166 * t7;
t203 = mrSges(6,1) * t166 - mrSges(6,2) * t163;
t311 = -t163 * t236 - t166 * t180 - mrSges(5,1) - t203;
t197 = -Ifges(6,2) * t166 - t290;
t310 = -Ifges(6,6) * qJD(5) / 0.2e1 + t197 * t294 + Ifges(7,5) * t301 + Ifges(7,6) * t338 + Ifges(7,3) * t337;
t201 = t162 * mrSges(7,1) + t165 * mrSges(7,2);
t286 = t118 * Ifges(7,4);
t56 = t117 * Ifges(7,2) + t147 * Ifges(7,6) - t286;
t308 = t201 * t28 - t162 * t56 / 0.2e1;
t307 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t306 = t342 * pkin(8) - mrSges(5,2) + mrSges(6,3) + t201;
t169 = qJD(2) ^ 2;
t305 = t60 / 0.2e1;
t304 = t61 / 0.2e1;
t303 = t115 / 0.2e1;
t300 = t118 / 0.2e1;
t6 = -qJDD(5) * pkin(5) - t8;
t297 = t163 * t6;
t293 = qJD(6) / 0.2e1;
t288 = Ifges(7,4) * t162;
t287 = Ifges(7,4) * t165;
t285 = t155 * t53;
t280 = t166 * t37;
t276 = (-qJD(2) * pkin(2) + t191) * t164;
t273 = t157 * t163;
t271 = t157 * t166;
t269 = t158 * t163;
t268 = t160 * t164;
t267 = t160 * t167;
t266 = t162 * t163;
t264 = t163 * t165;
t234 = mrSges(6,3) * t257;
t133 = -qJD(5) * mrSges(6,2) - t234;
t262 = t166 * t133;
t261 = pkin(2) * t270 + qJ(3) * t272;
t250 = qJD(6) * t163;
t248 = qJDD(2) * mrSges(4,1);
t247 = qJDD(2) * mrSges(5,1);
t246 = qJDD(2) * mrSges(5,2);
t244 = m(5) + t342;
t238 = Ifges(7,5) * t60 + Ifges(7,6) * t61 + Ifges(7,3) * t115;
t229 = pkin(3) * t270 + t261;
t221 = t166 * t245;
t216 = t249 / 0.2e1;
t215 = -t242 / 0.2e1;
t156 = sin(pkin(10));
t159 = cos(pkin(10));
t102 = t156 * t164 - t159 * t267;
t103 = t156 * t167 + t159 * t268;
t214 = -t102 * pkin(2) + qJ(3) * t103;
t104 = t156 * t267 + t159 * t164;
t105 = -t156 * t268 + t159 * t167;
t213 = -t104 * pkin(2) + qJ(3) * t105;
t48 = t102 * t155 + t103 * t158;
t212 = -t102 * t158 + t103 * t155;
t52 = t104 * t155 + t105 * t158;
t211 = -t104 * t158 + t105 * t155;
t209 = -t102 * pkin(3) + t214;
t208 = -t104 * pkin(3) + t213;
t199 = Ifges(7,1) * t165 - t288;
t198 = Ifges(7,1) * t162 + t287;
t196 = -Ifges(7,2) * t162 + t287;
t195 = Ifges(7,2) * t165 + t288;
t194 = -Ifges(6,5) * t166 + Ifges(6,6) * t163;
t193 = Ifges(7,5) * t165 - Ifges(7,6) * t162;
t192 = Ifges(7,5) * t162 + Ifges(7,6) * t165;
t190 = -t155 * t62 + t158 * t63;
t66 = -t160 * t163 - t166 * t317;
t27 = t162 * t89 + t165 * t66;
t26 = -t162 * t66 + t165 * t89;
t189 = -t163 * t36 + t280;
t65 = t160 * t166 - t163 * t317;
t185 = t252 * t28 + t297;
t183 = t53 * (-mrSges(6,1) * t163 - mrSges(6,2) * t166);
t178 = -t162 * t250 + t165 * t252;
t177 = t163 * t249 + t221;
t176 = -g(1) * t104 - g(2) * t102 + g(3) * t270;
t175 = -Ifges(7,5) * t163 - t166 * t199;
t174 = -Ifges(7,6) * t163 - t166 * t196;
t173 = -Ifges(7,3) * t163 - t166 * t193;
t171 = (-t163 * t37 - t166 * t36) * qJD(5) + t312;
t122 = t206 * qJD(2);
t121 = t203 * qJD(2);
t99 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t125;
t80 = qJD(2) * t89;
t79 = t317 * qJD(2);
t76 = mrSges(7,1) * t147 + mrSges(7,3) * t118;
t75 = -mrSges(7,2) * t147 + mrSges(7,3) * t117;
t72 = -mrSges(6,1) * t125 + mrSges(6,2) * t126;
t42 = t120 * t263 + t162 * t82;
t41 = -t120 * t265 + t165 * t82;
t35 = -t156 * t273 + t166 * t52;
t33 = t159 * t273 + t166 * t48;
t25 = -qJD(5) * t65 + t166 * t80;
t24 = -t160 * t254 + t163 * t80 - t252 * t317;
t18 = t122 * t162 + t165 * t36;
t17 = t122 * t165 - t162 * t36;
t14 = t60 * Ifges(7,1) + t61 * Ifges(7,4) + t115 * Ifges(7,5);
t13 = t60 * Ifges(7,4) + t61 * Ifges(7,2) + t115 * Ifges(7,6);
t4 = qJD(6) * t26 + t162 * t79 + t165 * t25;
t3 = -qJD(6) * t27 - t162 * t25 + t165 * t79;
t11 = [t79 * t121 + t25 * t133 + t26 * t38 + t27 * t39 + t3 * t76 + t4 * t75 + t66 * t99 + t89 * t72 + t279 * t65 - t278 * t24 + (qJD(2) * t80 - qJDD(2) * t317) * mrSges(5,2) + (qJD(2) * t79 + qJDD(2) * t89) * mrSges(5,1) + (-m(2) - t244 - t339) * g(3) + m(5) * (-t22 * t89 - t23 * t317 - t62 * t79 + t63 * t80) + m(6) * (t20 * t89 - t24 * t36 + t25 * t37 + t53 * t79 - t65 * t8 + t66 * t7) + m(7) * (t1 * t27 + t10 * t4 + t2 * t26 + t24 * t28 + t3 * t9 + t6 * t65) + (m(4) * (qJD(2) * t276 + t164 * t74 + t167 * t321) + m(3) * (t164 * t98 + t167 * t97) + (-t164 * t292 - t167 * t328) * t169 - t340 * qJDD(2)) * t157 + m(2) * qJDD(1) + (t333 * t339 - t336) * t160; (t173 * t337 + t174 * t338 + t175 * t301 + t183 + t194 * qJD(5) / 0.2e1) * qJD(5) + (-m(4) * t214 - m(5) * t209 - t342 * (pkin(4) * t212 + t209) + t328 * t103 + t292 * t102 + t311 * t212 + t306 * t48) * g(2) + (-m(4) * t213 - m(5) * t208 - t342 * (pkin(4) * t211 + t208) + t328 * t105 + t292 * t104 + t311 * t211 + t306 * t52) * g(1) + (-m(4) * t261 - m(5) * t229 - t342 * (t89 * pkin(4) + t229) + t340 * t157 - t306 * t317 + t311 * t89) * g(3) + t343 * t121 + (-Ifges(6,6) * qJDD(5) + Ifges(6,4) * t329 + Ifges(6,2) * t330 + Ifges(7,3) * t303 + Ifges(7,6) * t304 + Ifges(7,5) * t305 + t238 / 0.2e1 + t307) * t166 + (Ifges(6,1) * t329 + Ifges(6,4) * t330 - Ifges(6,5) * qJDD(5) + t120 * t279 - t193 * t303 - t196 * t304 - t199 * t305 + t56 * t216) * t163 + (t192 * t337 + t195 * t338 + t198 * t301) * t250 + t119 * t72 + (t1 * t266 + t10 * t177 + t178 * t9 + t2 * t264) * mrSges(7,3) + t99 * t274 + t125 * t197 / 0.2e1 + t126 * t200 / 0.2e1 + (-t189 * t81 - t53 * t78 + t119 * t20 + (t158 * t189 + t285) * qJD(3) + t171 * t120) * m(6) + t326 * t75 + (t1 * t42 + t10 * t326 + t120 * t185 + t2 * t41 + t28 * t318 + t327 * t9) * m(7) + t327 * t76 + t41 * t38 + t42 * t39 + t320 * t215 + (t252 * t36 + t254 * t37 - t312) * mrSges(6,3) - t314 * t252 / 0.2e1 + (t74 - t139 + t315) * mrSges(4,3) + t319 * t133 + (-t9 * mrSges(7,1) + t10 * mrSges(7,2) - t310) * t254 + t28 * (-mrSges(7,1) * t177 - mrSges(7,2) * t178) + t20 * t203 + (t97 + t138) * mrSges(3,1) + (qJD(3) * t190 + t127 * t22 + t128 * t23 + t62 * t78 - t63 * t81) * m(5) + (-t77 + t138) * mrSges(4,1) + (-qJD(2) * t78 + t155 * t243 - t22) * mrSges(5,1) + (-qJD(2) * t81 + t158 * t243 + t23) * mrSges(5,2) + (Ifges(3,3) + Ifges(5,3) + Ifges(4,2)) * qJDD(2) + t278 * (-t120 * t252 - t318) + t128 * t246 + pkin(2) * t248 + (-t98 + t139) * mrSges(3,2) + (-pkin(2) * t77 + qJ(3) * t74 + qJD(3) * t124 - (t124 * t167 + t276) * t260) * m(4) - t201 * t297 - t14 * t264 / 0.2e1 - t81 * t262 - t127 * t247 + t56 * t221 / 0.2e1 + (qJD(6) * t57 + t13) * t266 / 0.2e1; -t248 - t169 * mrSges(4,3) + t106 * t38 + t107 * t39 + t325 * t76 + t324 * t75 + (-t169 * mrSges(5,2) - t247 - t72) * t158 + (-t169 * mrSges(5,1) + t246 + t166 * t99 + t279 * t163 + (-t133 * t163 - t166 * t278) * qJD(5)) * t155 + (t155 * t171 - t158 * t20 + t176) * m(6) + (t155 * t23 + t158 * t22 + t176) * m(5) + (t176 - t321) * m(4) + (t1 * t107 + t10 * t324 + t106 * t2 + t155 * t185 + t325 * t9 + t176) * m(7) + ((t163 * t278 - t262) * t158 - t121 * t155 - m(6) * (t158 * t280 - t269 * t36 + t285) - t269 * t331 - m(5) * t190) * qJD(2); t336 + ((-t162 * t76 + t165 * t75 + t133) * qJD(5) + m(7) * (t10 * t253 - t245 * t9 - t6) + m(6) * (t8 + t277) - t279) * t166 + (t99 + (-t162 * t75 - t165 * t76) * qJD(6) - t278 * qJD(5) + m(7) * (qJD(5) * t28 + t334) + m(6) * (-qJD(5) * t36 + t7) + t313) * t163 + (t160 * g(3) + (g(1) * t156 - g(2) * t159) * t157) * t244; t165 * t13 / 0.2e1 + t334 * mrSges(7,3) + t162 * t14 / 0.2e1 + Ifges(6,6) * t125 + Ifges(6,5) * t126 + t192 * t303 + t195 * t304 + t198 * t305 + t6 * t202 + t194 * t215 + (-t235 + t278 - t331) * t37 - t18 * t75 - t17 * t76 + (t322 * t66 + t323 * t65) * g(3) + (t322 * t33 - t323 * (t159 * t271 - t163 * t48)) * g(2) + (t322 * t35 - t323 * (-t156 * t271 - t163 * t52)) * g(1) + t320 * t169 / 0.2e1 - pkin(5) * t16 + (m(7) * ((-t10 * t162 - t165 * t9) * qJD(6) + t205) - t76 * t249 - t75 * t251 + t313) * pkin(9) + t310 * t258 - t7 * mrSges(6,2) + t8 * mrSges(6,1) + (t199 * t301 + t308) * qJD(6) + (-t133 - t234) * t36 + (t173 * t294 + t193 * t293) * t147 + (-pkin(5) * t6 - t10 * t18 - t17 * t9) * m(7) + (t174 * t294 + t196 * t293) * t117 + Ifges(6,3) * qJDD(5) + (t175 * t300 - t183 - t9 * (-mrSges(7,1) * t163 + mrSges(7,3) * t263) - t10 * (mrSges(7,2) * t163 + mrSges(7,3) * t265)) * qJD(2) + (t314 / 0.2e1 + t308) * t257 + t57 * t216; -t28 * (-mrSges(7,1) * t118 + mrSges(7,2) * t117) + (Ifges(7,1) * t117 + t286) * t300 + t56 * t301 - t147 * (Ifges(7,5) * t117 + Ifges(7,6) * t118) / 0.2e1 - t9 * t75 + t10 * t76 - g(1) * ((-t162 * t35 + t165 * t211) * mrSges(7,1) + (-t162 * t211 - t165 * t35) * mrSges(7,2)) - g(2) * ((-t162 * t33 + t165 * t212) * mrSges(7,1) + (-t162 * t212 - t165 * t33) * mrSges(7,2)) - g(3) * (mrSges(7,1) * t26 - mrSges(7,2) * t27) + (-t10 * t118 + t117 * t9) * mrSges(7,3) + t238 - (Ifges(7,2) * t118 + t114 + t57) * t117 / 0.2e1 + t307;];
tau  = t11;
