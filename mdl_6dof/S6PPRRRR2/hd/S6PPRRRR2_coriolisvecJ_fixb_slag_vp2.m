% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PPRRRR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:03:22
% EndTime: 2019-03-08 19:03:36
% DurationCPUTime: 7.08s
% Computational Cost: add. (7569->521), mult. (19928->755), div. (0->0), fcn. (16419->14), ass. (0->258)
t187 = sin(qJ(5));
t191 = cos(qJ(5));
t192 = cos(qJ(4));
t246 = t191 * t192;
t185 = cos(pkin(6));
t172 = qJD(1) * t185 + qJD(2);
t181 = sin(pkin(7));
t184 = cos(pkin(7));
t193 = cos(qJ(3));
t183 = cos(pkin(13));
t182 = sin(pkin(6));
t241 = qJD(1) * t182;
t225 = t183 * t241;
t180 = sin(pkin(13));
t189 = sin(qJ(3));
t251 = t180 * t189;
t97 = t193 * (t172 * t181 + t184 * t225) - t241 * t251;
t248 = t183 * t184;
t199 = (t180 * t193 + t189 * t248) * t182;
t250 = t181 * t189;
t98 = qJD(1) * t199 + t172 * t250;
t54 = t187 * t98 + t246 * t97;
t188 = sin(qJ(4));
t218 = pkin(4) * t188 - pkin(10) * t192;
t162 = t218 * qJD(4);
t165 = -pkin(4) * t192 - pkin(10) * t188 - pkin(3);
t232 = qJD(5) * t191;
t233 = qJD(5) * t187;
t235 = qJD(4) * t191;
t88 = t187 * t162 + t165 * t232 + (-t188 * t235 - t192 * t233) * pkin(9);
t337 = t88 - t54;
t176 = pkin(9) * t246;
t208 = pkin(5) * t188 - pkin(11) * t246;
t236 = qJD(4) * t187;
t242 = t188 * pkin(9) * t236 + t191 * t162;
t282 = pkin(11) * t188;
t247 = t187 * t192;
t53 = t191 * t98 - t247 * t97;
t336 = -t53 + t208 * qJD(4) + (-t176 + (-t165 + t282) * t187) * qJD(5) + t242;
t234 = qJD(4) * t192;
t201 = t187 * t234 + t188 * t232;
t335 = pkin(11) * t201 - t337;
t304 = -pkin(11) - pkin(10);
t226 = qJD(5) * t304;
t238 = qJD(3) * t192;
t159 = t218 * qJD(3);
t131 = t172 * t184 - t181 * t225;
t96 = qJD(3) * pkin(9) + t98;
t64 = t131 * t192 - t188 * t96;
t51 = t187 * t159 + t191 * t64;
t334 = -t51 + (pkin(11) * t238 + t226) * t187;
t50 = t191 * t159 - t187 * t64;
t333 = -qJD(3) * t208 + t191 * t226 - t50;
t332 = -Ifges(5,1) / 0.2e1;
t179 = Ifges(5,4) * t238;
t331 = -t179 / 0.2e1;
t330 = qJD(4) / 0.2e1;
t190 = cos(qJ(6));
t186 = sin(qJ(6));
t240 = qJD(3) * t188;
t154 = -t187 * t240 + t235;
t254 = t131 * t188;
t65 = t192 * t96 + t254;
t63 = qJD(4) * pkin(10) + t65;
t87 = qJD(3) * t165 - t97;
t34 = t187 * t87 + t191 * t63;
t31 = pkin(11) * t154 + t34;
t259 = t186 * t31;
t175 = qJD(5) - t238;
t155 = t191 * t240 + t236;
t33 = -t187 * t63 + t191 * t87;
t30 = -pkin(11) * t155 + t33;
t29 = pkin(5) * t175 + t30;
t10 = t190 * t29 - t259;
t219 = t190 * t154 - t155 * t186;
t104 = Ifges(7,4) * t219;
t258 = t190 * t31;
t11 = t186 * t29 + t258;
t110 = t154 * t186 + t155 * t190;
t231 = qJD(3) * qJD(4);
t220 = t188 * t231;
t173 = Ifges(7,3) * t220;
t270 = Ifges(7,4) * t110;
t170 = qJD(6) + t175;
t288 = -t170 / 0.2e1;
t297 = -t110 / 0.2e1;
t299 = -t219 / 0.2e1;
t62 = -qJD(4) * pkin(4) - t64;
t52 = -pkin(5) * t154 + t62;
t57 = Ifges(7,1) * t110 + Ifges(7,5) * t170 + t104;
t329 = t173 + (Ifges(7,5) * t219 - Ifges(7,6) * t110) * t288 + (t10 * t219 + t11 * t110) * mrSges(7,3) + (-Ifges(7,2) * t110 + t104 + t57) * t299 - t52 * (mrSges(7,1) * t110 + mrSges(7,2) * t219) + (Ifges(7,1) * t219 - t270) * t297;
t222 = Ifges(5,5) * t330;
t153 = t191 * t165;
t111 = -t191 * t282 + t153 + (-pkin(9) * t187 - pkin(5)) * t192;
t133 = t187 * t165 + t176;
t121 = -t187 * t282 + t133;
t69 = t111 * t186 + t121 * t190;
t328 = -qJD(6) * t69 + t335 * t186 + t336 * t190;
t68 = t111 * t190 - t121 * t186;
t327 = qJD(6) * t68 + t336 * t186 - t335 * t190;
t168 = t304 * t187;
t169 = t304 * t191;
t124 = t168 * t190 + t169 * t186;
t326 = qJD(6) * t124 + t333 * t186 + t334 * t190;
t125 = t168 * t186 - t169 * t190;
t325 = -qJD(6) * t125 - t334 * t186 + t333 * t190;
t273 = Ifges(6,4) * t155;
t100 = Ifges(6,2) * t154 + Ifges(6,6) * t175 + t273;
t151 = Ifges(6,4) * t154;
t101 = Ifges(6,1) * t155 + Ifges(6,5) * t175 + t151;
t210 = t187 * t34 + t191 * t33;
t271 = Ifges(6,4) * t191;
t213 = -Ifges(6,2) * t187 + t271;
t272 = Ifges(6,4) * t187;
t215 = Ifges(6,1) * t191 - t272;
t216 = mrSges(6,1) * t187 + mrSges(6,2) * t191;
t268 = Ifges(6,6) * t187;
t269 = Ifges(6,5) * t191;
t284 = t191 / 0.2e1;
t285 = -t187 / 0.2e1;
t289 = t155 / 0.2e1;
t195 = -t210 * mrSges(6,3) + t101 * t284 + t100 * t285 + t62 * t216 + t154 * t213 / 0.2e1 + t215 * t289 + t175 * (-t268 + t269) / 0.2e1;
t95 = -qJD(3) * pkin(3) - t97;
t324 = -t95 * mrSges(5,2) + t64 * mrSges(5,3) + t240 * t332 - t195 - t222 + t331;
t200 = (t193 * t248 - t251) * t182;
t249 = t181 * t193;
t323 = t185 * t249 + t200;
t230 = qJD(4) * qJD(5);
t126 = t191 * t230 + (-t188 * t233 + t191 * t234) * qJD(3);
t93 = (qJD(1) * t200 + t172 * t249) * qJD(3);
t36 = qJD(4) * t64 + t192 * t93;
t82 = (t162 + t98) * qJD(3);
t9 = -qJD(5) * t34 - t187 * t36 + t191 * t82;
t6 = pkin(5) * t220 - pkin(11) * t126 + t9;
t127 = -qJD(3) * t201 - t187 * t230;
t8 = t187 * t82 + t191 * t36 + t87 * t232 - t233 * t63;
t7 = pkin(11) * t127 + t8;
t2 = qJD(6) * t10 + t186 * t6 + t190 * t7;
t3 = -qJD(6) * t11 - t186 * t7 + t190 * t6;
t48 = qJD(6) * t219 + t126 * t190 + t127 * t186;
t49 = -qJD(6) * t110 - t126 * t186 + t127 * t190;
t322 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t48 + Ifges(7,6) * t49;
t56 = Ifges(7,2) * t219 + Ifges(7,6) * t170 + t270;
t320 = t56 / 0.2e1;
t221 = -Ifges(5,6) * qJD(4) / 0.2e1;
t301 = m(7) * t52;
t66 = -mrSges(7,1) * t219 + mrSges(7,2) * t110;
t315 = t66 + t301;
t314 = qJD(4) * t65;
t209 = t186 * t187 - t190 * t191;
t138 = t209 * t188;
t313 = -t187 * t9 + t191 * t8;
t312 = -m(5) * t64 + m(6) * t62;
t311 = qJD(5) + qJD(6);
t308 = -t9 * mrSges(6,1) + t8 * mrSges(6,2) - Ifges(6,5) * t126 - Ifges(6,6) * t127 - t322;
t307 = 0.2e1 * m(5);
t306 = t48 / 0.2e1;
t305 = t49 / 0.2e1;
t298 = t219 / 0.2e1;
t296 = t110 / 0.2e1;
t295 = t126 / 0.2e1;
t294 = t127 / 0.2e1;
t157 = t186 * t191 + t187 * t190;
t137 = t157 * t188;
t293 = -t137 / 0.2e1;
t292 = -t138 / 0.2e1;
t291 = -t154 / 0.2e1;
t290 = -t155 / 0.2e1;
t287 = t170 / 0.2e1;
t286 = -t175 / 0.2e1;
t283 = pkin(5) * t187;
t37 = t188 * t93 + t314;
t113 = t185 * t250 + t199;
t139 = -t181 * t182 * t183 + t184 * t185;
t84 = t113 * t188 - t139 * t192;
t279 = t37 * t84;
t26 = -mrSges(7,1) * t49 + mrSges(7,2) * t48;
t86 = -mrSges(6,1) * t127 + mrSges(6,2) * t126;
t275 = t26 + t86;
t274 = Ifges(5,4) * t188;
t94 = qJD(3) * t98;
t265 = t323 * t94;
t141 = -t192 * t184 + t188 * t250;
t264 = t141 * t37;
t257 = t193 * t94;
t114 = t311 * t209;
t203 = t209 * t192;
t135 = qJD(3) * t203;
t245 = -t114 + t135;
t115 = t311 * t157;
t204 = t157 * t192;
t134 = qJD(3) * t204;
t244 = -t115 + t134;
t243 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t154 + mrSges(6,2) * t155 + mrSges(5,3) * t240;
t239 = qJD(3) * t189;
t237 = qJD(3) * t193;
t229 = -t66 - t243;
t228 = mrSges(5,3) * t238;
t224 = t181 * t239;
t223 = t181 * t237;
t217 = mrSges(6,1) * t191 - mrSges(6,2) * t187;
t214 = Ifges(6,1) * t187 + t271;
t212 = Ifges(6,2) * t191 + t272;
t211 = Ifges(6,5) * t187 + Ifges(6,6) * t191;
t85 = t113 * t192 + t139 * t188;
t44 = -t187 * t85 - t191 * t323;
t45 = -t187 * t323 + t191 * t85;
t22 = -t186 * t45 + t190 * t44;
t23 = t186 * t44 + t190 * t45;
t142 = t184 * t188 + t192 * t250;
t119 = -t142 * t187 - t191 * t249;
t207 = -t142 * t191 + t187 * t249;
t70 = t119 * t190 + t186 * t207;
t71 = t119 * t186 - t190 * t207;
t196 = t10 * mrSges(7,1) + t33 * mrSges(6,1) + t95 * mrSges(5,1) + t221 - (t192 * Ifges(5,2) + t274) * qJD(3) / 0.2e1 + t170 * Ifges(7,3) + t110 * Ifges(7,5) + t219 * Ifges(7,6) + t175 * Ifges(6,3) + t155 * Ifges(6,5) + t154 * Ifges(6,6) - t11 * mrSges(7,2) - t34 * mrSges(6,2) - t65 * mrSges(5,3);
t194 = qJD(3) ^ 2;
t178 = -pkin(5) * t191 - pkin(4);
t174 = Ifges(6,3) * t220;
t167 = -qJD(4) * mrSges(5,2) + t228;
t163 = (pkin(9) + t283) * t188;
t158 = (-mrSges(5,1) * t192 + mrSges(5,2) * t188) * qJD(3);
t149 = (mrSges(5,1) * t188 + mrSges(5,2) * t192) * t231;
t132 = -pkin(9) * t247 + t153;
t130 = pkin(5) * t201 + pkin(9) * t234;
t129 = mrSges(6,1) * t175 - mrSges(6,3) * t155;
t128 = -mrSges(6,2) * t175 + mrSges(6,3) * t154;
t118 = qJD(4) * t142 + t188 * t223;
t117 = -qJD(4) * t141 + t192 * t223;
t106 = -mrSges(6,2) * t220 + mrSges(6,3) * t127;
t105 = mrSges(6,1) * t220 - mrSges(6,3) * t126;
t103 = t113 * qJD(3);
t102 = t323 * qJD(3);
t91 = mrSges(7,1) * t170 - mrSges(7,3) * t110;
t90 = -mrSges(7,2) * t170 + mrSges(7,3) * t219;
t89 = -qJD(5) * t133 + t242;
t76 = t126 * Ifges(6,1) + t127 * Ifges(6,4) + Ifges(6,5) * t220;
t75 = t126 * Ifges(6,4) + t127 * Ifges(6,2) + Ifges(6,6) * t220;
t74 = -qJD(4) * t204 + t138 * t311;
t73 = -qJD(4) * t203 - t115 * t188;
t61 = qJD(5) * t207 - t117 * t187 + t191 * t224;
t60 = qJD(5) * t119 + t117 * t191 + t187 * t224;
t58 = t254 + (qJD(3) * t283 + t96) * t192;
t43 = -qJD(4) * t84 + t102 * t192;
t42 = qJD(4) * t85 + t102 * t188;
t40 = -mrSges(7,2) * t220 + mrSges(7,3) * t49;
t39 = mrSges(7,1) * t220 - mrSges(7,3) * t48;
t32 = -pkin(5) * t127 + t37;
t25 = t48 * Ifges(7,1) + t49 * Ifges(7,4) + Ifges(7,5) * t220;
t24 = t48 * Ifges(7,4) + t49 * Ifges(7,2) + Ifges(7,6) * t220;
t17 = -qJD(6) * t71 - t186 * t60 + t190 * t61;
t16 = qJD(6) * t70 + t186 * t61 + t190 * t60;
t15 = qJD(5) * t44 + t103 * t187 + t191 * t43;
t14 = -qJD(5) * t45 + t103 * t191 - t187 * t43;
t13 = t190 * t30 - t259;
t12 = -t186 * t30 - t258;
t5 = -qJD(6) * t23 + t14 * t190 - t15 * t186;
t4 = qJD(6) * t22 + t14 * t186 + t15 * t190;
t1 = [t103 * t158 + t44 * t105 + t45 * t106 - t323 * t149 + t15 * t128 + t14 * t129 + t43 * t167 + t22 * t39 + t23 * t40 + t4 * t90 + t5 * t91 + t275 * t84 - t229 * t42 + (-t103 * mrSges(4,1) - t102 * mrSges(4,2) + (-t188 * t85 + t192 * t84) * qJD(4) * mrSges(5,3)) * qJD(3) + m(7) * (t10 * t5 + t11 * t4 + t2 * t23 + t22 * t3 + t32 * t84 + t42 * t52) + m(4) * (t102 * t98 - t103 * t97 + t113 * t93 - t265) + m(5) * (t103 * t95 + t36 * t85 - t42 * t64 + t43 * t65 - t265 + t279) + m(6) * (t14 * t33 + t15 * t34 + t42 * t62 + t44 * t9 + t45 * t8 + t279); -t142 * mrSges(5,3) * t220 + t119 * t105 - t207 * t106 + t117 * t167 + t60 * t128 + t61 * t129 + t16 * t90 + t17 * t91 + t70 * t39 + t71 * t40 + (qJD(4) * t228 + t275) * t141 - t229 * t118 + m(7) * (t10 * t17 + t11 * t16 + t118 * t52 + t141 * t32 + t2 * t71 + t3 * t70) + m(5) * (t117 * t65 - t118 * t64 + t142 * t36 + t264) + m(6) * (t118 * t62 + t119 * t9 - t207 * t8 + t33 * t61 + t34 * t60 + t264) + ((-mrSges(4,2) * t194 - t149) * t193 + (-mrSges(4,1) * t194 + qJD(3) * t158) * t189 + m(4) * (t189 * t93 + t237 * t98 - t239 * t97 - t257) + m(5) * (t239 * t95 - t257)) * t181; (-t94 * mrSges(5,1) - t97 * t167 + t36 * mrSges(5,3) - t173 / 0.2e1 - t174 / 0.2e1 + (-t65 * t97 / 0.2e1 + pkin(9) * t36 / 0.2e1) * t307 + (0.3e1 / 0.2e1 * t179 + t222 + (t243 + t312) * pkin(9) - t324) * qJD(4) + t308) * t192 - m(6) * (t33 * t53 + t34 * t54) + (t76 * t284 + t75 * t285 + t94 * mrSges(5,2) + t215 * t295 + t213 * t294 + (mrSges(5,3) + t216) * t37 + (-t187 * t8 - t191 * t9) * mrSges(6,3) + (t62 * t217 + t212 * t291 + t214 * t290 + t211 * t286 - t191 * t100 / 0.2e1 + t101 * t285 + (t33 * t187 - t34 * t191) * mrSges(6,3)) * qJD(5) + (t196 + ((t269 / 0.2e1 - t268 / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(5,4)) * t188 + Ifges(7,5) * t292 + Ifges(7,6) * t293 + (-Ifges(6,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(5,2) - Ifges(7,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(5,1)) * t192) * qJD(3) + t221) * qJD(4) + (m(6) * t37 - t167 * qJD(4) + t86 + (t37 - t314) * m(5)) * pkin(9) + (t229 - t301 - t312) * t97) * t188 + m(6) * (t132 * t9 + t133 * t8 + t33 * t89 + t34 * t88) + t337 * t128 + (qJD(3) * t97 - t93) * mrSges(4,2) + (-Ifges(7,1) * t138 - Ifges(7,4) * t137) * t306 + t32 * (mrSges(7,1) * t137 - mrSges(7,2) * t138) + (-t10 * t73 + t11 * t74 - t137 * t2 + t138 * t3) * mrSges(7,3) + (-Ifges(7,4) * t138 - Ifges(7,2) * t137) * t305 + t327 * t90 + t328 * t91 + (t10 * t328 + t11 * t327 + t130 * t52 + t163 * t32 + t2 * t69 + t3 * t68) * m(7) + t68 * t39 + t69 * t40 + t73 * t57 / 0.2e1 + t52 * (-mrSges(7,1) * t74 + mrSges(7,2) * t73) + t130 * t66 + t132 * t105 + t133 * t106 + (t89 - t53) * t129 + (Ifges(7,5) * t73 + Ifges(7,6) * t74) * t287 + t25 * t292 + t24 * t293 + (Ifges(7,1) * t73 + Ifges(7,4) * t74) * t296 + (Ifges(7,4) * t73 + Ifges(7,2) * t74) * t298 + (-t95 * t98 / 0.2e1 - pkin(3) * t94 / 0.2e1) * t307 - pkin(3) * t149 - t98 * t158 + t163 * t26 + t74 * t320; (-mrSges(5,1) - t217) * t37 + (-t105 * t187 + t106 * t191) * pkin(10) - m(6) * (t33 * t50 + t34 * t51 + t62 * t65) + (t315 * t283 + (-m(6) * t210 - t187 * t128 - t191 * t129) * pkin(10) + t195) * qJD(5) + ((t331 + t222 + t324) * t192 + (-t196 + (t274 / 0.2e1 + (Ifges(5,2) / 0.2e1 + t332) * t192) * qJD(3) + t221 + (Ifges(7,5) * t157 - Ifges(7,6) * t209 + t211) * t330) * t188) * qJD(3) + t325 * t91 + t326 * t90 + (t10 * t325 + t11 * t326 + t124 * t3 + t125 * t2 + t178 * t32 - t52 * t58) * m(7) + (-Ifges(7,5) * t114 - Ifges(7,6) * t115) * t287 + (-Ifges(7,1) * t114 - Ifges(7,4) * t115) * t296 + (-Ifges(7,4) * t114 - Ifges(7,2) * t115) * t298 + (-t115 / 0.2e1 + t134 / 0.2e1) * t56 + (-Ifges(7,5) * t135 - Ifges(7,6) * t134) * t288 + (-Ifges(7,1) * t135 - Ifges(7,4) * t134) * t297 + m(6) * (-pkin(4) * t37 + pkin(10) * t313) + t313 * mrSges(6,3) + (-Ifges(7,4) * t135 - Ifges(7,2) * t134) * t299 + (-t114 / 0.2e1 + t135 / 0.2e1) * t57 + (Ifges(7,4) * t157 - Ifges(7,2) * t209) * t305 + (Ifges(7,1) * t157 - Ifges(7,4) * t209) * t306 + t32 * (mrSges(7,1) * t209 + mrSges(7,2) * t157) + (-t10 * t245 + t11 * t244 - t157 * t3 - t2 * t209) * mrSges(7,3) - t209 * t24 / 0.2e1 - t58 * t66 - pkin(4) * t86 + t124 * t39 + t125 * t40 - t51 * t128 - t50 * t129 - t36 * mrSges(5,2) + t75 * t284 + t212 * t294 + t214 * t295 + t157 * t25 / 0.2e1 - t64 * t167 + t178 * t26 + t187 * t76 / 0.2e1 - t243 * t65 + (-mrSges(7,1) * t244 + mrSges(7,2) * t245) * t52; (t186 * t40 + t190 * t39 + m(7) * (t186 * t2 + t190 * t3) - t315 * t155 + (-t186 * t91 + t190 * t90 + m(7) * (-t10 * t186 + t11 * t190)) * qJD(6)) * pkin(5) + (t154 * t33 + t155 * t34) * mrSges(6,3) - t308 + t110 * t320 - m(7) * (t10 * t12 + t11 * t13) + t174 - t13 * t90 - t12 * t91 - t33 * t128 + t34 * t129 + (Ifges(6,5) * t154 - Ifges(6,6) * t155) * t286 + t100 * t289 + (Ifges(6,1) * t154 - t273) * t290 - t62 * (mrSges(6,1) * t155 + mrSges(6,2) * t154) + (-Ifges(6,2) * t155 + t101 + t151) * t291 + t329; -t10 * t90 + t11 * t91 + t56 * t296 + t322 + t329;];
tauc  = t1(:);
