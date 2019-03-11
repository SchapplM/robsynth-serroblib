% Calculate time derivative of joint inertia matrix for
% S6RRPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 11:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR12_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR12_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR12_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR12_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR12_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:17:20
% EndTime: 2019-03-09 11:17:32
% DurationCPUTime: 5.50s
% Computational Cost: add. (7141->550), mult. (17822->804), div. (0->0), fcn. (16739->10), ass. (0->253)
t322 = Ifges(5,3) + Ifges(6,3);
t187 = sin(qJ(6));
t190 = cos(qJ(6));
t188 = sin(qJ(4));
t192 = -pkin(2) - pkin(9);
t251 = qJ(5) - t192;
t153 = t251 * t188;
t184 = sin(pkin(11));
t285 = cos(qJ(4));
t235 = t285 * t192;
t196 = -qJ(5) * t285 + t235;
t264 = cos(pkin(11));
t111 = -t153 * t264 + t184 * t196;
t208 = t264 * t285;
t143 = t184 * t188 - t208;
t144 = t184 * t285 + t188 * t264;
t176 = t188 * pkin(4) + qJ(3);
t95 = pkin(5) * t144 + pkin(10) * t143 + t176;
t55 = -t111 * t187 + t190 * t95;
t129 = qJD(4) * t196 - t188 * qJD(5);
t246 = qJD(4) * t188;
t195 = -qJD(5) * t285 + t246 * t251;
t79 = t129 * t264 + t184 * t195;
t136 = -qJD(4) * t208 + t184 * t246;
t137 = t144 * qJD(4);
t220 = qJD(4) * t285;
t167 = pkin(4) * t220 + qJD(3);
t83 = -pkin(5) * t136 + pkin(10) * t137 + t167;
t24 = qJD(6) * t55 + t187 * t83 + t190 * t79;
t56 = t111 * t190 + t187 * t95;
t25 = -qJD(6) * t56 - t187 * t79 + t190 * t83;
t203 = -t187 * t25 + t190 * t24;
t244 = qJD(6) * t190;
t245 = qJD(6) * t187;
t321 = -t55 * t244 - t56 * t245 + t203;
t252 = t190 * t137;
t198 = t143 * t245 - t252;
t219 = -t245 / 0.2e1;
t283 = pkin(4) * t184;
t174 = pkin(10) + t283;
t320 = m(7) * t174;
t319 = mrSges(4,2) - mrSges(3,1);
t259 = t143 * t187;
t225 = t259 / 0.2e1;
t318 = qJD(6) * t225 - t252 / 0.2e1;
t218 = t244 / 0.2e1;
t253 = t187 * t137;
t317 = t143 * t218 + t253 / 0.2e1;
t262 = t136 * t187;
t316 = -t144 * t244 + t262;
t261 = t136 * t190;
t315 = t144 * t245 + t261;
t185 = sin(pkin(6));
t189 = sin(qJ(2));
t255 = t185 * t189;
t169 = pkin(8) * t255;
t186 = cos(pkin(6));
t191 = cos(qJ(2));
t284 = pkin(1) * t191;
t236 = -pkin(2) - t284;
t103 = pkin(3) * t255 + t169 + (-pkin(9) + t236) * t186;
t217 = -qJ(3) * t189 - pkin(1);
t117 = (t191 * t192 + t217) * t185;
t172 = t186 * t189 * pkin(1);
t254 = t185 * t191;
t291 = pkin(3) + pkin(8);
t118 = (t254 * t291 + t172) * qJD(2);
t233 = qJD(2) * t255;
t165 = pkin(2) * t233;
t247 = qJD(3) * t189;
t98 = t165 + (-t247 + (pkin(9) * t189 - qJ(3) * t191) * qJD(2)) * t185;
t34 = t103 * t220 - t117 * t246 + t188 * t118 + t285 * t98;
t66 = t188 * t103 + t285 * t117;
t35 = -qJD(4) * t66 + t285 * t118 - t188 * t98;
t313 = t188 * t34 + t285 * t35;
t312 = t136 * t184 + t264 * t137;
t197 = -t186 * t285 + t188 * t254;
t65 = t285 * t103 - t117 * t188;
t48 = pkin(4) * t255 + qJ(5) * t197 + t65;
t138 = -t186 * t188 - t254 * t285;
t51 = qJ(5) * t138 + t66;
t23 = t184 * t48 + t264 * t51;
t19 = pkin(10) * t255 + t23;
t141 = pkin(8) * t254 + t172;
t125 = -t186 * qJ(3) - t141;
t116 = pkin(3) * t254 - t125;
t80 = -pkin(4) * t138 + t116;
t86 = -t138 * t264 - t184 * t197;
t87 = t184 * t138 - t197 * t264;
t38 = pkin(5) * t86 - pkin(10) * t87 + t80;
t10 = -t187 * t19 + t190 * t38;
t114 = qJD(4) * t138 + t188 * t233;
t115 = qJD(4) * t197 + t233 * t285;
t63 = t114 * t184 - t115 * t264;
t64 = t114 * t264 + t184 * t115;
t241 = t186 * t284;
t166 = qJD(2) * t241;
t179 = t186 * qJD(3);
t102 = -t233 * t291 + t166 + t179;
t68 = -pkin(4) * t115 + t102;
t15 = pkin(5) * t63 - pkin(10) * t64 + t68;
t248 = qJD(2) * t191;
t232 = t185 * t248;
t17 = pkin(4) * t232 - t114 * qJ(5) + qJD(5) * t197 + t35;
t21 = qJ(5) * t115 + qJD(5) * t138 + t34;
t6 = t184 * t17 + t264 * t21;
t4 = pkin(10) * t232 + t6;
t1 = qJD(6) * t10 + t15 * t187 + t190 * t4;
t11 = t187 * t38 + t19 * t190;
t2 = -qJD(6) * t11 + t15 * t190 - t187 * t4;
t311 = t1 * t190 - t187 * t2;
t22 = -t184 * t51 + t264 * t48;
t310 = t136 * t23 + t137 * t22 - t144 * t6;
t215 = t136 * (-t187 ^ 2 - t190 ^ 2);
t18 = -pkin(5) * t255 - t22;
t71 = -t187 * t87 + t190 * t255;
t72 = t187 * t255 + t190 * t87;
t39 = -mrSges(7,1) * t71 + mrSges(7,2) * t72;
t77 = mrSges(6,1) * t255 - mrSges(6,3) * t87;
t309 = m(7) * t18 + t39 - t77;
t154 = -mrSges(7,1) * t190 + mrSges(7,2) * t187;
t234 = t264 * pkin(4);
t175 = -t234 - pkin(5);
t308 = m(7) * t175 - mrSges(6,1) + t154;
t120 = -mrSges(5,2) * t255 + mrSges(5,3) * t138;
t85 = -mrSges(5,2) * t232 + mrSges(5,3) * t115;
t307 = m(5) * ((-t188 * t65 + t285 * t66) * qJD(4) + t313) + t120 * t220 + t188 * t85;
t36 = qJD(6) * t71 + t187 * t232 + t190 * t64;
t37 = -qJD(6) * t72 - t187 * t64 + t190 * t232;
t12 = -mrSges(7,1) * t37 + mrSges(7,2) * t36;
t5 = t17 * t264 - t184 * t21;
t3 = -pkin(5) * t232 - t5;
t54 = mrSges(6,1) * t232 - mrSges(6,3) * t64;
t306 = -m(6) * t5 + m(7) * t3 + t12 - t54;
t305 = 0.2e1 * m(6);
t304 = 0.2e1 * m(7);
t303 = -0.2e1 * pkin(1);
t302 = 2 * mrSges(4,1);
t301 = -2 * mrSges(3,3);
t300 = -0.2e1 * mrSges(6,3);
t299 = 0.2e1 * mrSges(6,3);
t126 = (-pkin(2) * t191 + t217) * t185;
t298 = -0.2e1 * t126;
t297 = m(6) * pkin(4);
t296 = t36 / 0.2e1;
t295 = t37 / 0.2e1;
t294 = t63 / 0.2e1;
t293 = t71 / 0.2e1;
t292 = t72 / 0.2e1;
t289 = -t136 / 0.2e1;
t177 = Ifges(7,5) * t244;
t288 = Ifges(7,6) * t219 + t177 / 0.2e1;
t287 = t187 / 0.2e1;
t286 = t190 / 0.2e1;
t279 = Ifges(4,6) + Ifges(3,4);
t277 = Ifges(5,4) * t188;
t276 = Ifges(7,4) * t187;
t275 = Ifges(7,4) * t190;
t274 = Ifges(7,6) * t187;
t110 = -t153 * t184 - t196 * t264;
t78 = t129 * t184 - t195 * t264;
t273 = t110 * t78;
t130 = -pkin(8) * t233 + t166;
t272 = t130 * mrSges(3,2);
t269 = t143 * mrSges(6,3);
t260 = t137 * t143;
t258 = t143 * t190;
t257 = t174 * t187;
t256 = t174 * t190;
t250 = -Ifges(6,5) * t137 + Ifges(6,6) * t136;
t249 = Ifges(3,5) * t232 + Ifges(4,5) * t233;
t131 = t141 * qJD(2);
t243 = 0.2e1 * t131;
t242 = 0.2e1 * t191;
t7 = Ifges(7,5) * t36 + Ifges(7,6) * t37 + Ifges(7,3) * t63;
t240 = Ifges(5,4) * t285;
t231 = t192 * t246;
t228 = t174 * t245;
t227 = t174 * t244;
t224 = -t258 / 0.2e1;
t33 = t63 * mrSges(6,1) + t64 * mrSges(6,2);
t89 = -t136 * mrSges(6,1) - t137 * mrSges(6,2);
t214 = mrSges(4,1) * t232;
t213 = mrSges(5,2) * t220;
t211 = t232 / 0.2e1;
t207 = mrSges(7,1) * t187 + mrSges(7,2) * t190;
t206 = Ifges(7,1) * t190 - t276;
t205 = -Ifges(7,2) * t187 + t275;
t204 = t10 * t187 - t11 * t190;
t100 = mrSges(7,1) * t144 + mrSges(7,3) * t258;
t99 = -mrSges(7,2) * t144 + mrSges(7,3) * t259;
t202 = t100 * t187 - t190 * t99;
t201 = t110 * t137 + t143 * t78;
t200 = Ifges(5,5) * t114 + Ifges(6,5) * t64 + Ifges(5,6) * t115 - Ifges(6,6) * t63 + t232 * t322;
t199 = t143 * t244 + t253;
t194 = (-t10 * t190 - t11 * t187) * qJD(6) + t311;
t40 = Ifges(7,5) * t198 + Ifges(7,6) * t199 - Ifges(7,3) * t136;
t160 = Ifges(5,1) * t285 - t277;
t159 = Ifges(7,1) * t187 + t275;
t158 = -Ifges(5,2) * t188 + t240;
t157 = Ifges(7,2) * t190 + t276;
t156 = Ifges(7,5) * t187 + Ifges(7,6) * t190;
t155 = t188 * mrSges(5,1) + mrSges(5,2) * t285;
t152 = (-Ifges(5,1) * t188 - t240) * qJD(4);
t151 = t206 * qJD(6);
t150 = (-Ifges(5,2) * t285 - t277) * qJD(4);
t149 = t205 * qJD(6);
t147 = (mrSges(5,1) * t285 - mrSges(5,2) * t188) * qJD(4);
t146 = t207 * qJD(6);
t145 = -mrSges(4,1) * t254 - mrSges(4,3) * t186;
t140 = -t169 + t241;
t127 = t186 * t236 + t169;
t122 = -t130 - t179;
t121 = mrSges(5,1) * t255 + mrSges(5,3) * t197;
t119 = t165 + (-qJ(3) * t248 - t247) * t185;
t106 = -Ifges(6,1) * t143 - Ifges(6,4) * t144;
t105 = -Ifges(6,4) * t143 - Ifges(6,2) * t144;
t104 = mrSges(6,1) * t144 - mrSges(6,2) * t143;
t93 = -mrSges(5,1) * t138 - mrSges(5,2) * t197;
t92 = t207 * t143;
t91 = -Ifges(6,1) * t137 + Ifges(6,4) * t136;
t90 = -Ifges(6,4) * t137 + Ifges(6,2) * t136;
t84 = mrSges(5,1) * t232 - mrSges(5,3) * t114;
t82 = -Ifges(5,1) * t197 + Ifges(5,4) * t138 + Ifges(5,5) * t255;
t81 = -Ifges(5,4) * t197 + Ifges(5,2) * t138 + Ifges(5,6) * t255;
t76 = -mrSges(6,2) * t255 - mrSges(6,3) * t86;
t75 = Ifges(7,5) * t144 - t143 * t206;
t74 = Ifges(7,6) * t144 - t143 * t205;
t73 = Ifges(7,3) * t144 + (-Ifges(7,5) * t190 + t274) * t143;
t70 = mrSges(7,2) * t136 + mrSges(7,3) * t199;
t69 = -mrSges(7,1) * t136 - mrSges(7,3) * t198;
t67 = -mrSges(5,1) * t115 + mrSges(5,2) * t114;
t58 = Ifges(5,1) * t114 + Ifges(5,4) * t115 + Ifges(5,5) * t232;
t57 = Ifges(5,4) * t114 + Ifges(5,2) * t115 + Ifges(5,6) * t232;
t53 = -mrSges(6,2) * t232 - mrSges(6,3) * t63;
t52 = -mrSges(7,1) * t199 + mrSges(7,2) * t198;
t50 = mrSges(6,1) * t86 + mrSges(6,2) * t87;
t47 = Ifges(6,1) * t87 - Ifges(6,4) * t86 + Ifges(6,5) * t255;
t46 = Ifges(6,4) * t87 - Ifges(6,2) * t86 + Ifges(6,6) * t255;
t44 = mrSges(7,1) * t86 - mrSges(7,3) * t72;
t43 = -mrSges(7,2) * t86 + mrSges(7,3) * t71;
t42 = Ifges(7,1) * t198 + Ifges(7,4) * t199 - Ifges(7,5) * t136;
t41 = Ifges(7,4) * t198 + Ifges(7,2) * t199 - Ifges(7,6) * t136;
t30 = Ifges(6,1) * t64 - Ifges(6,4) * t63 + Ifges(6,5) * t232;
t29 = Ifges(6,4) * t64 - Ifges(6,2) * t63 + Ifges(6,6) * t232;
t28 = Ifges(7,1) * t72 + Ifges(7,4) * t71 + Ifges(7,5) * t86;
t27 = Ifges(7,4) * t72 + Ifges(7,2) * t71 + Ifges(7,6) * t86;
t26 = Ifges(7,5) * t72 + Ifges(7,6) * t71 + Ifges(7,3) * t86;
t14 = -mrSges(7,2) * t63 + mrSges(7,3) * t37;
t13 = mrSges(7,1) * t63 - mrSges(7,3) * t36;
t9 = Ifges(7,1) * t36 + Ifges(7,4) * t37 + Ifges(7,5) * t63;
t8 = Ifges(7,4) * t36 + Ifges(7,2) * t37 + Ifges(7,6) * t63;
t16 = [(t1 * t11 + t10 * t2 + t18 * t3) * t304 + (t22 * t5 + t23 * t6 + t68 * t80) * t305 + (t7 - t29) * t86 + 0.2e1 * m(3) * (t130 * t141 - t131 * t140) + 0.2e1 * m(4) * (t119 * t126 + t122 * t125 + t127 * t131) + 0.2e1 * m(5) * (t102 * t116 + t34 * t66 + t35 * t65) + ((t119 * mrSges(4,2) + t130 * mrSges(3,3)) * t242 + (-0.2e1 * t119 * mrSges(4,3) + (mrSges(4,1) + mrSges(3,3)) * t243 + t200) * t189 + ((t125 * t302 + mrSges(4,2) * t298 + t141 * t301 + (Ifges(4,5) - (2 * Ifges(3,6))) * t186 + (mrSges(3,1) * t303 - 0.2e1 * t189 * t279) * t185) * t189 + (t127 * t302 + t140 * t301 + mrSges(4,3) * t298 - Ifges(5,5) * t197 + Ifges(6,5) * t87 + Ifges(5,6) * t138 - Ifges(6,6) * t86 + (-(2 * Ifges(4,4)) + Ifges(3,5)) * t186 + (mrSges(3,2) * t303 + t242 * t279) * t185 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) + (2 * Ifges(4,2)) - (2 * Ifges(4,3)) + t322) * t255) * t191) * qJD(2)) * t185 + (t319 * t243 + t249 - 0.2e1 * t272) * t186 - t197 * t58 + 0.2e1 * t122 * t145 + t138 * t57 + 0.2e1 * t10 * t13 + 0.2e1 * t11 * t14 + (t26 - t46) * t63 + 0.2e1 * t18 * t12 + t36 * t28 + t37 * t27 + 0.2e1 * t3 * t39 + 0.2e1 * t1 * t43 + 0.2e1 * t2 * t44 + 0.2e1 * t23 * t53 + 0.2e1 * t22 * t54 + t64 * t47 + 0.2e1 * t68 * t50 + t71 * t8 + t72 * t9 + 0.2e1 * t6 * t76 + 0.2e1 * t5 * t77 + 0.2e1 * t80 * t33 + 0.2e1 * t65 * t84 + 0.2e1 * t66 * t85 + t87 * t30 + 0.2e1 * t102 * t93 + t114 * t82 + t115 * t81 + 0.2e1 * t116 * t67 + 0.2e1 * t34 * t120 + 0.2e1 * t35 * t121; t249 - t272 + m(6) * (t111 * t6 + t167 * t80 + t176 * t68 + t23 * t79) + m(7) * (t1 * t56 + t10 * t25 + t11 * t24 + t2 * t55) + (t250 + (-Ifges(5,5) * t188 - Ifges(5,6) * t285) * qJD(4)) * t255 / 0.2e1 + (t40 / 0.2e1 - t90 / 0.2e1) * t86 + t42 * t292 + t41 * t293 + t73 * t294 + t74 * t295 + t75 * t296 + t26 * t289 + t5 * t269 + t84 * t235 + (Ifges(5,5) * t285 - Ifges(5,6) * t188) * t211 + t285 * t58 / 0.2e1 - t82 * t246 / 0.2e1 + (t7 / 0.2e1 - t29 / 0.2e1 - Ifges(6,6) * t211) * t144 - Ifges(3,6) * t233 - t121 * t231 - Ifges(4,4) * t232 - t81 * t220 / 0.2e1 + (-t30 / 0.2e1 - Ifges(6,5) * t211) * t143 - pkin(2) * t214 + (-m(4) * t125 + m(5) * t116 - t145 + t93) * qJD(3) + (-m(4) * t122 + m(5) * t102 - mrSges(4,1) * t233 + t67) * qJ(3) + (-m(4) * pkin(2) + t319) * t131 + t317 * t27 + t318 * t28 - t197 * t152 / 0.2e1 + (-t220 * t66 + t246 * t65 - t313) * mrSges(5,3) - t188 * t57 / 0.2e1 + t167 * t50 + t176 * t33 + t310 * mrSges(6,3) + t307 * t192 + t138 * t150 / 0.2e1 + (-m(6) * t22 + t309) * t78 + t102 * t155 + t115 * t158 / 0.2e1 + t114 * t160 / 0.2e1 + t116 * t147 - t137 * t47 / 0.2e1 + t136 * t46 / 0.2e1 + t306 * t110 + t9 * t224 + t8 * t225 + t24 * t43 + t25 * t44 + t18 * t52 + t55 * t13 + t56 * t14 + t10 * t69 + t11 * t70 + t79 * t76 + t80 * t89 + t87 * t91 / 0.2e1 - t3 * t92 + t1 * t99 + t2 * t100 + t68 * t104 - t63 * t105 / 0.2e1 + t64 * t106 / 0.2e1 + t111 * t53 - t122 * mrSges(4,3); t285 * t152 + 0.2e1 * qJ(3) * t147 + 0.2e1 * t25 * t100 + 0.2e1 * t167 * t104 + 0.2e1 * t110 * t52 - t188 * t150 + 0.2e1 * t176 * t89 + 0.2e1 * t24 * t99 + 0.2e1 * t55 * t69 + 0.2e1 * t56 * t70 - 0.2e1 * t78 * t92 + (-t158 * t285 - t188 * t160) * qJD(4) + (t111 * t79 + t167 * t176 + t273) * t305 + (t24 * t56 + t25 * t55 + t273) * t304 + (t300 * t79 + t40 - t90) * t144 + (t111 * t299 + t105 - t73) * t136 + (t110 * t300 + t187 * t74 - t190 * t75 - t106) * t137 + 0.2e1 * (mrSges(4,3) + t155 + (m(4) + m(5)) * qJ(3)) * qJD(3) + (t78 * t300 + t187 * t41 - t190 * t42 - t91 + (t187 * t75 + t190 * t74) * qJD(6)) * t143; -m(6) * t310 - t121 * t246 + t285 * t84 + m(4) * t131 + t214 + t316 * t44 - t315 * t43 + (m(7) * t204 - t76) * t136 + t306 * t143 + t309 * t137 + (m(7) * t194 - t13 * t187 + t14 * t190 + t53) * t144 + t307; t143 * t52 + (-t92 - 0.2e1 * t269) * t137 + t202 * t136 + m(7) * (-t56 * t261 + t55 * t262 + t201) + m(6) * (-t111 * t136 + t201) + (t136 * t299 - t187 * t69 + t190 * t70 + (-t190 * t100 - t187 * t99) * qJD(6) + m(7) * t321 + m(6) * t79) * t144; 0.2e1 * m(6) * (-t136 * t144 + t260) + 0.2e1 * m(7) * (t144 * t215 + t260); t28 * t218 + t27 * t219 + t200 + t151 * t292 + t149 * t293 + t156 * t294 + t157 * t295 + t159 * t296 + (t184 * t6 + t264 * t5) * t297 + t53 * t283 + t8 * t286 + t9 * t287 + t86 * t288 + t14 * t256 + t54 * t234 - t13 * t257 - t44 * t227 - t43 * t228 + m(7) * (t174 * t194 + t175 * t3) + (-t10 * t244 - t11 * t245 + t311) * mrSges(7,3) + t175 * t12 + t3 * t154 + t18 * t146 + t5 * mrSges(6,1) - t6 * mrSges(6,2) - t34 * mrSges(5,2) + t35 * mrSges(5,1); t250 + t75 * t218 + t74 * t219 + (t184 * t297 - mrSges(6,2)) * t79 + ((-t187 * t56 - t190 * t55) * qJD(6) + t203) * t320 + t41 * t286 + t42 * t287 + t144 * t288 + t156 * t289 + t70 * t256 - t69 * t257 - Ifges(5,5) * t246 - mrSges(5,1) * t231 - t100 * t227 - t99 * t228 - Ifges(5,6) * t220 - t192 * t213 + t312 * mrSges(6,3) * pkin(4) + t321 * mrSges(7,3) + t317 * t157 + t318 * t159 + t175 * t52 + (-t264 * t297 + t308) * t78 + t110 * t146 + t151 * t224 + t149 * t225; -mrSges(5,1) * t246 + t136 * mrSges(6,2) + t308 * t137 + t143 * t146 - t312 * t297 - t213 + (t320 + mrSges(7,3)) * t215; 0.2e1 * t146 * t175 + t149 * t190 + t151 * t187 + (-t157 * t187 + t159 * t190) * qJD(6); t190 * t13 + t187 * t14 + (-t187 * t44 + t190 * t43) * qJD(6) + m(7) * (-qJD(6) * t204 + t1 * t187 + t190 * t2) + m(6) * t68 + t33; t187 * t70 + t190 * t69 - t202 * qJD(6) + m(7) * (t187 * t24 + t190 * t25 + (-t187 * t55 + t190 * t56) * qJD(6)) + m(6) * t167 + t89; 0; 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t7; mrSges(7,1) * t25 - mrSges(7,2) * t24 + t40; mrSges(7,1) * t316 + mrSges(7,2) * t315; t177 + (t154 * t174 - t274) * qJD(6); -t146; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t16(1) t16(2) t16(4) t16(7) t16(11) t16(16); t16(2) t16(3) t16(5) t16(8) t16(12) t16(17); t16(4) t16(5) t16(6) t16(9) t16(13) t16(18); t16(7) t16(8) t16(9) t16(10) t16(14) t16(19); t16(11) t16(12) t16(13) t16(14) t16(15) t16(20); t16(16) t16(17) t16(18) t16(19) t16(20) t16(21);];
Mq  = res;
