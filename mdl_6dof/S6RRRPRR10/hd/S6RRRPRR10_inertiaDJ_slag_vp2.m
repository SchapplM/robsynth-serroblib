% Calculate time derivative of joint inertia matrix for
% S6RRRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR10_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR10_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR10_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR10_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR10_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:15:41
% EndTime: 2019-03-09 19:15:58
% DurationCPUTime: 7.64s
% Computational Cost: add. (7349->589), mult. (16542->850), div. (0->0), fcn. (14143->8), ass. (0->234)
t271 = Ifges(5,4) + Ifges(4,5);
t292 = -Ifges(5,2) - Ifges(4,3);
t204 = sin(qJ(2));
t203 = sin(qJ(3));
t208 = cos(qJ(2));
t248 = qJD(2) * t208;
t237 = t203 * t248;
t207 = cos(qJ(3));
t245 = qJD(3) * t207;
t212 = t204 * t245 + t237;
t247 = qJD(3) * t203;
t291 = Ifges(5,6) * t247 + t271 * t245;
t202 = sin(qJ(5));
t206 = cos(qJ(5));
t214 = t202 * t203 + t206 * t207;
t288 = qJD(3) - qJD(5);
t104 = t288 * t214;
t215 = t202 * t207 - t203 * t206;
t63 = t104 * t204 - t215 * t248;
t131 = t215 * t204;
t64 = t131 * t288 + t214 * t248;
t290 = Ifges(6,5) * t64 + Ifges(6,6) * t63;
t279 = pkin(8) - pkin(9);
t175 = t279 * t203;
t176 = t279 * t207;
t112 = t202 * t175 + t206 * t176;
t264 = Ifges(5,5) * t203;
t219 = Ifges(5,1) * t207 + t264;
t266 = Ifges(4,4) * t203;
t220 = Ifges(4,1) * t207 - t266;
t289 = (t219 + t220) * qJD(3);
t209 = -pkin(3) - pkin(4);
t166 = -qJ(4) * t202 + t206 * t209;
t197 = t203 * qJ(4);
t216 = -t207 * pkin(3) - t197;
t287 = qJD(5) + qJD(6);
t201 = sin(qJ(6));
t205 = cos(qJ(6));
t145 = -t201 * t202 + t205 * t206;
t102 = t287 * t145;
t147 = t201 * t206 + t202 * t205;
t103 = t287 * t147;
t230 = -t103 * mrSges(7,1) - t102 * mrSges(7,2);
t286 = -(mrSges(6,1) * t202 + mrSges(6,2) * t206) * qJD(5) + t230;
t236 = t207 * t248;
t249 = qJD(2) * t204;
t285 = -t212 * Ifges(5,6) - t271 * t236 + t292 * t249;
t284 = 2 * m(4);
t283 = 2 * m(6);
t282 = 2 * m(7);
t281 = -2 * pkin(1);
t280 = 2 * pkin(7);
t275 = pkin(7) * t203;
t162 = -pkin(5) + t166;
t167 = t206 * qJ(4) + t202 * t209;
t100 = t162 * t205 - t167 * t201;
t129 = t206 * qJD(4) + qJD(5) * t166;
t130 = -t202 * qJD(4) - qJD(5) * t167;
t46 = qJD(6) * t100 + t129 * t205 + t130 * t201;
t274 = t46 * mrSges(7,2);
t272 = qJD(2) / 0.2e1;
t270 = Ifges(6,3) + Ifges(7,3);
t132 = t214 * t204;
t78 = -t131 * t205 - t132 * t201;
t23 = qJD(6) * t78 + t201 * t63 + t205 * t64;
t79 = -t131 * t201 + t132 * t205;
t24 = -qJD(6) * t79 - t201 * t64 + t205 * t63;
t269 = Ifges(7,5) * t23 + Ifges(7,6) * t24;
t105 = t288 * t215;
t94 = t201 * t215 - t205 * t214;
t31 = qJD(6) * t94 + t104 * t205 - t105 * t201;
t95 = -t201 * t214 - t205 * t215;
t32 = -qJD(6) * t95 - t104 * t201 - t105 * t205;
t268 = Ifges(7,5) * t31 + Ifges(7,6) * t32;
t169 = -pkin(2) * t208 - pkin(8) * t204 - pkin(1);
t186 = t208 * t275;
t200 = t208 * pkin(3);
t90 = pkin(4) * t208 + t186 + t200 + (-pkin(9) * t204 - t169) * t207;
t255 = t207 * t208;
t187 = pkin(7) * t255;
t123 = t203 * t169 + t187;
t113 = -qJ(4) * t208 + t123;
t257 = t203 * t204;
t96 = pkin(9) * t257 + t113;
t49 = t202 * t90 + t206 * t96;
t267 = Ifges(6,5) * t104 - Ifges(6,6) * t105;
t265 = Ifges(4,4) * t207;
t263 = Ifges(5,5) * t207;
t262 = Ifges(4,6) * t207;
t261 = t129 * mrSges(6,2);
t260 = t130 * mrSges(6,1);
t259 = t208 * Ifges(4,6);
t256 = t204 * t207;
t126 = -Ifges(5,4) * t208 + t204 * t219;
t127 = -Ifges(4,5) * t208 + t204 * t220;
t254 = t126 + t127;
t165 = (pkin(2) * t204 - pkin(8) * t208) * qJD(2);
t253 = t203 * t165 + t169 * t245;
t244 = qJD(4) * t207;
t252 = qJ(4) * t236 + t204 * t244;
t251 = qJ(4) * t245 + t203 * qJD(4);
t246 = qJD(3) * t204;
t243 = qJD(5) * t202;
t242 = qJD(5) * t206;
t241 = qJD(6) * t201;
t240 = qJD(6) * t205;
t168 = -pkin(2) + t216;
t239 = t209 * t203;
t238 = -pkin(3) - t275;
t235 = t203 * t246;
t171 = -Ifges(5,3) * t207 + t264;
t172 = Ifges(4,2) * t207 + t266;
t233 = t171 / 0.2e1 - t172 / 0.2e1;
t173 = Ifges(5,1) * t203 - t263;
t174 = Ifges(4,1) * t203 + t265;
t232 = t173 / 0.2e1 + t174 / 0.2e1;
t101 = t162 * t201 + t167 * t205;
t47 = -qJD(6) * t101 - t129 * t201 + t130 * t205;
t42 = t47 * mrSges(7,1);
t231 = t42 - t274;
t48 = -t202 * t96 + t206 * t90;
t111 = t206 * t175 - t176 * t202;
t122 = t169 * t207 - t186;
t141 = t207 * pkin(4) - t168;
t229 = -pkin(7) + t239;
t217 = Ifges(5,3) * t203 + t263;
t124 = -Ifges(5,6) * t208 + t204 * t217;
t218 = -Ifges(4,2) * t203 + t265;
t125 = t204 * t218 - t259;
t228 = t124 - t125 + t259;
t227 = qJD(3) * t187 - t165 * t207 + t169 * t247;
t33 = pkin(5) * t208 - pkin(10) * t132 + t48;
t34 = -pkin(10) * t131 + t49;
t17 = -t201 * t34 + t205 * t33;
t50 = pkin(9) * t235 + (-pkin(9) * t255 + (-pkin(4) + t238) * t204) * qJD(2) + t227;
t189 = qJ(4) * t249;
t51 = t189 + (-pkin(7) * qJD(2) + pkin(9) * qJD(3)) * t256 + (-qJD(4) + (-pkin(7) * qJD(3) + pkin(9) * qJD(2)) * t203) * t208 + t253;
t13 = -qJD(5) * t49 - t202 * t51 + t206 * t50;
t7 = -pkin(5) * t249 - pkin(10) * t64 + t13;
t12 = t202 * t50 + t206 * t51 + t90 * t242 - t243 * t96;
t8 = pkin(10) * t63 + t12;
t2 = qJD(6) * t17 + t201 * t7 + t205 * t8;
t18 = t201 * t33 + t205 * t34;
t3 = -qJD(6) * t18 - t201 * t8 + t205 * t7;
t226 = -t3 * mrSges(7,1) + t2 * mrSges(7,2) - t269;
t163 = t279 * t247;
t164 = qJD(3) * t176;
t65 = -t206 * t163 + t202 * t164 + t175 * t242 - t176 * t243;
t37 = -pkin(10) * t105 + t65;
t66 = -qJD(5) * t112 + t163 * t202 + t206 * t164;
t38 = -pkin(10) * t104 + t66;
t86 = pkin(10) * t215 + t111;
t87 = -pkin(10) * t214 + t112;
t40 = -t201 * t87 + t205 * t86;
t10 = qJD(6) * t40 + t201 * t38 + t205 * t37;
t41 = t201 * t86 + t205 * t87;
t11 = -qJD(6) * t41 - t201 * t37 + t205 * t38;
t225 = -t11 * mrSges(7,1) + t10 * mrSges(7,2) - t268;
t224 = -t207 * mrSges(4,1) + t203 * mrSges(4,2);
t223 = mrSges(4,1) * t203 + mrSges(4,2) * t207;
t170 = -t207 * mrSges(5,1) - t203 * mrSges(5,3);
t222 = mrSges(5,1) * t203 - mrSges(5,3) * t207;
t121 = qJD(3) * t239 + t251;
t185 = qJ(4) * t256;
t110 = t204 * t229 + t185;
t213 = -t235 + t236;
t118 = mrSges(5,2) * t236 + (-mrSges(5,1) * qJD(2) - mrSges(5,2) * t247) * t204;
t211 = t13 * mrSges(6,1) - t12 * mrSges(6,2) - t226 + t290;
t210 = t66 * mrSges(6,1) - t65 * mrSges(6,2) - t225 + t267;
t73 = (-t207 * t249 - t208 * t247) * pkin(7) + t253;
t60 = (t207 * t209 - t197) * t246 + t229 * t248 + t252;
t161 = -mrSges(5,2) * t257 - mrSges(5,3) * t208;
t160 = mrSges(5,1) * t208 + mrSges(5,2) * t256;
t159 = -mrSges(4,1) * t208 - mrSges(4,3) * t256;
t158 = mrSges(4,2) * t208 - mrSges(4,3) * t257;
t155 = t218 * qJD(3);
t154 = t217 * qJD(3);
t153 = t223 * qJD(3);
t152 = t222 * qJD(3);
t142 = (-mrSges(7,1) * t201 - mrSges(7,2) * t205) * qJD(6) * pkin(5);
t138 = t222 * t204;
t133 = pkin(3) * t247 - t251;
t128 = -t185 + (pkin(3) * t203 + pkin(7)) * t204;
t120 = -mrSges(5,2) * t212 + mrSges(5,3) * t249;
t119 = -mrSges(4,2) * t249 - mrSges(4,3) * t212;
t117 = mrSges(4,1) * t249 - mrSges(4,3) * t213;
t116 = mrSges(6,1) * t208 - mrSges(6,3) * t132;
t115 = -mrSges(6,2) * t208 - mrSges(6,3) * t131;
t114 = -t122 + t200;
t109 = pkin(5) * t214 + t141;
t108 = -Ifges(6,1) * t215 - Ifges(6,4) * t214;
t107 = -Ifges(6,4) * t215 - Ifges(6,2) * t214;
t106 = mrSges(6,1) * t214 - mrSges(6,2) * t215;
t93 = mrSges(4,1) * t212 + mrSges(4,2) * t213;
t92 = mrSges(5,1) * t212 - mrSges(5,3) * t213;
t85 = mrSges(6,1) * t131 + mrSges(6,2) * t132;
t84 = -t174 * t246 + (Ifges(4,5) * t204 + t208 * t220) * qJD(2);
t83 = -t173 * t246 + (Ifges(5,4) * t204 + t208 * t219) * qJD(2);
t82 = -t172 * t246 + (Ifges(4,6) * t204 + t208 * t218) * qJD(2);
t81 = -t171 * t246 + (Ifges(5,6) * t204 + t208 * t217) * qJD(2);
t77 = Ifges(6,1) * t132 - Ifges(6,4) * t131 + Ifges(6,5) * t208;
t76 = Ifges(6,4) * t132 - Ifges(6,2) * t131 + Ifges(6,6) * t208;
t75 = pkin(5) * t131 + t110;
t74 = t249 * t275 - t227;
t72 = pkin(3) * t212 + pkin(7) * t248 + qJ(4) * t235 - t252;
t71 = mrSges(7,1) * t208 - mrSges(7,3) * t79;
t70 = -mrSges(7,2) * t208 + mrSges(7,3) * t78;
t69 = t238 * t249 + t227;
t68 = pkin(5) * t105 + t121;
t67 = -qJD(4) * t208 + t189 + t73;
t59 = Ifges(6,1) * t104 - Ifges(6,4) * t105;
t58 = Ifges(6,4) * t104 - Ifges(6,2) * t105;
t57 = mrSges(6,1) * t105 + mrSges(6,2) * t104;
t56 = -mrSges(6,1) * t249 - mrSges(6,3) * t64;
t55 = mrSges(6,2) * t249 + mrSges(6,3) * t63;
t54 = Ifges(7,1) * t95 + Ifges(7,4) * t94;
t53 = Ifges(7,4) * t95 + Ifges(7,2) * t94;
t52 = -mrSges(7,1) * t94 + mrSges(7,2) * t95;
t39 = -mrSges(7,1) * t78 + mrSges(7,2) * t79;
t36 = Ifges(7,1) * t79 + Ifges(7,4) * t78 + Ifges(7,5) * t208;
t35 = Ifges(7,4) * t79 + Ifges(7,2) * t78 + Ifges(7,6) * t208;
t28 = -mrSges(6,1) * t63 + mrSges(6,2) * t64;
t27 = -pkin(5) * t63 + t60;
t26 = Ifges(6,1) * t64 + Ifges(6,4) * t63 - Ifges(6,5) * t249;
t25 = Ifges(6,4) * t64 + Ifges(6,2) * t63 - Ifges(6,6) * t249;
t20 = mrSges(7,2) * t249 + mrSges(7,3) * t24;
t19 = -mrSges(7,1) * t249 - mrSges(7,3) * t23;
t16 = Ifges(7,1) * t31 + Ifges(7,4) * t32;
t15 = Ifges(7,4) * t31 + Ifges(7,2) * t32;
t14 = -mrSges(7,1) * t32 + mrSges(7,2) * t31;
t6 = -mrSges(7,1) * t24 + mrSges(7,2) * t23;
t5 = Ifges(7,1) * t23 + Ifges(7,4) * t24 - Ifges(7,5) * t249;
t4 = Ifges(7,4) * t23 + Ifges(7,2) * t24 - Ifges(7,6) * t249;
t1 = [0.2e1 * m(5) * (t113 * t67 + t114 * t69 + t128 * t72) + (t93 * t280 + (t83 + t84) * t207 + (t81 - t82) * t203 + (t228 * t207 + (t208 * t271 - t254) * t203) * qJD(3) + (-Ifges(6,5) * t132 + Ifges(6,6) * t131 - Ifges(7,5) * t79 - Ifges(7,6) * t78 + (mrSges(3,1) * t281) + (-(2 * Ifges(3,4)) + t271 * t207 + (-Ifges(4,6) + Ifges(5,6)) * t203) * t204 + ((pkin(7) ^ 2 * t284) + t223 * t280 + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) - (2 * Ifges(6,3)) - (2 * Ifges(7,3)) + t292) * t208) * qJD(2)) * t204 + (t17 * t3 + t18 * t2 + t27 * t75) * t282 + (t110 * t60 + t12 * t49 + t13 * t48) * t283 + (t122 * t74 + t123 * t73) * t284 + (((mrSges(3,2) * t281) + 0.2e1 * Ifges(3,4) * t208 + t203 * t228 + t207 * t254) * qJD(2) + t269 + t285 + t290) * t208 + 0.2e1 * t17 * t19 + 0.2e1 * t18 * t20 + t24 * t35 + t23 * t36 + 0.2e1 * t27 * t39 + 0.2e1 * t49 * t55 + 0.2e1 * t48 * t56 + 0.2e1 * t2 * t70 + 0.2e1 * t3 * t71 + 0.2e1 * t75 * t6 + t63 * t76 + t64 * t77 + t78 * t4 + t79 * t5 + 0.2e1 * t60 * t85 + 0.2e1 * t110 * t28 + 0.2e1 * t12 * t115 + 0.2e1 * t13 * t116 + 0.2e1 * t114 * t118 + 0.2e1 * t113 * t120 + 0.2e1 * t122 * t117 + 0.2e1 * t123 * t119 + 0.2e1 * t128 * t92 - t131 * t25 + t132 * t26 + 0.2e1 * t72 * t138 + 0.2e1 * t73 * t158 + 0.2e1 * t74 * t159 + 0.2e1 * t69 * t160 + 0.2e1 * t67 * t161; (-t104 * t48 - t105 * t49 - t12 * t214 + t13 * t215) * mrSges(6,3) - t214 * t25 / 0.2e1 + (t262 * t272 - Ifges(3,6) * qJD(2) + (mrSges(3,2) * qJD(2) + t153) * pkin(7) - (-Ifges(6,5) * t215 + Ifges(7,5) * t95 - Ifges(6,6) * t214 + Ifges(7,6) * t94) * qJD(2) / 0.2e1 + (t154 / 0.2e1 - t155 / 0.2e1 - t232 * qJD(3) + t271 * t272) * t203 + (t289 / 0.2e1 - Ifges(5,6) * t272 + qJD(3) * t233) * t207) * t204 - t215 * t26 / 0.2e1 + (t268 + t267) * t208 / 0.2e1 + m(6) * (t110 * t121 + t111 * t13 + t112 * t12 + t141 * t60 + t48 * t66 + t49 * t65) + m(7) * (t10 * t18 + t109 * t27 + t11 * t17 + t2 * t41 + t3 * t40 + t68 * t75) - t291 * t208 / 0.2e1 + m(5) * (t128 * t133 + t168 * t72) + (Ifges(3,5) + t232 * t207 + t233 * t203 + (-m(4) * pkin(2) - mrSges(3,1) + t224) * pkin(7)) * t248 + (t67 * mrSges(5,2) + t73 * mrSges(4,3) - t81 / 0.2e1 + t82 / 0.2e1 + (t126 / 0.2e1 + t127 / 0.2e1 + t114 * mrSges(5,2) - t122 * mrSges(4,3)) * qJD(3)) * t207 + (-t17 * t31 + t18 * t32 + t2 * t94 - t3 * t95) * mrSges(7,3) + t32 * t35 / 0.2e1 + t31 * t36 / 0.2e1 + t40 * t19 + t41 * t20 + t27 * t52 + t24 * t53 / 0.2e1 + t23 * t54 / 0.2e1 + t68 * t39 + t10 * t70 + t11 * t71 + t75 * t14 + t78 * t15 / 0.2e1 + t79 * t16 / 0.2e1 + ((t119 + t120) * t207 + (-t117 + t118) * t203 + ((-t159 + t160) * t207 + (-t158 - t161) * t203) * qJD(3) + m(4) * (-t122 * t245 - t123 * t247 - t203 * t74 + t207 * t73) + m(5) * (-t113 * t247 + t114 * t245 + t203 * t69 + t207 * t67)) * pkin(8) - pkin(2) * t93 + t94 * t4 / 0.2e1 + t95 * t5 / 0.2e1 + t104 * t77 / 0.2e1 - t105 * t76 / 0.2e1 + t60 * t106 + t63 * t107 / 0.2e1 + t64 * t108 / 0.2e1 + t109 * t6 + t110 * t57 + t111 * t56 + t112 * t55 + t65 * t115 + t66 * t116 + t121 * t85 - t131 * t58 / 0.2e1 + t132 * t59 / 0.2e1 + (t69 * mrSges(5,2) - t74 * mrSges(4,3) + t83 / 0.2e1 + t84 / 0.2e1 + (-t113 * mrSges(5,2) - t123 * mrSges(4,3) - t125 / 0.2e1 + t124 / 0.2e1 + t259 / 0.2e1) * qJD(3)) * t203 + t133 * t138 + t141 * t28 + t128 * t152 + t168 * t92 + t72 * t170; -0.2e1 * pkin(2) * t153 + t104 * t108 - t105 * t107 + 0.2e1 * t121 * t106 + 0.2e1 * t109 * t14 + 0.2e1 * t141 * t57 - t214 * t58 - t215 * t59 + t94 * t15 + 0.2e1 * t168 * t152 + t95 * t16 + t31 * t54 + t32 * t53 + 0.2e1 * t68 * t52 + (-t154 + t155) * t207 + t289 * t203 + (t111 * t66 + t112 * t65 + t121 * t141) * t283 + (t10 * t41 + t109 * t68 + t11 * t40) * t282 + ((t173 + t174) * t207 + (t171 - t172) * t203) * qJD(3) + 0.2e1 * (m(5) * t168 + t170) * t133 + 0.2e1 * (t10 * t94 - t11 * t95 - t31 * t40 + t32 * t41) * mrSges(7,3) + 0.2e1 * (-t104 * t111 - t105 * t112 - t214 * t65 + t215 * t66) * mrSges(6,3); -t285 - t211 - Ifges(4,6) * t237 + t67 * mrSges(5,3) - t69 * mrSges(5,1) + t46 * t70 + t47 * t71 - t73 * mrSges(4,2) + t74 * mrSges(4,1) + t100 * t19 + t101 * t20 - pkin(3) * t118 + qJ(4) * t120 + t129 * t115 + t130 * t116 + qJD(4) * t161 + t166 * t56 + t167 * t55 + (t270 * qJD(2) + (-t203 * t271 - t262) * qJD(3)) * t204 + m(6) * (t12 * t167 + t129 * t49 + t13 * t166 + t130 * t48) + m(5) * (-pkin(3) * t69 + qJ(4) * t67 + qJD(4) * t113) + m(7) * (t100 * t3 + t101 * t2 + t17 * t47 + t18 * t46); (qJD(3) * t216 + t244) * mrSges(5,2) - t210 + m(6) * (t111 * t130 + t112 * t129 + t166 * t66 + t167 * t65) + m(7) * (t10 * t101 + t100 * t11 + t40 * t47 + t41 * t46) + (-t100 * t31 + t101 * t32 + t46 * t94 - t47 * t95) * mrSges(7,3) + (-t104 * t166 - t105 * t167 - t129 * t214 + t130 * t215) * mrSges(6,3) - Ifges(4,6) * t247 + (m(5) * t244 + (m(5) * t216 + t170 + t224) * qJD(3)) * pkin(8) + t291; -0.2e1 * t260 + 0.2e1 * t261 + 0.2e1 * t274 - 0.2e1 * t42 + 0.2e1 * (m(5) * qJ(4) + mrSges(5,3)) * qJD(4) + (t100 * t47 + t101 * t46) * t282 + (t129 * t167 + t130 * t166) * t283; t102 * t70 - t103 * t71 + t145 * t19 + t147 * t20 + t202 * t55 + t206 * t56 + (t115 * t206 - t116 * t202) * qJD(5) + m(7) * (t102 * t18 - t103 * t17 + t145 * t3 + t147 * t2) + m(6) * (t12 * t202 + t13 * t206 + (-t202 * t48 + t206 * t49) * qJD(5)) + m(5) * t69 + t118; (m(5) * pkin(8) + mrSges(5,2)) * t245 + m(7) * (t10 * t147 + t102 * t41 - t103 * t40 + t11 * t145) + m(6) * (t202 * t65 + t206 * t66 + (-t111 * t202 + t112 * t206) * qJD(5)) + (t102 * t94 + t103 * t95 - t145 * t31 + t147 * t32) * mrSges(7,3) + (-t206 * t104 - t202 * t105 + (-t202 * t215 - t206 * t214) * qJD(5)) * mrSges(6,3); m(7) * (-t100 * t103 + t101 * t102 + t145 * t47 + t147 * t46) + m(6) * (t129 * t202 + t130 * t206 + (-t166 * t202 + t167 * t206) * qJD(5)) - t286; (t102 * t147 - t103 * t145) * t282; -t270 * t249 + (m(7) * (-t17 * t241 + t18 * t240 + t2 * t201 + t205 * t3) + t70 * t240 + t201 * t20 - t71 * t241 + t205 * t19) * pkin(5) + t211; (m(7) * (t10 * t201 + t11 * t205 + (-t201 * t40 + t205 * t41) * qJD(6)) + (t201 * t32 - t205 * t31 + (t201 * t95 + t205 * t94) * qJD(6)) * mrSges(7,3)) * pkin(5) + t210; t260 - t261 + (m(7) * (-t100 * t241 + t101 * t240 + t201 * t46 + t205 * t47) + mrSges(7,2) * t240 + mrSges(7,1) * t241) * pkin(5) + t231; m(7) * (t102 * t201 - t103 * t205 + (-t145 * t201 + t147 * t205) * qJD(6)) * pkin(5) + t286; 0.2e1 * t142; -Ifges(7,3) * t249 - t226; -t225; t231; t230; t142; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
