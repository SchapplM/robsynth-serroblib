% Calculate time derivative of joint inertia matrix for
% S6RRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 16:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:47:47
% EndTime: 2019-03-09 16:47:55
% DurationCPUTime: 5.38s
% Computational Cost: add. (7268->525), mult. (17653->747), div. (0->0), fcn. (15844->8), ass. (0->216)
t280 = -Ifges(4,3) - Ifges(5,3);
t278 = -mrSges(6,1) - mrSges(7,1);
t277 = Ifges(7,4) + Ifges(6,5);
t279 = -Ifges(7,2) - Ifges(6,3);
t276 = Ifges(7,6) - Ifges(6,6);
t207 = sin(qJ(3));
t208 = sin(qJ(2));
t209 = cos(qJ(3));
t241 = qJD(3) * t209;
t210 = cos(qJ(2));
t244 = qJD(2) * t210;
t217 = t207 * t244 + t208 * t241;
t205 = cos(pkin(10));
t193 = pkin(3) * t205 + pkin(4);
t206 = sin(qJ(5));
t204 = sin(pkin(10));
t266 = pkin(3) * t204;
t267 = cos(qJ(5));
t159 = t206 * t193 + t267 * t266;
t180 = -pkin(2) * t210 - t208 * pkin(8) - pkin(1);
t248 = t209 * t210;
t192 = pkin(7) * t248;
t142 = t207 * t180 + t192;
t167 = t204 * t209 + t205 * t207;
t161 = t167 * qJD(3);
t252 = t205 * t209;
t224 = t204 * t207 - t252;
t109 = -t161 * t208 - t224 * t244;
t245 = qJD(2) * t208;
t240 = qJD(4) * t209;
t177 = (pkin(2) * t208 - pkin(8) * t210) * qJD(2);
t265 = pkin(7) * t207;
t246 = t209 * t177 + t245 * t265;
t75 = -t208 * t240 + (pkin(3) * t208 - qJ(4) * t248) * qJD(2) + (-t192 + (qJ(4) * t208 - t180) * t207) * qJD(3) + t246;
t247 = t207 * t177 + t180 * t241;
t250 = t208 * t209;
t85 = (-qJD(2) * pkin(7) - qJ(4) * qJD(3)) * t250 + (-qJD(4) * t208 + (-pkin(7) * qJD(3) - qJ(4) * qJD(2)) * t210) * t207 + t247;
t28 = -t204 * t85 + t205 * t75;
t21 = pkin(4) * t245 - pkin(9) * t109 + t28;
t242 = qJD(3) * t208;
t108 = -t167 * t244 + t224 * t242;
t29 = t204 * t75 + t205 * t85;
t23 = pkin(9) * t108 + t29;
t153 = t224 * t208;
t169 = t209 * t180;
t124 = -qJ(4) * t250 + t169 + (-pkin(3) - t265) * t210;
t251 = t207 * t208;
t130 = -qJ(4) * t251 + t142;
t86 = t205 * t124 - t204 * t130;
t60 = -pkin(4) * t210 + t153 * pkin(9) + t86;
t152 = t167 * t208;
t87 = t204 * t124 + t205 * t130;
t65 = -pkin(9) * t152 + t87;
t262 = t206 * t60 + t267 * t65;
t6 = -qJD(5) * t262 - t206 * t23 + t21 * t267;
t236 = t209 * t244;
t275 = -Ifges(4,5) * t236 - Ifges(5,5) * t109 - Ifges(5,6) * t108 + t280 * t245;
t274 = 2 * m(4);
t273 = 2 * m(5);
t272 = 2 * m(6);
t271 = 2 * m(7);
t270 = -0.2e1 * pkin(1);
t269 = 0.2e1 * pkin(7);
t268 = -t207 / 0.2e1;
t263 = -qJ(4) - pkin(8);
t104 = -t206 * t152 - t153 * t267;
t93 = -mrSges(6,1) * t210 - t104 * mrSges(6,3);
t94 = mrSges(7,1) * t210 + t104 * mrSges(7,2);
t261 = t94 - t93;
t260 = Ifges(4,4) * t207;
t259 = Ifges(4,4) * t209;
t258 = Ifges(4,6) * t207;
t232 = qJD(5) * t267;
t238 = t206 * t266;
t146 = -qJD(5) * t238 + t193 * t232;
t140 = qJD(6) + t146;
t257 = t140 * mrSges(7,3);
t256 = t146 * mrSges(6,2);
t147 = t159 * qJD(5);
t181 = t263 * t207;
t170 = t204 * t181;
t182 = t263 * t209;
t132 = -t205 * t182 + t170;
t112 = -pkin(9) * t224 + t132;
t131 = t181 * t205 + t182 * t204;
t216 = -pkin(9) * t167 + t131;
t214 = t267 * t216;
t63 = t206 * t112 - t214;
t255 = t147 * t63;
t254 = t210 * Ifges(4,6);
t123 = t167 * t267 - t206 * t224;
t253 = t123 * t147;
t184 = Ifges(4,1) * t207 + t259;
t249 = t209 * t184;
t230 = qJD(3) * t263;
t111 = t205 * (t207 * t230 + t240) + t204 * (-qJD(4) * t207 + t209 * t230);
t178 = pkin(3) * t251 + t208 * pkin(7);
t243 = qJD(3) * t207;
t239 = qJD(5) * t206;
t201 = pkin(7) * t244;
t200 = pkin(3) * t243;
t138 = pkin(3) * t217 + t201;
t194 = -pkin(3) * t209 - pkin(2);
t110 = -t167 * qJD(4) + (t252 * t263 - t170) * qJD(3);
t162 = t224 * qJD(3);
t211 = t162 * pkin(9) + t110;
t90 = -pkin(9) * t161 + t111;
t17 = qJD(5) * t214 - t112 * t239 + t206 * t211 + t267 * t90;
t64 = t112 * t267 + t206 * t216;
t18 = qJD(5) * t64 + t206 * t90 - t211 * t267;
t237 = t64 * t17 + t18 * t63;
t235 = t207 * t242;
t220 = -t152 * t267 + t206 * t153;
t47 = qJD(5) * t220 + t206 * t108 + t109 * t267;
t48 = qJD(5) * t104 - t108 * t267 + t206 * t109;
t14 = t48 * mrSges(6,1) + t47 * mrSges(6,2);
t219 = -t206 * t167 - t224 * t267;
t76 = qJD(5) * t219 - t206 * t161 - t162 * t267;
t77 = qJD(5) * t123 + t161 * t267 - t206 * t162;
t31 = t77 * mrSges(6,1) + t76 * mrSges(6,2);
t13 = t48 * mrSges(7,1) - t47 * mrSges(7,3);
t30 = t77 * mrSges(7,1) - t76 * mrSges(7,3);
t231 = (2 * Ifges(3,4)) + t258;
t135 = pkin(4) * t161 + t200;
t66 = -t108 * mrSges(5,1) + t109 * mrSges(5,2);
t116 = t161 * mrSges(5,1) - t162 * mrSges(5,2);
t125 = pkin(4) * t152 + t178;
t38 = -mrSges(7,1) * t245 + t47 * mrSges(7,2);
t229 = -mrSges(4,1) * t209 + mrSges(4,2) * t207;
t228 = mrSges(4,1) * t207 + mrSges(4,2) * t209;
t227 = Ifges(4,1) * t209 - t260;
t226 = -Ifges(4,2) * t207 + t259;
t183 = Ifges(4,2) * t209 + t260;
t225 = Ifges(4,5) * t207 + Ifges(4,6) * t209;
t143 = pkin(4) * t224 + t194;
t88 = -pkin(4) * t108 + t138;
t223 = t245 * t279 - t276 * t48 - t277 * t47;
t26 = -t206 * t65 + t267 * t60;
t5 = t206 * t21 + t267 * t23 + t60 * t232 - t239 * t65;
t158 = t193 * t267 - t238;
t218 = -t235 + t236;
t71 = Ifges(7,6) * t77;
t72 = Ifges(6,6) * t77;
t73 = Ifges(6,5) * t76;
t74 = Ifges(7,4) * t76;
t215 = t71 - t72 + t73 + t74 + t278 * t18 + (-mrSges(6,2) + mrSges(7,3)) * t17;
t2 = qJ(6) * t245 - qJD(6) * t210 + t5;
t3 = -pkin(5) * t245 - t6;
t213 = t6 * mrSges(6,1) - t3 * mrSges(7,1) - t5 * mrSges(6,2) + t2 * mrSges(7,3) - t223;
t203 = qJD(6) * mrSges(7,3);
t199 = Ifges(4,5) * t241;
t176 = -mrSges(4,1) * t210 - mrSges(4,3) * t250;
t175 = mrSges(4,2) * t210 - mrSges(4,3) * t251;
t174 = t227 * qJD(3);
t173 = t226 * qJD(3);
t172 = t228 * qJD(3);
t157 = Ifges(5,5) * t162;
t156 = Ifges(5,6) * t161;
t154 = -pkin(5) - t158;
t151 = qJ(6) + t159;
t150 = -Ifges(4,5) * t210 + t208 * t227;
t149 = t208 * t226 - t254;
t141 = -t210 * t265 + t169;
t137 = -mrSges(4,2) * t245 - mrSges(4,3) * t217;
t136 = mrSges(4,1) * t245 - mrSges(4,3) * t218;
t134 = -mrSges(5,1) * t210 + t153 * mrSges(5,3);
t133 = mrSges(5,2) * t210 - t152 * mrSges(5,3);
t129 = Ifges(5,1) * t167 - Ifges(5,4) * t224;
t128 = Ifges(5,4) * t167 - Ifges(5,2) * t224;
t127 = mrSges(5,1) * t224 + mrSges(5,2) * t167;
t121 = mrSges(4,1) * t217 + mrSges(4,2) * t218;
t118 = -Ifges(5,1) * t162 - Ifges(5,4) * t161;
t117 = -Ifges(5,4) * t162 - Ifges(5,2) * t161;
t115 = -t184 * t242 + (Ifges(4,5) * t208 + t210 * t227) * qJD(2);
t114 = -t183 * t242 + (Ifges(4,6) * t208 + t210 * t226) * qJD(2);
t113 = mrSges(5,1) * t152 - mrSges(5,2) * t153;
t100 = -Ifges(5,1) * t153 - Ifges(5,4) * t152 - Ifges(5,5) * t210;
t99 = -Ifges(5,4) * t153 - Ifges(5,2) * t152 - Ifges(5,6) * t210;
t98 = -qJD(3) * t142 + t246;
t97 = (-t209 * t245 - t210 * t243) * pkin(7) + t247;
t96 = mrSges(5,1) * t245 - mrSges(5,3) * t109;
t95 = -mrSges(5,2) * t245 + mrSges(5,3) * t108;
t92 = mrSges(6,2) * t210 + mrSges(6,3) * t220;
t91 = mrSges(7,2) * t220 - mrSges(7,3) * t210;
t84 = Ifges(6,1) * t123 + Ifges(6,4) * t219;
t83 = Ifges(7,1) * t123 - Ifges(7,5) * t219;
t82 = Ifges(6,4) * t123 + Ifges(6,2) * t219;
t81 = Ifges(7,5) * t123 - Ifges(7,3) * t219;
t80 = -mrSges(6,1) * t219 + mrSges(6,2) * t123;
t79 = -mrSges(7,1) * t219 - mrSges(7,3) * t123;
t61 = -pkin(5) * t219 - qJ(6) * t123 + t143;
t59 = -mrSges(6,1) * t220 + mrSges(6,2) * t104;
t58 = -mrSges(7,1) * t220 - mrSges(7,3) * t104;
t57 = Ifges(5,1) * t109 + Ifges(5,4) * t108 + Ifges(5,5) * t245;
t56 = Ifges(5,4) * t109 + Ifges(5,2) * t108 + Ifges(5,6) * t245;
t53 = Ifges(6,1) * t104 + Ifges(6,4) * t220 - Ifges(6,5) * t210;
t52 = Ifges(7,1) * t104 - Ifges(7,4) * t210 - Ifges(7,5) * t220;
t51 = Ifges(6,4) * t104 + Ifges(6,2) * t220 - Ifges(6,6) * t210;
t50 = Ifges(7,5) * t104 - Ifges(7,6) * t210 - Ifges(7,3) * t220;
t49 = -pkin(5) * t220 - qJ(6) * t104 + t125;
t39 = -mrSges(6,2) * t245 - mrSges(6,3) * t48;
t37 = mrSges(6,1) * t245 - mrSges(6,3) * t47;
t36 = -mrSges(7,2) * t48 + mrSges(7,3) * t245;
t35 = Ifges(6,1) * t76 - Ifges(6,4) * t77;
t34 = Ifges(7,1) * t76 + Ifges(7,5) * t77;
t33 = Ifges(6,4) * t76 - Ifges(6,2) * t77;
t32 = Ifges(7,5) * t76 + Ifges(7,3) * t77;
t25 = t210 * pkin(5) - t26;
t24 = -qJ(6) * t210 + t262;
t19 = pkin(5) * t77 - qJ(6) * t76 - qJD(6) * t123 + t135;
t12 = Ifges(6,1) * t47 - Ifges(6,4) * t48 + Ifges(6,5) * t245;
t11 = Ifges(7,1) * t47 + Ifges(7,4) * t245 + Ifges(7,5) * t48;
t10 = Ifges(6,4) * t47 - Ifges(6,2) * t48 + Ifges(6,6) * t245;
t9 = Ifges(7,5) * t47 + Ifges(7,6) * t245 + Ifges(7,3) * t48;
t7 = pkin(5) * t48 - qJ(6) * t47 - qJD(6) * t104 + t88;
t1 = [(t50 - t51) * t48 + (t52 + t53) * t47 + ((mrSges(3,2) * t270 - t207 * t149 + t209 * t150 + t210 * t231) * qJD(2) + t223 + t275) * t210 - (t9 - t10) * t220 + (t121 * t269 - t207 * t114 + t209 * t115 + (-t209 * t149 - t207 * t150 + t210 * t225) * qJD(3) + (-Ifges(5,5) * t153 - Ifges(5,6) * t152 + mrSges(3,1) * t270 + (Ifges(4,5) * t209 - t231) * t208 + t277 * t104 - t276 * t220 + (pkin(7) ^ 2 * t274 + t228 * t269 + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) + t279 + t280) * t210) * qJD(2)) * t208 + (t125 * t88 + t26 * t6 + t262 * t5) * t272 + 0.2e1 * t262 * t39 + (t11 + t12) * t104 + (t2 * t24 + t25 * t3 + t49 * t7) * t271 + (t138 * t178 + t28 * t86 + t29 * t87) * t273 + (t141 * t98 + t142 * t97) * t274 + 0.2e1 * t24 * t36 + 0.2e1 * t26 * t37 + 0.2e1 * t25 * t38 + 0.2e1 * t49 * t13 + 0.2e1 * t7 * t58 + 0.2e1 * t88 * t59 + 0.2e1 * t2 * t91 + 0.2e1 * t5 * t92 + 0.2e1 * t6 * t93 + 0.2e1 * t3 * t94 + 0.2e1 * t87 * t95 + 0.2e1 * t86 * t96 + t108 * t99 + t109 * t100 + 0.2e1 * t125 * t14 + 0.2e1 * t29 * t133 + 0.2e1 * t28 * t134 + 0.2e1 * t138 * t113 + 0.2e1 * t141 * t136 + 0.2e1 * t142 * t137 - t152 * t56 - t153 * t57 + 0.2e1 * t97 * t175 + 0.2e1 * t98 * t176 + 0.2e1 * t178 * t66; (t11 / 0.2e1 + t12 / 0.2e1) * t123 + m(7) * (t17 * t24 + t18 * t25 + t19 * t49 + t2 * t64 + t3 * t63 + t61 * t7) + (t209 * t174 / 0.2e1 + t173 * t268 - Ifges(3,6) * qJD(2) + (-t209 * t183 / 0.2e1 + t184 * t268) * qJD(3) + (qJD(2) * mrSges(3,2) + t172) * pkin(7) + (Ifges(5,5) * t167 - Ifges(5,6) * t224 + t277 * t123 - t219 * t276 + t225) * qJD(2) / 0.2e1) * t208 + (-t161 * t87 + t162 * t86 - t167 * t28 - t224 * t29) * mrSges(5,3) - t224 * t56 / 0.2e1 + (-t199 / 0.2e1 - t73 / 0.2e1 + t72 / 0.2e1 - t74 / 0.2e1 - t71 / 0.2e1 + t157 / 0.2e1 + t156 / 0.2e1 + (Ifges(3,5) + t249 / 0.2e1 + t183 * t268 + (-mrSges(3,1) + t229) * pkin(7)) * qJD(2)) * t210 - (t32 / 0.2e1 - t33 / 0.2e1) * t220 - (t9 / 0.2e1 - t10 / 0.2e1) * t219 + (-t123 * t6 + t219 * t5 - t26 * t76 - t262 * t77) * mrSges(6,3) + (t123 * t3 + t2 * t219 - t24 * t77 + t25 * t76) * mrSges(7,2) + (pkin(8) * t137 + t97 * mrSges(4,3) + t114 / 0.2e1) * t209 + (t81 / 0.2e1 - t82 / 0.2e1) * t48 + (-pkin(8) * t136 - t98 * mrSges(4,3) + t115 / 0.2e1) * t207 + (t83 / 0.2e1 + t84 / 0.2e1) * t47 + t261 * t18 + ((-pkin(8) * t176 - t141 * mrSges(4,3) + t150 / 0.2e1) * t209 + (-pkin(8) * t175 + pkin(3) * t113 - t142 * mrSges(4,3) - t149 / 0.2e1 + t254 / 0.2e1) * t207) * qJD(3) + (t50 / 0.2e1 - t51 / 0.2e1) * t77 + m(5) * (t110 * t86 + t111 * t87 + t131 * t28 + t132 * t29 + t138 * t194 + t178 * t200) + (t36 + t39) * t64 + (t38 - t37) * t63 + (t91 + t92) * t17 + (t34 / 0.2e1 + t35 / 0.2e1) * t104 + (-pkin(2) * t201 + (-t207 * t98 + t209 * t97 + (-t141 * t209 - t142 * t207) * qJD(3)) * pkin(8)) * m(4) + (t52 / 0.2e1 + t53 / 0.2e1) * t76 + m(6) * (t125 * t135 + t143 * t88 + t17 * t262 - t18 * t26 + t5 * t64 - t6 * t63) + t49 * t30 + t19 * t58 + t61 * t13 + t7 * t79 + t88 * t80 - pkin(2) * t121 + t125 * t31 + t108 * t128 / 0.2e1 + t109 * t129 / 0.2e1 + t131 * t96 + t132 * t95 + t111 * t133 + t110 * t134 + t135 * t59 + t138 * t127 + t143 * t14 - t152 * t117 / 0.2e1 - t153 * t118 / 0.2e1 - t161 * t99 / 0.2e1 - t162 * t100 / 0.2e1 + t167 * t57 / 0.2e1 + t178 * t116 + t194 * t66; -0.2e1 * pkin(2) * t172 + 0.2e1 * t194 * t116 - t224 * t117 + t167 * t118 - t161 * t128 - t162 * t129 + 0.2e1 * t135 * t80 + 0.2e1 * t143 * t31 + t209 * t173 + t207 * t174 + 0.2e1 * t19 * t79 + 0.2e1 * t61 * t30 + (t81 - t82) * t77 + (t83 + t84) * t76 + (t34 + t35) * t123 - (t32 - t33) * t219 + (t249 + (0.2e1 * pkin(3) * t127 - t183) * t207) * qJD(3) + (t110 * t131 + t111 * t132 + t194 * t200) * t273 + (t135 * t143 + t237) * t272 + (t19 * t61 + t237) * t271 + 0.2e1 * (-t110 * t167 - t111 * t224 + t131 * t162 - t132 * t161) * mrSges(5,3) + 0.2e1 * (mrSges(6,3) + mrSges(7,2)) * (t123 * t18 + t17 * t219 + t63 * t76 - t64 * t77); -t275 + t213 + m(7) * (t140 * t24 + t147 * t25 + t151 * t2 + t154 * t3) + (t204 * t95 + m(5) * (t204 * t29 + t205 * t28) + t205 * t96) * pkin(3) + t261 * t147 - Ifges(4,5) * t235 - t217 * Ifges(4,6) + m(6) * (t146 * t262 - t147 * t26 + t158 * t6 + t159 * t5) + t28 * mrSges(5,1) - t29 * mrSges(5,2) - t97 * mrSges(4,2) + t98 * mrSges(4,1) + t140 * t91 + t146 * t92 + t151 * t36 + t154 * t38 + t158 * t37 + t159 * t39; t215 + (pkin(8) * t229 - t258) * qJD(3) + (m(5) * (t110 * t205 + t111 * t204) + (-t161 * t204 + t162 * t205) * mrSges(5,3)) * pkin(3) - t156 - t157 + t199 + m(7) * (t140 * t64 + t151 * t17 + t154 * t18 + t255) + m(6) * (t146 * t64 - t158 * t18 + t159 * t17 + t255) + (t146 * t219 - t158 * t76 - t159 * t77 + t253) * mrSges(6,3) + (t140 * t219 - t151 * t77 + t154 * t76 + t253) * mrSges(7,2) + t110 * mrSges(5,1) - t111 * mrSges(5,2); -0.2e1 * t256 + 0.2e1 * t257 + 0.2e1 * t278 * t147 + (t140 * t151 + t147 * t154) * t271 + (t146 * t159 - t147 * t158) * t272; m(5) * t138 + m(6) * t88 + m(7) * t7 + t13 + t14 + t66; m(5) * t200 + m(6) * t135 + m(7) * t19 + t116 + t30 + t31; 0; 0; t213 + m(7) * (-pkin(5) * t3 + qJ(6) * t2 + qJD(6) * t24) + qJ(6) * t36 - pkin(5) * t38 + qJD(6) * t91; m(7) * (-pkin(5) * t18 + qJ(6) * t17 + qJD(6) * t64) + (-pkin(5) * t76 - qJ(6) * t77 + qJD(6) * t219) * mrSges(7,2) + t215; t257 + m(7) * (qJ(6) * t140 + qJD(6) * t151) + t203 - t256 + (-m(7) * pkin(5) + t278) * t147; 0; 0.2e1 * m(7) * qJ(6) * qJD(6) + 0.2e1 * t203; m(7) * t3 + t38; m(7) * t18 + t76 * mrSges(7,2); m(7) * t147; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
