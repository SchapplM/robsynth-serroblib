% Calculate time derivative of joint inertia matrix for
% S6RRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 22:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR8_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR8_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR8_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR8_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:38:54
% EndTime: 2019-03-09 22:39:07
% DurationCPUTime: 6.02s
% Computational Cost: add. (7784->590), mult. (18151->838), div. (0->0), fcn. (15664->8), ass. (0->227)
t285 = -mrSges(5,1) - mrSges(6,1);
t282 = Ifges(5,5) + Ifges(6,4);
t284 = -Ifges(6,2) - Ifges(5,3);
t283 = Ifges(5,6) - Ifges(6,6);
t207 = sin(qJ(2));
t206 = sin(qJ(3));
t211 = cos(qJ(2));
t247 = qJD(2) * t211;
t234 = t206 * t247;
t210 = cos(qJ(3));
t244 = qJD(3) * t210;
t219 = t207 * t244 + t234;
t204 = sin(qJ(6));
t208 = cos(qJ(6));
t205 = sin(qJ(4));
t209 = cos(qJ(4));
t159 = t205 * t206 - t209 * t210;
t149 = t159 * t207;
t174 = -pkin(2) * t211 - pkin(8) * t207 - pkin(1);
t158 = t210 * t174;
t254 = t207 * t210;
t267 = pkin(7) * t206;
t115 = -pkin(9) * t254 + t158 + (-pkin(3) - t267) * t211;
t252 = t210 * t211;
t187 = pkin(7) * t252;
t141 = t206 * t174 + t187;
t255 = t206 * t207;
t130 = -pkin(9) * t255 + t141;
t65 = t115 * t209 - t205 * t130;
t55 = t211 * pkin(4) - t65;
t42 = pkin(5) * t211 + pkin(10) * t149 + t55;
t160 = t205 * t210 + t206 * t209;
t148 = t160 * t207;
t66 = t205 * t115 + t209 * t130;
t54 = -qJ(5) * t211 + t66;
t45 = pkin(10) * t148 + t54;
t22 = -t204 * t45 + t208 * t42;
t242 = qJD(4) * t209;
t243 = qJD(4) * t205;
t170 = (pkin(2) * t207 - pkin(8) * t211) * qJD(2);
t248 = qJD(2) * t207;
t250 = t210 * t170 + t248 * t267;
t60 = (pkin(3) * t207 - pkin(9) * t252) * qJD(2) + (-t187 + (pkin(9) * t207 - t174) * t206) * qJD(3) + t250;
t246 = qJD(3) * t206;
t91 = t206 * t170 + t174 * t244 + (-t210 * t248 - t211 * t246) * pkin(7);
t75 = -pkin(9) * t219 + t91;
t18 = -t115 * t243 - t130 * t242 - t205 * t75 + t209 * t60;
t212 = -pkin(4) - pkin(5);
t278 = qJD(3) + qJD(4);
t81 = -t148 * t278 - t159 * t247;
t7 = -pkin(10) * t81 + t212 * t248 - t18;
t17 = t115 * t242 - t130 * t243 + t205 * t60 + t209 * t75;
t14 = qJ(5) * t248 - qJD(5) * t211 + t17;
t233 = t210 * t247;
t245 = qJD(3) * t207;
t237 = t206 * t245;
t220 = t233 - t237;
t82 = -t243 * t255 + (t254 * t278 + t234) * t209 + t220 * t205;
t8 = pkin(10) * t82 + t14;
t2 = qJD(6) * t22 + t204 * t7 + t208 * t8;
t23 = t204 * t42 + t208 * t45;
t3 = -qJD(6) * t23 - t204 * t8 + t208 * t7;
t280 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t270 = -pkin(9) - pkin(8);
t179 = t270 * t206;
t180 = t270 * t210;
t132 = t205 * t179 - t209 * t180;
t279 = -Ifges(4,5) * t233 - Ifges(4,3) * t248;
t277 = (-mrSges(5,2) * t209 + t205 * t285) * pkin(3) * qJD(4);
t276 = 2 * m(4);
t275 = 2 * m(5);
t274 = 2 * m(6);
t273 = 2 * m(7);
t272 = -0.2e1 * pkin(1);
t271 = 0.2e1 * pkin(7);
t268 = -t206 / 0.2e1;
t192 = -pkin(3) * t209 - pkin(4);
t188 = -pkin(5) + t192;
t189 = pkin(3) * t205 + qJ(5);
t142 = t188 * t208 - t189 * t204;
t185 = pkin(3) * t242 + qJD(5);
t239 = pkin(3) * t243;
t99 = qJD(6) * t142 + t185 * t208 + t204 * t239;
t266 = t99 * mrSges(7,2);
t263 = Ifges(4,4) * t206;
t262 = Ifges(4,4) * t210;
t261 = Ifges(4,6) * t206;
t260 = Ifges(4,6) * t211;
t143 = t188 * t204 + t189 * t208;
t100 = -qJD(6) * t143 - t185 * t204 + t208 * t239;
t258 = t100 * mrSges(7,1);
t172 = -t204 * qJ(5) + t208 * t212;
t146 = t208 * qJD(5) + qJD(6) * t172;
t257 = t146 * mrSges(7,2);
t173 = t208 * qJ(5) + t204 * t212;
t147 = -t204 * qJD(5) - qJD(6) * t173;
t256 = t147 * mrSges(7,1);
t178 = Ifges(4,1) * t206 + t262;
t253 = t210 * t178;
t135 = -mrSges(5,1) * t211 + mrSges(5,3) * t149;
t136 = mrSges(6,1) * t211 - mrSges(6,2) * t149;
t251 = -t135 + t136;
t249 = (-mrSges(7,1) * t204 - mrSges(7,2) * t208) * qJD(6);
t171 = pkin(3) * t255 + t207 * pkin(7);
t97 = t148 * t208 + t149 * t204;
t28 = qJD(6) * t97 + t204 * t82 + t208 * t81;
t98 = t148 * t204 - t149 * t208;
t29 = -qJD(6) * t98 - t204 * t81 + t208 * t82;
t241 = Ifges(7,5) * t28 + Ifges(7,6) * t29 - Ifges(7,3) * t248;
t240 = pkin(3) * t246;
t201 = pkin(7) * t247;
t139 = pkin(3) * t219 + t201;
t193 = -pkin(3) * t210 - pkin(2);
t238 = qJD(3) * t270;
t131 = -t209 * t179 - t180 * t205;
t235 = t131 * t243;
t231 = (2 * Ifges(3,4)) + t261;
t169 = t206 * t238;
t229 = t210 * t238;
t85 = t209 * t169 + t179 * t242 + t180 * t243 + t205 * t229;
t86 = qJD(4) * t132 + t169 * t205 - t209 * t229;
t230 = t131 * t86 + t132 * t85;
t63 = -mrSges(6,1) * t248 + t81 * mrSges(6,2);
t228 = -qJ(5) * t149 - t171;
t123 = t278 * t160;
t46 = pkin(10) * t123 + t85;
t122 = t278 * t159;
t47 = pkin(10) * t122 + t86;
t105 = -pkin(10) * t160 + t131;
t106 = pkin(10) * t159 + t132;
t49 = t105 * t208 - t106 * t204;
t10 = qJD(6) * t49 + t204 * t47 + t208 * t46;
t50 = t105 * t204 + t106 * t208;
t11 = -qJD(6) * t50 - t204 * t46 + t208 * t47;
t113 = t159 * t204 + t160 * t208;
t41 = -qJD(6) * t113 + t122 * t204 + t123 * t208;
t38 = Ifges(7,6) * t41;
t112 = t159 * t208 - t160 * t204;
t40 = qJD(6) * t112 - t122 * t208 + t123 * t204;
t39 = Ifges(7,5) * t40;
t227 = t11 * mrSges(7,1) - t10 * mrSges(7,2) + t38 + t39;
t226 = -mrSges(4,1) * t210 + mrSges(4,2) * t206;
t225 = mrSges(4,1) * t206 + mrSges(4,2) * t210;
t224 = Ifges(4,1) * t210 - t263;
t223 = -Ifges(4,2) * t206 + t262;
t177 = Ifges(4,2) * t210 + t263;
t222 = Ifges(4,5) * t206 + Ifges(4,6) * t210;
t221 = qJ(5) * t160 - t193;
t218 = t248 * t284 - t282 * t81 + t283 * t82 + t241;
t217 = qJ(5) * t81 - qJD(5) * t149 - t139;
t216 = -qJ(5) * t122 + qJD(5) * t160 - t240;
t116 = Ifges(6,6) * t123;
t117 = Ifges(5,6) * t123;
t118 = Ifges(5,5) * t122;
t119 = Ifges(6,4) * t122;
t214 = t116 - t117 - t118 - t119 - t227 + t285 * t86 + (-mrSges(5,2) + mrSges(6,3)) * t85;
t15 = -pkin(4) * t248 - t18;
t213 = t18 * mrSges(5,1) - t15 * mrSges(6,1) - t17 * mrSges(5,2) + t14 * mrSges(6,3) - t218 - t280;
t200 = Ifges(4,5) * t244;
t168 = -mrSges(4,1) * t211 - mrSges(4,3) * t254;
t167 = mrSges(4,2) * t211 - mrSges(4,3) * t255;
t166 = t224 * qJD(3);
t165 = t223 * qJD(3);
t164 = t225 * qJD(3);
t145 = -Ifges(4,5) * t211 + t207 * t224;
t144 = t207 * t223 - t260;
t140 = -t211 * t267 + t158;
t138 = -mrSges(4,2) * t248 - mrSges(4,3) * t219;
t137 = mrSges(4,1) * t248 - mrSges(4,3) * t220;
t134 = mrSges(5,2) * t211 - mrSges(5,3) * t148;
t133 = -mrSges(6,2) * t148 - mrSges(6,3) * t211;
t129 = Ifges(5,1) * t160 - Ifges(5,4) * t159;
t128 = Ifges(6,1) * t160 + Ifges(6,5) * t159;
t127 = Ifges(5,4) * t160 - Ifges(5,2) * t159;
t126 = Ifges(6,5) * t160 + Ifges(6,3) * t159;
t125 = mrSges(5,1) * t159 + mrSges(5,2) * t160;
t124 = mrSges(6,1) * t159 - mrSges(6,3) * t160;
t111 = mrSges(4,1) * t219 + mrSges(4,2) * t220;
t110 = pkin(4) * t159 - t221;
t104 = mrSges(5,1) * t148 - mrSges(5,2) * t149;
t103 = mrSges(6,1) * t148 + mrSges(6,3) * t149;
t102 = -t178 * t245 + (Ifges(4,5) * t207 + t211 * t224) * qJD(2);
t101 = -t177 * t245 + (Ifges(4,6) * t207 + t211 * t223) * qJD(2);
t96 = -Ifges(5,1) * t149 - Ifges(5,4) * t148 - Ifges(5,5) * t211;
t95 = -Ifges(6,1) * t149 - Ifges(6,4) * t211 + Ifges(6,5) * t148;
t94 = -Ifges(5,4) * t149 - Ifges(5,2) * t148 - Ifges(5,6) * t211;
t93 = -Ifges(6,5) * t149 - Ifges(6,6) * t211 + Ifges(6,3) * t148;
t92 = -qJD(3) * t141 + t250;
t90 = t159 * t212 + t221;
t89 = mrSges(7,1) * t211 - mrSges(7,3) * t98;
t88 = -mrSges(7,2) * t211 + mrSges(7,3) * t97;
t87 = pkin(4) * t148 - t228;
t74 = -Ifges(5,1) * t122 - Ifges(5,4) * t123;
t73 = -Ifges(6,1) * t122 + Ifges(6,5) * t123;
t72 = -Ifges(5,4) * t122 - Ifges(5,2) * t123;
t71 = -Ifges(6,5) * t122 + Ifges(6,3) * t123;
t70 = mrSges(5,1) * t123 - mrSges(5,2) * t122;
t69 = mrSges(6,1) * t123 + mrSges(6,3) * t122;
t64 = -mrSges(5,2) * t248 - mrSges(5,3) * t82;
t62 = mrSges(5,1) * t248 - mrSges(5,3) * t81;
t61 = -mrSges(6,2) * t82 + mrSges(6,3) * t248;
t59 = Ifges(7,1) * t113 + Ifges(7,4) * t112;
t58 = Ifges(7,4) * t113 + Ifges(7,2) * t112;
t57 = -mrSges(7,1) * t112 + mrSges(7,2) * t113;
t56 = t148 * t212 + t228;
t52 = pkin(4) * t123 - t216;
t48 = -mrSges(7,1) * t97 + mrSges(7,2) * t98;
t44 = Ifges(7,1) * t98 + Ifges(7,4) * t97 + Ifges(7,5) * t211;
t43 = Ifges(7,4) * t98 + Ifges(7,2) * t97 + Ifges(7,6) * t211;
t37 = t123 * t212 + t216;
t36 = mrSges(5,1) * t82 + mrSges(5,2) * t81;
t35 = mrSges(6,1) * t82 - mrSges(6,3) * t81;
t34 = Ifges(5,1) * t81 - Ifges(5,4) * t82 + Ifges(5,5) * t248;
t33 = Ifges(6,1) * t81 + Ifges(6,4) * t248 + Ifges(6,5) * t82;
t32 = Ifges(5,4) * t81 - Ifges(5,2) * t82 + Ifges(5,6) * t248;
t31 = Ifges(6,5) * t81 + Ifges(6,6) * t248 + Ifges(6,3) * t82;
t30 = pkin(4) * t82 - t217;
t25 = mrSges(7,2) * t248 + mrSges(7,3) * t29;
t24 = -mrSges(7,1) * t248 - mrSges(7,3) * t28;
t21 = Ifges(7,1) * t40 + Ifges(7,4) * t41;
t20 = Ifges(7,4) * t40 + Ifges(7,2) * t41;
t19 = -mrSges(7,1) * t41 + mrSges(7,2) * t40;
t12 = t212 * t82 + t217;
t6 = -mrSges(7,1) * t29 + mrSges(7,2) * t28;
t5 = Ifges(7,1) * t28 + Ifges(7,4) * t29 - Ifges(7,5) * t248;
t4 = Ifges(7,4) * t28 + Ifges(7,2) * t29 - Ifges(7,6) * t248;
t1 = [(t111 * t271 - t206 * t101 + t210 * t102 + (-t210 * t144 - t206 * t145 + t211 * t222) * qJD(3) + (-Ifges(7,5) * t98 - Ifges(7,6) * t97 + mrSges(3,1) * t272 + (Ifges(4,5) * t210 - t231) * t207 - t282 * t149 - t283 * t148 + (pkin(7) ^ 2 * t276 + t225 * t271 + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3) - Ifges(7,3) + t284) * t211) * qJD(2)) * t207 + ((mrSges(3,2) * t272 - t206 * t144 + t210 * t145 + t211 * t231) * qJD(2) + t218 + t279) * t211 - (t33 + t34) * t149 + (t31 - t32) * t148 + (t12 * t56 + t2 * t23 + t22 * t3) * t273 + (t14 * t54 + t15 * t55 + t30 * t87) * t274 + (t139 * t171 + t17 * t66 + t18 * t65) * t275 + (t140 * t92 + t141 * t91) * t276 + (t95 + t96) * t81 + (t93 - t94) * t82 + 0.2e1 * t22 * t24 + 0.2e1 * t23 * t25 + t29 * t43 + t28 * t44 + 0.2e1 * t12 * t48 + 0.2e1 * t56 * t6 + 0.2e1 * t54 * t61 + 0.2e1 * t55 * t63 + 0.2e1 * t65 * t62 + 0.2e1 * t66 * t64 + 0.2e1 * t87 * t35 + 0.2e1 * t2 * t88 + 0.2e1 * t3 * t89 + t97 * t4 + t98 * t5 + 0.2e1 * t30 * t103 + 0.2e1 * t14 * t133 + 0.2e1 * t17 * t134 + 0.2e1 * t18 * t135 + 0.2e1 * t15 * t136 + 0.2e1 * t139 * t104 + 0.2e1 * t140 * t137 + 0.2e1 * t141 * t138 + 0.2e1 * t91 * t167 + 0.2e1 * t92 * t168 + 0.2e1 * t171 * t36; (pkin(8) * t138 + t91 * mrSges(4,3) + t101 / 0.2e1) * t210 + (t210 * t166 / 0.2e1 + t165 * t268 - qJD(2) * (Ifges(7,5) * t113 + Ifges(7,6) * t112) / 0.2e1 - Ifges(3,6) * qJD(2) + (-t210 * t177 / 0.2e1 + t178 * t268) * qJD(3) + (mrSges(3,2) * qJD(2) + t164) * pkin(7) + (-t159 * t283 + t282 * t160 + t222) * qJD(2) / 0.2e1) * t207 + (t128 / 0.2e1 + t129 / 0.2e1) * t81 + (t122 * t65 - t123 * t66 - t159 * t17 - t160 * t18) * mrSges(5,3) + (-t122 * t55 - t123 * t54 - t14 * t159 + t15 * t160) * mrSges(6,2) - (t73 / 0.2e1 + t74 / 0.2e1) * t149 + (t71 / 0.2e1 - t72 / 0.2e1) * t148 + (t133 + t134) * t85 + (-pkin(2) * t201 + (-t92 * t206 + t91 * t210 + (-t140 * t210 - t141 * t206) * qJD(3)) * pkin(8)) * m(4) + (t61 + t64) * t132 + (t63 - t62) * t131 + (-pkin(8) * t137 - t92 * mrSges(4,3) + t102 / 0.2e1) * t206 + m(6) * (t110 * t30 + t131 * t15 + t132 * t14 + t52 * t87 + t54 * t85 + t55 * t86) + m(7) * (t10 * t23 + t11 * t22 + t12 * t90 + t2 * t50 + t3 * t49 + t37 * t56) + (t112 * t2 - t113 * t3 - t22 * t40 + t23 * t41) * mrSges(7,3) + (t126 / 0.2e1 - t127 / 0.2e1) * t82 + m(5) * (-t131 * t18 + t132 * t17 + t139 * t193 + t171 * t240 - t65 * t86 + t66 * t85) + t251 * t86 + ((-pkin(8) * t168 - t140 * mrSges(4,3) + t145 / 0.2e1) * t210 + (-pkin(8) * t167 - t141 * mrSges(4,3) + pkin(3) * t104 - t144 / 0.2e1 + t260 / 0.2e1) * t206) * qJD(3) - (t95 / 0.2e1 + t96 / 0.2e1) * t122 + (t93 / 0.2e1 - t94 / 0.2e1) * t123 + (-t200 / 0.2e1 + t39 / 0.2e1 + t38 / 0.2e1 + t118 / 0.2e1 + t117 / 0.2e1 + t119 / 0.2e1 - t116 / 0.2e1 + (Ifges(3,5) + t253 / 0.2e1 + t177 * t268 + (-mrSges(3,1) + t226) * pkin(7)) * qJD(2)) * t211 + (t31 / 0.2e1 - t32 / 0.2e1) * t159 + (t33 / 0.2e1 + t34 / 0.2e1) * t160 + t41 * t43 / 0.2e1 + t40 * t44 / 0.2e1 + t37 * t48 + t49 * t24 + t50 * t25 + t56 * t19 + t12 * t57 + t29 * t58 / 0.2e1 + t28 * t59 / 0.2e1 + t87 * t69 + t10 * t88 + t11 * t89 + t90 * t6 + t97 * t20 / 0.2e1 + t98 * t21 / 0.2e1 + t52 * t103 + t110 * t35 - pkin(2) * t111 + t112 * t4 / 0.2e1 + t113 * t5 / 0.2e1 + t30 * t124 + t139 * t125 + t171 * t70 + t193 * t36; -0.2e1 * pkin(2) * t164 + 0.2e1 * t110 * t69 + t112 * t20 + t113 * t21 + 0.2e1 * t52 * t124 + t210 * t165 + t206 * t166 + 0.2e1 * t90 * t19 + 0.2e1 * t193 * t70 + 0.2e1 * t37 * t57 + t40 * t59 + t41 * t58 + (t73 + t74) * t160 + (t71 - t72) * t159 + (t126 - t127) * t123 - (t128 + t129) * t122 + (t253 + (0.2e1 * pkin(3) * t125 - t177) * t206) * qJD(3) + (t193 * t240 + t230) * t275 + (t110 * t52 + t230) * t274 + (t10 * t50 + t11 * t49 + t37 * t90) * t273 + 0.2e1 * (t10 * t112 - t11 * t113 - t40 * t49 + t41 * t50) * mrSges(7,3) + 0.2e1 * (mrSges(5,3) + mrSges(6,2)) * (-t122 * t131 - t123 * t132 - t159 * t85 + t160 * t86); t213 + ((m(5) * t18 + t62 + (m(5) * t66 + t134) * qJD(4)) * t209 + (m(5) * t17 + t64 + (-m(5) * t65 + m(6) * t55 + t251) * qJD(4)) * t205) * pkin(3) - Ifges(4,5) * t237 + m(6) * (t14 * t189 + t15 * t192 + t185 * t54) + m(7) * (t100 * t22 + t142 * t3 + t143 * t2 + t23 * t99) - t219 * Ifges(4,6) - t91 * mrSges(4,2) + t92 * mrSges(4,1) + t99 * t88 + t100 * t89 + t142 * t24 + t143 * t25 + t185 * t133 + t189 * t61 + t192 * t63 - t279; t200 + t214 + (t160 * mrSges(6,2) * t243 + m(5) * (t132 * t242 + t205 * t85 - t209 * t86 + t235) + m(6) * t235 + (t209 * t122 - t205 * t123 + (-t159 * t209 + t160 * t205) * qJD(4)) * mrSges(5,3)) * pkin(3) + (pkin(8) * t226 - t261) * qJD(3) + m(6) * (t132 * t185 + t189 * t85 + t192 * t86) + (-t122 * t192 - t123 * t189 - t159 * t185) * mrSges(6,2) + m(7) * (t10 * t143 + t100 * t49 + t11 * t142 + t50 * t99) + (-t100 * t113 + t112 * t99 - t142 * t40 + t143 * t41) * mrSges(7,3); -0.2e1 * t258 + 0.2e1 * t266 + 0.2e1 * t185 * mrSges(6,3) + 0.2e1 * t277 + (t100 * t142 + t143 * t99) * t273 + (t185 * t189 + t192 * t239) * t274; t213 + m(7) * (t146 * t23 + t147 * t22 + t172 * t3 + t173 * t2) + m(6) * (-pkin(4) * t15 + qJ(5) * t14 + qJD(5) * t54) + qJ(5) * t61 - pkin(4) * t63 + qJD(5) * t133 + t146 * t88 + t147 * t89 + t172 * t24 + t173 * t25; t214 + m(7) * (t10 * t173 + t11 * t172 + t146 * t50 + t147 * t49) + m(6) * (-pkin(4) * t86 + qJ(5) * t85 + qJD(5) * t132) + (t112 * t146 - t113 * t147 - t172 * t40 + t173 * t41) * mrSges(7,3) + (pkin(4) * t122 - qJ(5) * t123 - qJD(5) * t159) * mrSges(6,2); (qJD(5) + t185) * mrSges(6,3) + (t146 + t99) * mrSges(7,2) + (-t147 - t100) * mrSges(7,1) + t277 + m(7) * (t100 * t172 + t142 * t147 + t143 * t146 + t173 * t99) + m(6) * (-pkin(4) * t239 + qJ(5) * t185 + qJD(5) * t189); (t146 * t173 + t147 * t172) * t273 + 0.2e1 * t257 - 0.2e1 * t256 + 0.2e1 * (m(6) * qJ(5) + mrSges(6,3)) * qJD(5); t204 * t25 + t208 * t24 + (-t204 * t89 + t208 * t88) * qJD(6) + m(7) * (t2 * t204 + t208 * t3 + (-t204 * t22 + t208 * t23) * qJD(6)) + m(6) * t15 + t63; -t122 * mrSges(6,2) + m(7) * (t10 * t204 + t11 * t208 + (-t204 * t49 + t208 * t50) * qJD(6)) + m(6) * t86 + (t204 * t41 - t208 * t40 + (t112 * t208 + t113 * t204) * qJD(6)) * mrSges(7,3); m(7) * (t100 * t208 + t204 * t99 + (-t142 * t204 + t143 * t208) * qJD(6)) + m(6) * t239 - t249; m(7) * (t146 * t204 + t147 * t208 + (-t172 * t204 + t173 * t208) * qJD(6)) - t249; 0; t241 + t280; t227; t258 - t266; t256 - t257; t249; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
