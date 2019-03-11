% Calculate time derivative of joint inertia matrix for
% S6RRRRRP6
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
% Datum: 2019-03-10 01:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP6_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP6_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP6_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP6_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP6_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP6_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:26:57
% EndTime: 2019-03-10 01:27:13
% DurationCPUTime: 7.72s
% Computational Cost: add. (10653->580), mult. (25242->829), div. (0->0), fcn. (22475->8), ass. (0->228)
t291 = qJD(3) + qJD(4);
t208 = sin(qJ(4));
t212 = cos(qJ(4));
t209 = sin(qJ(3));
t283 = -pkin(9) - pkin(8);
t249 = t283 * t209;
t213 = cos(qJ(3));
t293 = t283 * t213;
t140 = t208 * t293 + t212 * t249;
t297 = -mrSges(6,1) - mrSges(7,1);
t296 = Ifges(7,4) + Ifges(6,5);
t299 = -Ifges(7,2) - Ifges(6,3);
t295 = Ifges(7,6) - Ifges(6,6);
t210 = sin(qJ(2));
t258 = qJD(3) * t213;
t214 = cos(qJ(2));
t261 = qJD(2) * t214;
t224 = t209 * t261 + t210 * t258;
t298 = -t208 * t249 + t212 * t293;
t169 = t208 * t213 + t209 * t212;
t135 = t291 * t169;
t221 = -t169 * pkin(10) + t140;
t99 = t291 * t140;
t294 = -t135 * pkin(10) + qJD(5) * t221 + t99;
t229 = t208 * t209 - t212 * t213;
t155 = t229 * t210;
t180 = -pkin(2) * t214 - t210 * pkin(8) - pkin(1);
t167 = t213 * t180;
t267 = t210 * t213;
t280 = pkin(7) * t209;
t129 = -pkin(9) * t267 + t167 + (-pkin(3) - t280) * t214;
t265 = t213 * t214;
t193 = pkin(7) * t265;
t149 = t209 * t180 + t193;
t268 = t209 * t210;
t139 = -pkin(9) * t268 + t149;
t87 = t208 * t129 + t212 * t139;
t242 = t213 * t261;
t262 = qJD(2) * t210;
t292 = -Ifges(4,5) * t242 - Ifges(4,3) * t262;
t96 = -t135 * t210 - t229 * t261;
t97 = t155 * t291 - t169 * t261;
t290 = -Ifges(5,5) * t96 - Ifges(5,6) * t97 - Ifges(5,3) * t262;
t207 = sin(qJ(5));
t177 = (pkin(2) * t210 - pkin(8) * t214) * qJD(2);
t263 = t213 * t177 + t262 * t280;
t83 = (pkin(3) * t210 - pkin(9) * t265) * qJD(2) + (-t193 + (pkin(9) * t210 - t180) * t209) * qJD(3) + t263;
t260 = qJD(3) * t209;
t105 = t209 * t177 + t180 * t258 + (-t213 * t262 - t214 * t260) * pkin(7);
t93 = -pkin(9) * t224 + t105;
t26 = -qJD(4) * t87 - t208 * t93 + t212 * t83;
t21 = pkin(4) * t262 - pkin(10) * t96 + t26;
t211 = cos(qJ(5));
t256 = qJD(4) * t212;
t257 = qJD(4) * t208;
t25 = t129 * t256 - t139 * t257 + t208 * t83 + t212 * t93;
t23 = pkin(10) * t97 + t25;
t86 = t212 * t129 - t208 * t139;
t66 = -pkin(4) * t214 + t155 * pkin(10) + t86;
t154 = t169 * t210;
t73 = -pkin(10) * t154 + t87;
t278 = t207 * t66 + t211 * t73;
t6 = -qJD(5) * t278 - t207 * t23 + t21 * t211;
t289 = 2 * m(4);
t288 = 2 * m(5);
t287 = 2 * m(6);
t286 = 2 * m(7);
t285 = -0.2e1 * pkin(1);
t284 = 0.2e1 * pkin(7);
t282 = -t209 / 0.2e1;
t281 = m(7) * t207;
t277 = Ifges(4,4) * t209;
t276 = Ifges(4,4) * t213;
t275 = Ifges(4,6) * t209;
t274 = Ifges(4,6) * t214;
t196 = pkin(3) * t212 + pkin(4);
t253 = pkin(3) * t207 * t208;
t254 = qJD(5) * t211;
t121 = t196 * t254 + t211 * pkin(3) * t256 + (-qJD(4) - qJD(5)) * t253;
t273 = t121 * mrSges(6,2);
t255 = qJD(5) * t207;
t269 = t208 * t211;
t122 = t196 * t255 + (t208 * t254 + (t207 * t212 + t269) * qJD(4)) * pkin(3);
t272 = t122 * mrSges(7,1);
t115 = -pkin(10) * t229 - t298;
t71 = t207 * t115 - t211 * t221;
t271 = t122 * t71;
t128 = t169 * t211 - t207 * t229;
t270 = t122 * t128;
t182 = Ifges(4,1) * t209 + t276;
t266 = t213 * t182;
t110 = -t154 * t207 - t155 * t211;
t103 = -mrSges(6,1) * t214 - t110 * mrSges(6,3);
t104 = mrSges(7,1) * t214 + t110 * mrSges(7,2);
t264 = -t103 + t104;
t158 = pkin(3) * t269 + t207 * t196;
t178 = pkin(3) * t268 + t210 * pkin(7);
t259 = qJD(3) * t210;
t252 = pkin(4) * t255;
t251 = pkin(4) * t254;
t204 = pkin(7) * t261;
t203 = pkin(3) * t260;
t247 = t71 * t255;
t146 = t224 * pkin(3) + t204;
t197 = -pkin(3) * t213 - pkin(2);
t100 = t291 * t298;
t134 = t291 * t229;
t215 = t134 * pkin(10) + t100;
t19 = -t115 * t255 + t207 * t215 + t211 * t294;
t20 = t115 * t254 + t207 * t294 - t211 * t215;
t72 = t211 * t115 + t207 * t221;
t246 = t72 * t19 + t20 * t71;
t245 = t209 * t259;
t241 = (2 * Ifges(3,4)) + t275;
t117 = pkin(4) * t135 + t203;
t132 = pkin(4) * t154 + t178;
t231 = -t211 * t154 + t155 * t207;
t46 = qJD(5) * t231 + t207 * t97 + t211 * t96;
t39 = -mrSges(7,1) * t262 + t46 * mrSges(7,2);
t238 = -mrSges(4,1) * t213 + mrSges(4,2) * t209;
t237 = mrSges(4,1) * t209 + mrSges(4,2) * t213;
t236 = Ifges(4,1) * t213 - t277;
t235 = -Ifges(4,2) * t209 + t276;
t181 = Ifges(4,2) * t213 + t277;
t234 = Ifges(4,5) * t209 + Ifges(4,6) * t213;
t35 = -t207 * t73 + t211 * t66;
t230 = -t169 * t207 - t211 * t229;
t74 = -pkin(4) * t97 + t146;
t150 = pkin(4) * t229 + t197;
t47 = qJD(5) * t110 + t207 * t96 - t211 * t97;
t228 = t299 * t262 - t295 * t47 - t296 * t46;
t157 = t196 * t211 - t253;
t5 = t207 * t21 + t211 * t23 + t66 * t254 - t255 * t73;
t119 = qJD(6) + t121;
t116 = t119 * mrSges(7,3);
t118 = t122 * mrSges(6,1);
t227 = t116 - t118 - t272 - t273;
t226 = (-mrSges(5,1) * t208 - mrSges(5,2) * t212) * qJD(4) * pkin(3);
t225 = t242 - t245;
t57 = qJD(5) * t128 - t134 * t207 + t211 * t135;
t52 = Ifges(7,6) * t57;
t53 = Ifges(6,6) * t57;
t56 = qJD(5) * t230 - t134 * t211 - t135 * t207;
t54 = Ifges(6,5) * t56;
t55 = Ifges(7,4) * t56;
t223 = t52 - t53 + t54 + t55 + t297 * t20 + (-mrSges(6,2) + mrSges(7,3)) * t19;
t2 = qJ(6) * t262 - qJD(6) * t214 + t5;
t3 = -pkin(5) * t262 - t6;
t222 = t6 * mrSges(6,1) - t3 * mrSges(7,1) - t5 * mrSges(6,2) + t2 * mrSges(7,3) - t228;
t190 = qJD(6) + t251;
t183 = t190 * mrSges(7,3);
t218 = -mrSges(6,2) * t251 + t252 * t297 + t183;
t130 = Ifges(5,6) * t135;
t131 = Ifges(5,5) * t134;
t217 = t100 * mrSges(5,1) - t99 * mrSges(5,2) - t130 - t131 + t223;
t216 = t26 * mrSges(5,1) - t25 * mrSges(5,2) + t222 - t290;
t206 = qJD(6) * mrSges(7,3);
t202 = Ifges(4,5) * t258;
t195 = -pkin(4) * t211 - pkin(5);
t194 = pkin(4) * t207 + qJ(6);
t176 = -mrSges(4,1) * t214 - mrSges(4,3) * t267;
t175 = mrSges(4,2) * t214 - mrSges(4,3) * t268;
t174 = t236 * qJD(3);
t173 = t235 * qJD(3);
t172 = t237 * qJD(3);
t156 = -pkin(5) - t157;
t153 = qJ(6) + t158;
t152 = -Ifges(4,5) * t214 + t210 * t236;
t151 = t210 * t235 - t274;
t148 = -t214 * t280 + t167;
t145 = -mrSges(4,2) * t262 - mrSges(4,3) * t224;
t144 = mrSges(4,1) * t262 - mrSges(4,3) * t225;
t143 = -mrSges(5,1) * t214 + t155 * mrSges(5,3);
t142 = mrSges(5,2) * t214 - t154 * mrSges(5,3);
t138 = Ifges(5,1) * t169 - Ifges(5,4) * t229;
t137 = Ifges(5,4) * t169 - Ifges(5,2) * t229;
t136 = mrSges(5,1) * t229 + mrSges(5,2) * t169;
t125 = mrSges(4,1) * t224 + mrSges(4,2) * t225;
t114 = mrSges(5,1) * t154 - mrSges(5,2) * t155;
t113 = -t182 * t259 + (Ifges(4,5) * t210 + t214 * t236) * qJD(2);
t112 = -t181 * t259 + (Ifges(4,6) * t210 + t214 * t235) * qJD(2);
t108 = -Ifges(5,1) * t155 - Ifges(5,4) * t154 - Ifges(5,5) * t214;
t107 = -Ifges(5,4) * t155 - Ifges(5,2) * t154 - Ifges(5,6) * t214;
t106 = -qJD(3) * t149 + t263;
t102 = mrSges(6,2) * t214 + mrSges(6,3) * t231;
t101 = mrSges(7,2) * t231 - mrSges(7,3) * t214;
t91 = -Ifges(5,1) * t134 - Ifges(5,4) * t135;
t90 = -Ifges(5,4) * t134 - Ifges(5,2) * t135;
t89 = mrSges(5,1) * t135 - mrSges(5,2) * t134;
t85 = -mrSges(5,2) * t262 + mrSges(5,3) * t97;
t84 = mrSges(5,1) * t262 - mrSges(5,3) * t96;
t82 = Ifges(6,1) * t128 + Ifges(6,4) * t230;
t81 = Ifges(7,1) * t128 - Ifges(7,5) * t230;
t80 = Ifges(6,4) * t128 + Ifges(6,2) * t230;
t79 = Ifges(7,5) * t128 - Ifges(7,3) * t230;
t78 = -mrSges(6,1) * t230 + mrSges(6,2) * t128;
t77 = -mrSges(7,1) * t230 - mrSges(7,3) * t128;
t70 = -pkin(5) * t230 - qJ(6) * t128 + t150;
t68 = -mrSges(6,1) * t231 + mrSges(6,2) * t110;
t67 = -mrSges(7,1) * t231 - mrSges(7,3) * t110;
t62 = Ifges(6,1) * t110 + Ifges(6,4) * t231 - Ifges(6,5) * t214;
t61 = Ifges(7,1) * t110 - Ifges(7,4) * t214 - Ifges(7,5) * t231;
t60 = Ifges(6,4) * t110 + Ifges(6,2) * t231 - Ifges(6,6) * t214;
t59 = Ifges(7,5) * t110 - Ifges(7,6) * t214 - Ifges(7,3) * t231;
t51 = -pkin(5) * t231 - qJ(6) * t110 + t132;
t50 = -mrSges(5,1) * t97 + mrSges(5,2) * t96;
t49 = Ifges(5,1) * t96 + Ifges(5,4) * t97 + Ifges(5,5) * t262;
t48 = Ifges(5,4) * t96 + Ifges(5,2) * t97 + Ifges(5,6) * t262;
t40 = -mrSges(6,2) * t262 - mrSges(6,3) * t47;
t38 = mrSges(6,1) * t262 - mrSges(6,3) * t46;
t37 = -mrSges(7,2) * t47 + mrSges(7,3) * t262;
t34 = pkin(5) * t214 - t35;
t33 = -qJ(6) * t214 + t278;
t32 = Ifges(6,1) * t56 - Ifges(6,4) * t57;
t31 = Ifges(7,1) * t56 + Ifges(7,5) * t57;
t30 = Ifges(6,4) * t56 - Ifges(6,2) * t57;
t29 = Ifges(7,5) * t56 + Ifges(7,3) * t57;
t28 = mrSges(6,1) * t57 + mrSges(6,2) * t56;
t27 = mrSges(7,1) * t57 - mrSges(7,3) * t56;
t15 = pkin(5) * t57 - qJ(6) * t56 - qJD(6) * t128 + t117;
t14 = mrSges(6,1) * t47 + mrSges(6,2) * t46;
t13 = mrSges(7,1) * t47 - mrSges(7,3) * t46;
t12 = Ifges(6,1) * t46 - Ifges(6,4) * t47 + Ifges(6,5) * t262;
t11 = Ifges(7,1) * t46 + Ifges(7,4) * t262 + Ifges(7,5) * t47;
t10 = Ifges(6,4) * t46 - Ifges(6,2) * t47 + Ifges(6,6) * t262;
t9 = Ifges(7,5) * t46 + Ifges(7,6) * t262 + Ifges(7,3) * t47;
t7 = pkin(5) * t47 - qJ(6) * t46 - qJD(6) * t110 + t74;
t1 = [(t132 * t74 + t278 * t5 + t35 * t6) * t287 + 0.2e1 * t278 * t40 + (t125 * t284 - t209 * t112 + t213 * t113 + (-t213 * t151 - t209 * t152 + t214 * t234) * qJD(3) + (-Ifges(5,5) * t155 - Ifges(5,6) * t154 + mrSges(3,1) * t285 + (Ifges(4,5) * t213 - t241) * t210 + t296 * t110 - t295 * t231 + (pkin(7) ^ 2 * t289 + t237 * t284 + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3) - Ifges(5,3) + t299) * t214) * qJD(2)) * t210 + ((mrSges(3,2) * t285 - t209 * t151 + t213 * t152 + t214 * t241) * qJD(2) + t228 + t290 + t292) * t214 + (t59 - t60) * t47 + (t2 * t33 + t3 * t34 + t51 * t7) * t286 + (t146 * t178 + t25 * t87 + t26 * t86) * t288 + (t149 * t105 + t148 * t106) * t289 + 0.2e1 * t105 * t175 + 0.2e1 * t106 * t176 + 0.2e1 * t178 * t50 + (t61 + t62) * t46 - t154 * t48 - t155 * t49 + 0.2e1 * t25 * t142 + 0.2e1 * t26 * t143 + 0.2e1 * t146 * t114 + 0.2e1 * t148 * t144 + 0.2e1 * t149 * t145 + 0.2e1 * t132 * t14 + 0.2e1 * t2 * t101 + 0.2e1 * t5 * t102 + 0.2e1 * t6 * t103 + 0.2e1 * t3 * t104 + t97 * t107 + t96 * t108 + 0.2e1 * t86 * t84 + 0.2e1 * t87 * t85 + 0.2e1 * t74 * t68 + 0.2e1 * t7 * t67 + 0.2e1 * t51 * t13 + 0.2e1 * t33 * t37 + 0.2e1 * t35 * t38 + 0.2e1 * t34 * t39 - (t9 - t10) * t231 + (t11 + t12) * t110; -t298 * t85 + m(5) * (t100 * t86 + t140 * t26 + t146 * t197 + t178 * t203 - t25 * t298 + t87 * t99) + (-pkin(8) * t144 - t106 * mrSges(4,3) + t113 / 0.2e1) * t209 + m(7) * (t15 * t51 + t19 * t33 + t2 * t72 + t20 * t34 + t3 * t71 + t7 * t70) + m(6) * (t117 * t132 + t150 * t74 + t19 * t278 - t20 * t35 + t5 * t72 - t6 * t71) + (pkin(8) * t145 + t105 * mrSges(4,3) + t112 / 0.2e1) * t213 + ((t152 / 0.2e1 - pkin(8) * t176 - t148 * mrSges(4,3)) * t213 + (-t151 / 0.2e1 + t274 / 0.2e1 - pkin(8) * t175 - t149 * mrSges(4,3) + pkin(3) * t114) * t209) * qJD(3) + t264 * t20 + (t37 + t40) * t72 + (-t38 + t39) * t71 + (t101 + t102) * t19 + (t79 / 0.2e1 - t80 / 0.2e1) * t47 + (t81 / 0.2e1 + t82 / 0.2e1) * t46 + (t32 / 0.2e1 + t31 / 0.2e1) * t110 + t197 * t50 + t169 * t49 / 0.2e1 + t178 * t89 - t154 * t90 / 0.2e1 - t155 * t91 / 0.2e1 + t99 * t142 + t100 * t143 + t146 * t136 + t150 * t14 + t132 * t28 - t134 * t108 / 0.2e1 - t135 * t107 / 0.2e1 + t97 * t137 / 0.2e1 + t96 * t138 / 0.2e1 + t140 * t84 + t117 * t68 - pkin(2) * t125 + t7 * t77 + t74 * t78 + t15 * t67 + t70 * t13 + t51 * t27 + (-pkin(2) * t204 + (t105 * t213 - t106 * t209 + (-t148 * t213 - t149 * t209) * qJD(3)) * pkin(8)) * m(4) + (t131 / 0.2e1 + t130 / 0.2e1 - t54 / 0.2e1 + t53 / 0.2e1 - t55 / 0.2e1 - t52 / 0.2e1 - t202 / 0.2e1 + (Ifges(3,5) + t266 / 0.2e1 + t181 * t282 + (-mrSges(3,1) + t238) * pkin(7)) * qJD(2)) * t214 - t229 * t48 / 0.2e1 + (t213 * t174 / 0.2e1 + t173 * t282 - Ifges(3,6) * qJD(2) + (-t213 * t181 / 0.2e1 + t182 * t282) * qJD(3) + (qJD(2) * mrSges(3,2) + t172) * pkin(7) + (Ifges(5,5) * t169 - Ifges(5,6) * t229 + t296 * t128 - t230 * t295 + t234) * qJD(2) / 0.2e1) * t210 + (t134 * t86 - t135 * t87 - t169 * t26 - t229 * t25) * mrSges(5,3) - (t29 / 0.2e1 - t30 / 0.2e1) * t231 + (-t128 * t6 + t230 * t5 - t278 * t57 - t35 * t56) * mrSges(6,3) + (t128 * t3 + t2 * t230 - t33 * t57 + t34 * t56) * mrSges(7,2) - (t9 / 0.2e1 - t10 / 0.2e1) * t230 + (t11 / 0.2e1 + t12 / 0.2e1) * t128 + (t61 / 0.2e1 + t62 / 0.2e1) * t56 + (t59 / 0.2e1 - t60 / 0.2e1) * t57; -0.2e1 * pkin(2) * t172 + 0.2e1 * t117 * t78 - t134 * t138 - t135 * t137 + 0.2e1 * t15 * t77 + 0.2e1 * t150 * t28 - t229 * t90 + t169 * t91 + t213 * t173 + t209 * t174 + 0.2e1 * t197 * t89 + 0.2e1 * t70 * t27 + (t79 - t80) * t57 + (t81 + t82) * t56 + (t32 + t31) * t128 - (t29 - t30) * t230 + (t266 + (0.2e1 * pkin(3) * t136 - t181) * t209) * qJD(3) + (t100 * t140 + t197 * t203 - t298 * t99) * t288 + (t117 * t150 + t246) * t287 + (t15 * t70 + t246) * t286 + 0.2e1 * (-t100 * t169 + t134 * t140 + t135 * t298 - t229 * t99) * mrSges(5,3) + 0.2e1 * (mrSges(6,3) + mrSges(7,2)) * (t128 * t20 + t19 * t230 + t56 * t71 - t57 * t72); t156 * t39 + t157 * t38 + t158 * t40 + t153 * t37 + t119 * t101 + t121 * t102 - t105 * mrSges(4,2) + t106 * mrSges(4,1) - Ifges(4,5) * t245 + t264 * t122 + (m(5) * (t208 * t25 + t212 * t26 + t256 * t87 - t257 * t86) + t212 * t84 + t208 * t85 + t142 * t256 - t143 * t257) * pkin(3) + t216 - t224 * Ifges(4,6) + m(7) * (t119 * t33 + t122 * t34 + t153 * t2 + t156 * t3) + m(6) * (t121 * t278 - t122 * t35 + t157 * t6 + t158 * t5) - t292; t202 + (m(5) * (t100 * t212 + t208 * t99 + (-t140 * t208 - t212 * t298) * qJD(4)) + (t212 * t134 - t208 * t135 + (t169 * t208 - t212 * t229) * qJD(4)) * mrSges(5,3)) * pkin(3) + t217 + m(7) * (t119 * t72 + t153 * t19 + t156 * t20 + t271) + m(6) * (t121 * t72 - t157 * t20 + t158 * t19 + t271) + (t121 * t230 - t157 * t56 - t158 * t57 + t270) * mrSges(6,3) + (t119 * t230 - t153 * t57 + t156 * t56 + t270) * mrSges(7,2) + (pkin(8) * t238 - t275) * qJD(3); -0.2e1 * t272 - 0.2e1 * t273 + 0.2e1 * t116 - 0.2e1 * t118 + 0.2e1 * t226 + (t121 * t158 - t122 * t157) * t287 + (t119 * t153 + t122 * t156) * t286; t194 * t37 + t195 * t39 + t190 * t101 + m(7) * (t190 * t33 + t194 * t2 + t195 * t3) + (m(6) * (t207 * t5 + t211 * t6) + t211 * t38 + t207 * t40 + (t211 * t102 + t264 * t207 + m(6) * (-t207 * t35 + t211 * t278) + t34 * t281) * qJD(5)) * pkin(4) + t216; m(7) * (t19 * t194 + t190 * t72 + t195 * t20) + (t128 * mrSges(7,2) * t255 + m(7) * t247 + m(6) * (t19 * t207 - t20 * t211 + t254 * t72 + t247) + (-t207 * t57 - t211 * t56 + (t128 * t207 + t211 * t230) * qJD(5)) * mrSges(6,3)) * pkin(4) + t217 + (t190 * t230 - t194 * t57 + t195 * t56) * mrSges(7,2); t183 + m(7) * (t119 * t194 + t122 * t195 + t153 * t190) + t226 + (m(6) * (t121 * t207 - t122 * t211) + (-t211 * mrSges(6,2) + t297 * t207 + m(6) * (-t157 * t207 + t158 * t211) + t156 * t281) * qJD(5)) * pkin(4) + t227; 0.2e1 * m(7) * (t190 * t194 + t195 * t252) + 0.2e1 * t218; t222 + qJD(6) * t101 + qJ(6) * t37 - pkin(5) * t39 + m(7) * (-pkin(5) * t3 + qJ(6) * t2 + qJD(6) * t33); m(7) * (-pkin(5) * t20 + qJ(6) * t19 + qJD(6) * t72) + (-pkin(5) * t56 - qJ(6) * t57 + qJD(6) * t230) * mrSges(7,2) + t223; t206 + m(7) * (-pkin(5) * t122 + qJ(6) * t119 + qJD(6) * t153) + t227; t206 + m(7) * (-pkin(5) * t252 + qJ(6) * t190 + qJD(6) * t194) + t218; 0.2e1 * m(7) * qJ(6) * qJD(6) + 0.2e1 * t206; m(7) * t3 + t39; m(7) * t20 + t56 * mrSges(7,2); m(7) * t122; m(7) * t252; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
