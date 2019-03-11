% Calculate time derivative of joint inertia matrix for
% S6RRRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-03-09 21:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:58:31
% EndTime: 2019-03-09 20:58:43
% DurationCPUTime: 5.82s
% Computational Cost: add. (7874->539), mult. (19423->771), div. (0->0), fcn. (16961->8), ass. (0->221)
t284 = -mrSges(6,1) - mrSges(7,1);
t283 = Ifges(7,4) + Ifges(6,5);
t282 = Ifges(7,6) - Ifges(6,6);
t207 = sin(qJ(3));
t208 = sin(qJ(2));
t210 = cos(qJ(3));
t245 = qJD(3) * t210;
t211 = cos(qJ(2));
t248 = qJD(2) * t211;
t218 = t207 * t248 + t208 * t245;
t286 = -Ifges(7,2) - Ifges(5,3) - Ifges(6,3);
t206 = sin(qJ(4));
t209 = cos(qJ(4));
t271 = -pkin(9) - pkin(8);
t239 = t271 * t207;
t281 = t271 * t210;
t285 = -t206 * t239 + t209 * t281;
t138 = t206 * t281 + t209 * t239;
t221 = t206 * t207 - t209 * t210;
t157 = t221 * t208;
t180 = -pkin(2) * t211 - t208 * pkin(8) - pkin(1);
t168 = t210 * t180;
t254 = t208 * t210;
t268 = pkin(7) * t207;
t127 = -pkin(9) * t254 + t168 + (-pkin(3) - t268) * t211;
t252 = t210 * t211;
t191 = pkin(7) * t252;
t147 = t207 * t180 + t191;
t255 = t207 * t208;
t137 = -pkin(9) * t255 + t147;
t90 = t206 * t127 + t209 * t137;
t232 = t210 * t248;
t249 = qJD(2) * t208;
t280 = -Ifges(4,5) * t232 - Ifges(4,3) * t249;
t279 = qJD(3) + qJD(4);
t278 = 2 * m(4);
t277 = 2 * m(5);
t276 = 2 * m(6);
t275 = 2 * m(7);
t274 = -0.2e1 * pkin(1);
t273 = 0.2e1 * pkin(7);
t272 = m(6) * pkin(4);
t270 = -t207 / 0.2e1;
t178 = (pkin(2) * t208 - pkin(8) * t211) * qJD(2);
t250 = t210 * t178 + t249 * t268;
t84 = (pkin(3) * t208 - pkin(9) * t252) * qJD(2) + (-t191 + (pkin(9) * t208 - t180) * t207) * qJD(3) + t250;
t247 = qJD(3) * t207;
t108 = t207 * t178 + t180 * t245 + (-t210 * t249 - t211 * t247) * pkin(7);
t96 = -pkin(9) * t218 + t108;
t26 = -qJD(4) * t90 - t206 * t96 + t209 * t84;
t170 = t206 * t210 + t207 * t209;
t133 = t279 * t170;
t99 = -t133 * t208 - t221 * t248;
t10 = pkin(4) * t249 - qJ(5) * t99 + qJD(5) * t157 + t26;
t100 = t157 * t279 - t170 * t248;
t156 = t170 * t208;
t243 = qJD(4) * t209;
t244 = qJD(4) * t206;
t25 = t127 * t243 - t137 * t244 + t206 * t84 + t209 * t96;
t12 = qJ(5) * t100 - qJD(5) * t156 + t25;
t205 = sin(pkin(10));
t257 = cos(pkin(10));
t6 = t205 * t10 + t257 * t12;
t269 = pkin(4) * t205;
t132 = t279 * t221;
t88 = -t132 * t257 - t205 * t133;
t267 = t88 * mrSges(7,2);
t89 = t209 * t127 - t206 * t137;
t63 = -pkin(4) * t211 + t157 * qJ(5) + t89;
t68 = -qJ(5) * t156 + t90;
t30 = t205 * t63 + t257 * t68;
t265 = Ifges(4,4) * t207;
t264 = Ifges(4,4) * t210;
t263 = Ifges(4,6) * t207;
t262 = Ifges(4,6) * t211;
t261 = pkin(3) * qJD(4);
t241 = pkin(3) * t243;
t242 = pkin(3) * t244;
t155 = -t205 * t242 + t257 * t241;
t149 = qJD(6) + t155;
t260 = t149 * mrSges(7,3);
t230 = t257 * t206;
t154 = (t205 * t209 + t230) * t261;
t117 = -qJ(5) * t221 - t285;
t214 = -t170 * qJ(5) + t138;
t66 = t117 * t205 - t214 * t257;
t259 = t154 * t66;
t258 = t155 * mrSges(6,2);
t126 = t170 * t257 - t205 * t221;
t256 = t126 * t154;
t182 = Ifges(4,1) * t207 + t264;
t253 = t210 * t182;
t113 = -t205 * t156 - t157 * t257;
t106 = -mrSges(6,1) * t211 - t113 * mrSges(6,3);
t107 = mrSges(7,1) * t211 + t113 * mrSges(7,2);
t251 = t107 - t106;
t194 = pkin(3) * t209 + pkin(4);
t159 = pkin(3) * t230 + t205 * t194;
t179 = pkin(3) * t255 + t208 * pkin(7);
t246 = qJD(3) * t208;
t202 = pkin(7) * t248;
t201 = pkin(3) * t247;
t144 = pkin(3) * t218 + t202;
t195 = -pkin(3) * t210 - pkin(2);
t216 = t285 * qJD(3);
t103 = qJD(4) * t285 + t216;
t220 = t132 * qJ(5) - t170 * qJD(5);
t227 = qJD(4) * t281;
t228 = qJD(4) * t239;
t102 = qJD(3) * t138 + t206 * t227 + t209 * t228;
t55 = -t133 * qJ(5) - qJD(5) * t221 + t102;
t21 = t205 * t55 - t257 * (t103 + t220);
t22 = t257 * t55 + (-t206 * t228 + t209 * t227 + t216 + t220) * t205;
t67 = t117 * t257 + t205 * t214;
t237 = t21 * t66 + t67 * t22;
t236 = t257 * pkin(4);
t235 = t207 * t246;
t50 = -t100 * t257 + t205 * t99;
t51 = t205 * t100 + t257 * t99;
t18 = t50 * mrSges(6,1) + t51 * mrSges(6,2);
t87 = -t132 * t205 + t133 * t257;
t32 = t87 * mrSges(6,1) + t88 * mrSges(6,2);
t17 = t50 * mrSges(7,1) - t51 * mrSges(7,3);
t31 = t87 * mrSges(7,1) - t88 * mrSges(7,3);
t231 = (2 * Ifges(3,4)) + t263;
t119 = pkin(4) * t133 + t201;
t130 = pkin(4) * t156 + t179;
t40 = -mrSges(7,1) * t249 + t51 * mrSges(7,2);
t226 = -mrSges(4,1) * t210 + mrSges(4,2) * t207;
t225 = mrSges(4,1) * t207 + mrSges(4,2) * t210;
t224 = Ifges(4,1) * t210 - t265;
t223 = -Ifges(4,2) * t207 + t264;
t181 = Ifges(4,2) * t210 + t265;
t222 = Ifges(4,5) * t207 + Ifges(4,6) * t210;
t148 = pkin(4) * t221 + t195;
t69 = -pkin(4) * t100 + t144;
t5 = t10 * t257 - t205 * t12;
t29 = -t205 * t68 + t257 * t63;
t158 = -t205 * t206 * pkin(3) + t194 * t257;
t219 = t232 - t235;
t217 = -Ifges(5,5) * t99 - Ifges(5,6) * t100 + t249 * t286 - t282 * t50 - t283 * t51;
t128 = Ifges(5,6) * t133;
t129 = Ifges(5,5) * t132;
t80 = Ifges(7,6) * t87;
t81 = Ifges(6,6) * t87;
t82 = Ifges(6,5) * t88;
t83 = Ifges(7,4) * t88;
t213 = t103 * mrSges(5,1) - t102 * mrSges(5,2) - t128 - t129 + t80 - t81 + t82 + t83 + (-mrSges(6,2) + mrSges(7,3)) * t22 + t284 * t21;
t2 = qJ(6) * t249 - qJD(6) * t211 + t6;
t3 = -pkin(5) * t249 - t5;
t212 = t26 * mrSges(5,1) + t5 * mrSges(6,1) - t3 * mrSges(7,1) - t25 * mrSges(5,2) - t6 * mrSges(6,2) + t2 * mrSges(7,3) - t217;
t204 = qJD(6) * mrSges(7,3);
t200 = Ifges(4,5) * t245;
t193 = -t236 - pkin(5);
t192 = qJ(6) + t269;
t177 = -mrSges(4,1) * t211 - mrSges(4,3) * t254;
t176 = mrSges(4,2) * t211 - mrSges(4,3) * t255;
t175 = t224 * qJD(3);
t174 = t223 * qJD(3);
t173 = t225 * qJD(3);
t153 = -pkin(5) - t158;
t152 = qJ(6) + t159;
t151 = -Ifges(4,5) * t211 + t208 * t224;
t150 = t208 * t223 - t262;
t146 = -t211 * t268 + t168;
t143 = -mrSges(4,2) * t249 - mrSges(4,3) * t218;
t142 = mrSges(4,1) * t249 - mrSges(4,3) * t219;
t141 = -mrSges(5,1) * t211 + t157 * mrSges(5,3);
t140 = mrSges(5,2) * t211 - t156 * mrSges(5,3);
t136 = Ifges(5,1) * t170 - Ifges(5,4) * t221;
t135 = Ifges(5,4) * t170 - Ifges(5,2) * t221;
t134 = mrSges(5,1) * t221 + mrSges(5,2) * t170;
t125 = t170 * t205 + t221 * t257;
t123 = mrSges(4,1) * t218 + mrSges(4,2) * t219;
t118 = mrSges(5,1) * t156 - mrSges(5,2) * t157;
t116 = -t182 * t246 + (Ifges(4,5) * t208 + t211 * t224) * qJD(2);
t115 = -t181 * t246 + (Ifges(4,6) * t208 + t211 * t223) * qJD(2);
t112 = t156 * t257 - t157 * t205;
t111 = -Ifges(5,1) * t157 - Ifges(5,4) * t156 - Ifges(5,5) * t211;
t110 = -Ifges(5,4) * t157 - Ifges(5,2) * t156 - Ifges(5,6) * t211;
t109 = -qJD(3) * t147 + t250;
t105 = mrSges(6,2) * t211 - t112 * mrSges(6,3);
t104 = -t112 * mrSges(7,2) - mrSges(7,3) * t211;
t94 = -Ifges(5,1) * t132 - Ifges(5,4) * t133;
t93 = -Ifges(5,4) * t132 - Ifges(5,2) * t133;
t92 = mrSges(5,1) * t133 - mrSges(5,2) * t132;
t86 = -mrSges(5,2) * t249 + mrSges(5,3) * t100;
t85 = mrSges(5,1) * t249 - mrSges(5,3) * t99;
t77 = Ifges(6,1) * t126 - Ifges(6,4) * t125;
t76 = Ifges(7,1) * t126 + Ifges(7,5) * t125;
t75 = Ifges(6,4) * t126 - Ifges(6,2) * t125;
t74 = Ifges(7,5) * t126 + Ifges(7,3) * t125;
t73 = mrSges(6,1) * t125 + mrSges(6,2) * t126;
t72 = mrSges(7,1) * t125 - mrSges(7,3) * t126;
t65 = pkin(5) * t125 - qJ(6) * t126 + t148;
t62 = mrSges(6,1) * t112 + mrSges(6,2) * t113;
t61 = mrSges(7,1) * t112 - mrSges(7,3) * t113;
t59 = Ifges(6,1) * t113 - Ifges(6,4) * t112 - Ifges(6,5) * t211;
t58 = Ifges(7,1) * t113 - Ifges(7,4) * t211 + Ifges(7,5) * t112;
t57 = Ifges(6,4) * t113 - Ifges(6,2) * t112 - Ifges(6,6) * t211;
t56 = Ifges(7,5) * t113 - Ifges(7,6) * t211 + Ifges(7,3) * t112;
t53 = pkin(5) * t112 - qJ(6) * t113 + t130;
t52 = -mrSges(5,1) * t100 + mrSges(5,2) * t99;
t42 = Ifges(5,1) * t99 + Ifges(5,4) * t100 + Ifges(5,5) * t249;
t41 = Ifges(5,4) * t99 + Ifges(5,2) * t100 + Ifges(5,6) * t249;
t39 = mrSges(6,1) * t249 - mrSges(6,3) * t51;
t38 = -mrSges(6,2) * t249 - mrSges(6,3) * t50;
t37 = -mrSges(7,2) * t50 + mrSges(7,3) * t249;
t36 = Ifges(6,1) * t88 - Ifges(6,4) * t87;
t35 = Ifges(7,1) * t88 + Ifges(7,5) * t87;
t34 = Ifges(6,4) * t88 - Ifges(6,2) * t87;
t33 = Ifges(7,5) * t88 + Ifges(7,3) * t87;
t28 = t211 * pkin(5) - t29;
t27 = -qJ(6) * t211 + t30;
t23 = pkin(5) * t87 - qJ(6) * t88 - qJD(6) * t126 + t119;
t16 = Ifges(6,1) * t51 - Ifges(6,4) * t50 + Ifges(6,5) * t249;
t15 = Ifges(7,1) * t51 + Ifges(7,4) * t249 + Ifges(7,5) * t50;
t14 = Ifges(6,4) * t51 - Ifges(6,2) * t50 + Ifges(6,6) * t249;
t13 = Ifges(7,5) * t51 + Ifges(7,6) * t249 + Ifges(7,3) * t50;
t7 = pkin(5) * t50 - qJ(6) * t51 - qJD(6) * t113 + t69;
t1 = [((mrSges(3,2) * t274 - t207 * t150 + t210 * t151 + t211 * t231) * qJD(2) + t217 + t280) * t211 + (t56 - t57) * t50 + 0.2e1 * t108 * t176 + 0.2e1 * t109 * t177 + 0.2e1 * t179 * t52 - t157 * t42 + (t123 * t273 - t207 * t115 + t210 * t116 + (-t210 * t150 - t207 * t151 + t211 * t222) * qJD(3) + (-Ifges(5,5) * t157 - Ifges(5,6) * t156 + mrSges(3,1) * t274 + (Ifges(4,5) * t210 - t231) * t208 + t283 * t113 + t282 * t112 + (pkin(7) ^ 2 * t278 + t225 * t273 + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3) + t286) * t211) * qJD(2)) * t208 - t156 * t41 + 0.2e1 * t25 * t140 + 0.2e1 * t26 * t141 + 0.2e1 * t144 * t118 + 0.2e1 * t146 * t142 + 0.2e1 * t147 * t143 + 0.2e1 * t130 * t18 + (t58 + t59) * t51 + (t2 * t27 + t28 * t3 + t53 * t7) * t275 + (t130 * t69 + t29 * t5 + t30 * t6) * t276 + (t144 * t179 + t25 * t90 + t26 * t89) * t277 + (t147 * t108 + t146 * t109) * t278 + (t13 - t14) * t112 + (t15 + t16) * t113 + 0.2e1 * t27 * t37 + 0.2e1 * t30 * t38 + 0.2e1 * t29 * t39 + 0.2e1 * t28 * t40 + 0.2e1 * t53 * t17 + 0.2e1 * t7 * t61 + 0.2e1 * t69 * t62 + 0.2e1 * t89 * t85 + 0.2e1 * t90 * t86 + 0.2e1 * t2 * t104 + 0.2e1 * t6 * t105 + 0.2e1 * t5 * t106 + 0.2e1 * t3 * t107 + t100 * t110 + t99 * t111; (-pkin(8) * t142 - t109 * mrSges(4,3) + t116 / 0.2e1) * t207 + (t40 - t39) * t66 + (t74 / 0.2e1 - t75 / 0.2e1) * t50 + (t76 / 0.2e1 + t77 / 0.2e1) * t51 + (pkin(8) * t143 + t108 * mrSges(4,3) + t115 / 0.2e1) * t210 + (-t200 / 0.2e1 + t129 / 0.2e1 + t128 / 0.2e1 - t82 / 0.2e1 + t81 / 0.2e1 - t83 / 0.2e1 - t80 / 0.2e1 + (Ifges(3,5) + t181 * t270 + t253 / 0.2e1 + (-mrSges(3,1) + t226) * pkin(7)) * qJD(2)) * t211 + (t104 + t105) * t22 + (-pkin(2) * t202 + (t108 * t210 - t109 * t207 + (-t146 * t210 - t147 * t207) * qJD(3)) * pkin(8)) * m(4) + (t35 / 0.2e1 + t36 / 0.2e1) * t113 + (t33 / 0.2e1 - t34 / 0.2e1) * t112 + (-t125 * t6 - t126 * t5 - t29 * t88 - t30 * t87) * mrSges(6,3) + (-t125 * t2 + t126 * t3 - t27 * t87 + t28 * t88) * mrSges(7,2) + (t132 * t89 - t133 * t90 - t170 * t26 - t221 * t25) * mrSges(5,3) - t221 * t41 / 0.2e1 + (t210 * t175 / 0.2e1 + t174 * t270 - Ifges(3,6) * qJD(2) + (t182 * t270 - t210 * t181 / 0.2e1) * qJD(3) + (mrSges(3,2) * qJD(2) + t173) * pkin(7) + (Ifges(5,5) * t170 - Ifges(5,6) * t221 + t282 * t125 + t283 * t126 + t222) * qJD(2) / 0.2e1) * t208 + (t37 + t38) * t67 + m(6) * (t119 * t130 + t148 * t69 - t21 * t29 + t22 * t30 - t5 * t66 + t6 * t67) + m(7) * (t2 * t67 + t21 * t28 + t22 * t27 + t23 * t53 + t3 * t66 + t65 * t7) + t195 * t52 + t170 * t42 / 0.2e1 + t179 * t92 + (t58 / 0.2e1 + t59 / 0.2e1) * t88 - t157 * t94 / 0.2e1 - t156 * t93 / 0.2e1 + t138 * t85 + t102 * t140 + t103 * t141 + t144 * t134 + t148 * t18 + t130 * t32 + (t56 / 0.2e1 - t57 / 0.2e1) * t87 - t132 * t111 / 0.2e1 - t133 * t110 / 0.2e1 + t100 * t135 / 0.2e1 + t99 * t136 / 0.2e1 + m(5) * (t102 * t90 + t103 * t89 + t138 * t26 + t144 * t195 + t179 * t201 - t25 * t285) - t285 * t86 + t251 * t21 + ((t151 / 0.2e1 - pkin(8) * t177 - t146 * mrSges(4,3)) * t210 + (-pkin(8) * t176 - t147 * mrSges(4,3) + pkin(3) * t118 - t150 / 0.2e1 + t262 / 0.2e1) * t207) * qJD(3) + (t13 / 0.2e1 - t14 / 0.2e1) * t125 + (t15 / 0.2e1 + t16 / 0.2e1) * t126 + t53 * t31 + t23 * t61 + t65 * t17 + t7 * t72 + t69 * t73 + t119 * t62 - pkin(2) * t123; -0.2e1 * pkin(2) * t173 + 0.2e1 * t119 * t73 - t132 * t136 - t133 * t135 + 0.2e1 * t148 * t32 - t221 * t93 + t170 * t94 + t210 * t174 + t207 * t175 + 0.2e1 * t195 * t92 + 0.2e1 * t23 * t72 + 0.2e1 * t65 * t31 + (t76 + t77) * t88 + (t74 - t75) * t87 + (t35 + t36) * t126 + (t33 - t34) * t125 + (t253 + (0.2e1 * pkin(3) * t134 - t181) * t207) * qJD(3) + (-t102 * t285 + t103 * t138 + t195 * t201) * t277 + (t119 * t148 + t237) * t276 + (t23 * t65 + t237) * t275 + 0.2e1 * (-t102 * t221 - t103 * t170 + t132 * t138 + t133 * t285) * mrSges(5,3) + 0.2e1 * (mrSges(6,3) + mrSges(7,2)) * (-t125 * t22 + t126 * t21 + t66 * t88 - t67 * t87); -t218 * Ifges(4,6) + t212 + t158 * t39 + t159 * t38 + t149 * t104 + t152 * t37 + t153 * t40 + t155 * t105 - Ifges(4,5) * t235 + m(6) * (-t154 * t29 + t155 * t30 + t158 * t5 + t159 * t6) + m(7) * (t149 * t27 + t152 * t2 + t153 * t3 + t154 * t28) + t251 * t154 + (-t141 * t244 + m(5) * (t206 * t25 + t209 * t26 + t243 * t90 - t244 * t89) + t209 * t85 + t206 * t86 + t140 * t243) * pkin(3) - t108 * mrSges(4,2) + t109 * mrSges(4,1) - t280; t200 + (pkin(8) * t226 - t263) * qJD(3) + t213 + m(6) * (t155 * t67 - t158 * t21 + t159 * t22 + t259) + m(7) * (t149 * t67 + t152 * t22 + t153 * t21 + t259) + (-t125 * t155 - t158 * t88 - t159 * t87 + t256) * mrSges(6,3) + (-t125 * t149 - t152 * t87 + t153 * t88 + t256) * mrSges(7,2) + (m(5) * (t102 * t206 + t103 * t209 + (-t138 * t206 - t209 * t285) * qJD(4)) + (t209 * t132 - t206 * t133 + (t170 * t206 - t209 * t221) * qJD(4)) * mrSges(5,3)) * pkin(3); -0.2e1 * t258 + 0.2e1 * t260 + (-t154 * t158 + t155 * t159) * t276 + (t149 * t152 + t153 * t154) * t275 + 0.2e1 * t284 * t154 + 0.2e1 * (-mrSges(5,1) * t206 - mrSges(5,2) * t209) * t261; t212 + m(7) * (qJD(6) * t27 + t192 * t2 + t193 * t3) + t38 * t269 + t39 * t236 + t192 * t37 + t193 * t40 + (t205 * t6 + t257 * t5) * t272 + qJD(6) * t104; t213 + t193 * t267 + (t205 * t22 - t21 * t257) * t272 + m(7) * (qJD(6) * t67 + t192 * t22 + t193 * t21) + (-t236 * t88 - t269 * t87) * mrSges(6,3) + (-qJD(6) * t125 - t192 * t87) * mrSges(7,2); t155 * t205 * t272 + t204 + m(7) * (qJD(6) * t152 + t149 * t192) + t260 - t258 - mrSges(5,2) * t241 - mrSges(5,1) * t242 + (m(7) * t193 - t257 * t272 + t284) * t154; 0.2e1 * m(7) * qJD(6) * t192 + 0.2e1 * t204; m(6) * t69 + m(7) * t7 + t17 + t18; m(6) * t119 + m(7) * t23 + t31 + t32; 0; 0; 0; m(7) * t3 + t40; m(7) * t21 + t267; m(7) * t154; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
