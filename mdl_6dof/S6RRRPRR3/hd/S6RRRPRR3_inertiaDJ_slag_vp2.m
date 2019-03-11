% Calculate time derivative of joint inertia matrix for
% S6RRRPRR3
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
% Datum: 2019-03-09 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:10:57
% EndTime: 2019-03-09 18:11:06
% DurationCPUTime: 4.30s
% Computational Cost: add. (6998->393), mult. (14709->548), div. (0->0), fcn. (13830->8), ass. (0->177)
t136 = cos(qJ(6));
t287 = -t136 / 0.2e1;
t132 = sin(qJ(6));
t157 = mrSges(7,1) * t132 + mrSges(7,2) * t136;
t100 = t157 * qJD(6);
t133 = sin(qJ(5));
t137 = cos(qJ(5));
t134 = sin(qJ(3));
t247 = cos(qJ(2));
t252 = -pkin(8) - pkin(7);
t151 = t252 * t247;
t135 = sin(qJ(2));
t246 = cos(qJ(3));
t176 = t246 * t135;
t98 = t252 * t176;
t69 = -t134 * t151 - t98;
t96 = t134 * t247 + t176;
t142 = -t96 * pkin(9) + t69;
t201 = t134 * t135;
t149 = t246 * t247 - t201;
t273 = t252 * t201;
t99 = t246 * t151;
t70 = -t99 + t273;
t56 = -pkin(9) * t149 + t70;
t36 = t133 * t56 - t137 * t142;
t208 = qJD(5) * t36;
t147 = qJD(2) * t151;
t267 = qJD(2) + qJD(3);
t46 = -qJD(3) * t99 - t246 * t147 + t267 * t273;
t67 = t267 * t149;
t274 = -t67 * pkin(9) + t46;
t45 = qJD(2) * t98 - t69 * qJD(3) + t134 * t147;
t68 = t267 * t96;
t35 = t68 * pkin(9) + t45;
t12 = t274 * t133 + t137 * t35 - t208;
t193 = qJD(6) * t132;
t123 = Ifges(7,6) * t193;
t192 = qJD(6) * t136;
t120 = -t247 * pkin(2) - pkin(1);
t57 = -pkin(3) * t149 - t96 * qJ(4) + t120;
t49 = pkin(4) * t149 - t57;
t59 = t133 * t96 + t137 * t149;
t60 = -t133 * t149 + t137 * t96;
t20 = pkin(5) * t59 - pkin(10) * t60 + t49;
t37 = t133 * t142 + t137 * t56;
t19 = t132 * t20 + t136 * t37;
t204 = qJD(6) * t19;
t196 = qJD(2) * t135;
t186 = pkin(2) * t196;
t39 = t68 * pkin(3) - t67 * qJ(4) - t96 * qJD(4) + t186;
t21 = -pkin(4) * t68 - t39;
t31 = -qJD(5) * t59 + t133 * t68 + t137 * t67;
t206 = qJD(5) * t60;
t32 = t133 * t67 - t137 * t68 + t206;
t8 = pkin(5) * t32 - pkin(10) * t31 + t21;
t3 = -t12 * t132 + t136 * t8 - t204;
t242 = t132 * t3;
t18 = -t132 * t37 + t136 * t20;
t283 = -t18 * t192 - t19 * t193;
t286 = -(Ifges(6,6) - Ifges(7,5) * t132 / 0.2e1 + Ifges(7,6) * t287) * t32 - t12 * mrSges(6,2) + Ifges(6,5) * t31 + (-t242 + t283) * mrSges(7,3) + t36 * t100 + t59 * (Ifges(7,5) * t192 - t123) / 0.2e1;
t285 = -mrSges(4,1) - mrSges(5,1);
t109 = -t136 * mrSges(7,1) + mrSges(7,2) * t132;
t209 = mrSges(6,1) - t109;
t197 = t132 ^ 2 + t136 ^ 2;
t168 = mrSges(7,3) * t197;
t284 = mrSges(6,2) - t168;
t263 = -m(7) * pkin(5) - t209;
t270 = t136 * t31 - t60 * t193;
t102 = Ifges(7,4) * t192 - Ifges(7,2) * t193;
t177 = t60 * t192;
t154 = t132 * t31 + t177;
t7 = Ifges(7,1) * t270 - Ifges(7,4) * t154 + t32 * Ifges(7,5);
t281 = t60 * t102 / 0.2e1 - t7 / 0.2e1;
t280 = -0.2e1 * t68;
t279 = 0.2e1 * t96;
t277 = Ifges(3,1) - Ifges(3,2);
t240 = t31 * mrSges(6,3);
t9 = mrSges(7,1) * t154 + mrSges(7,2) * t270;
t276 = t240 + t9;
t40 = t157 * t60;
t275 = t60 * mrSges(6,3) + t40;
t163 = t197 * t137;
t187 = t246 * pkin(2);
t119 = -t187 - pkin(3);
t116 = -pkin(4) + t119;
t245 = pkin(2) * t134;
t117 = qJ(4) + t245;
t74 = t133 * t116 + t137 * t117;
t138 = -pkin(3) - pkin(4);
t107 = t137 * qJ(4) + t133 * t138;
t223 = t132 * t60;
t41 = -mrSges(7,2) * t59 - mrSges(7,3) * t223;
t216 = t136 * t60;
t42 = mrSges(7,1) * t59 - mrSges(7,3) * t216;
t272 = -t42 * t192 - t41 * t193;
t103 = Ifges(7,1) * t192 - Ifges(7,4) * t193;
t203 = t132 * t103;
t271 = t136 * t102 + t203;
t266 = m(7) * pkin(10) + mrSges(7,3);
t227 = pkin(2) * qJD(3);
t264 = (-mrSges(4,2) * t246 + t285 * t134) * t227;
t219 = t136 * t19;
t236 = t59 * mrSges(6,3);
t260 = t136 * t41 - t132 * t42 - t236 + m(7) * (-t132 * t18 + t219);
t259 = 2 * m(5);
t258 = 2 * m(6);
t257 = 0.2e1 * m(7);
t254 = -0.2e1 * t100;
t253 = t31 / 0.2e1;
t205 = qJD(6) * t18;
t2 = t12 * t136 + t132 * t8 + t205;
t251 = t2 * mrSges(7,3);
t230 = Ifges(7,4) * t132;
t111 = Ifges(7,2) * t136 + t230;
t249 = -t111 / 0.2e1;
t207 = qJD(5) * t37;
t13 = t133 * t35 - t137 * t274 + t207;
t243 = t13 * t36;
t241 = t136 * t2;
t239 = t32 * mrSges(6,3);
t160 = qJD(3) * t187;
t115 = t160 + qJD(4);
t185 = t134 * t227;
t55 = t74 * qJD(5) + t115 * t133 - t137 * t185;
t238 = t36 * t55;
t76 = t133 * qJD(4) + t107 * qJD(5);
t237 = t36 * t76;
t234 = t67 * mrSges(5,2);
t229 = Ifges(7,4) * t136;
t228 = Ifges(7,5) * t136;
t72 = -pkin(10) + t74;
t222 = t132 * t72;
t106 = -t133 * qJ(4) + t137 * t138;
t75 = t137 * qJD(4) + qJD(5) * t106;
t221 = t132 * t75;
t15 = -mrSges(7,2) * t32 - mrSges(7,3) * t154;
t220 = t136 * t15;
t215 = t137 * t55;
t214 = t137 * t76;
t211 = t60 * t103;
t112 = Ifges(7,1) * t132 + t229;
t210 = t60 * t112;
t202 = t132 * t111;
t199 = t136 * t112;
t194 = qJD(5) * t137;
t195 = qJD(5) * t133;
t198 = -mrSges(6,1) * t195 - mrSges(6,2) * t194;
t173 = t70 * t45 + t46 * t69;
t172 = qJD(2) * t247;
t170 = t193 / 0.2e1;
t169 = -t192 / 0.2e1;
t167 = t75 * t197;
t73 = t116 * t137 - t117 * t133;
t54 = qJD(5) * t73 + t115 * t137 + t133 * t185;
t166 = t197 * t54;
t105 = -pkin(10) + t107;
t165 = t105 * t197;
t164 = t133 * t197;
t162 = t96 * t185;
t161 = 0.2e1 * t209;
t158 = t241 - t242;
t14 = mrSges(7,1) * t32 - mrSges(7,3) * t270;
t156 = -t132 * t14 + t220;
t152 = 0.2e1 * t284;
t91 = t111 * t193;
t150 = t111 * t170 + t91 / 0.2e1 + 0.2e1 * t112 * t169 - t271;
t146 = t203 + (qJD(6) * t112 + t102) * t136 - t91;
t144 = (-qJD(5) * t168 + t100) * t137 - t109 * t195 - t198;
t143 = t158 + t283;
t140 = t270 * Ifges(7,5) - Ifges(7,6) * t154 + Ifges(7,3) * t32;
t27 = t59 * Ifges(7,6) + (-Ifges(7,2) * t132 + t229) * t60;
t28 = t59 * Ifges(7,5) + (Ifges(7,1) * t136 - t230) * t60;
t6 = Ifges(7,4) * t270 - Ifges(7,2) * t154 + t32 * Ifges(7,6);
t139 = t28 * t169 + t202 * t253 + t111 * t177 / 0.2e1 - t31 * t199 / 0.2e1 + (Ifges(5,6) - Ifges(4,6)) * t68 + (Ifges(5,4) + Ifges(4,5)) * t67 + t285 * t46 + (mrSges(5,3) - mrSges(4,2)) * t45 + (t211 + t6) * t287 + (t27 + t210) * t170 + t281 * t132 + t209 * t13 - t286;
t104 = pkin(5) - t106;
t71 = pkin(5) - t73;
t1 = [(Ifges(5,3) + Ifges(4,2)) * t149 * t280 + (Ifges(5,1) + Ifges(4,1)) * t67 * t279 + 0.2e1 * (-mrSges(4,1) * t149 + mrSges(4,2) * t96) * t186 + 0.2e1 * t39 * (-mrSges(5,1) * t149 - mrSges(5,3) * t96) + 0.2e1 * (Ifges(4,4) - Ifges(5,5)) * (t149 * t67 - t68 * t96) + (mrSges(4,3) + mrSges(5,2)) * (0.2e1 * t149 * t45 + t46 * t279 + t70 * t280 + 0.2e1 * t67 * t69) - 0.2e1 * pkin(1) * (mrSges(3,1) * t135 + mrSges(3,2) * t247) * qJD(2) + t270 * t28 - t154 * t27 + (0.2e1 * Ifges(3,4) * t247 + t277 * t135) * t172 + (-0.2e1 * Ifges(3,4) * t135 + t277 * t247) * t196 + (t18 * t3 + t19 * t2 + t243) * t257 + (t12 * t37 + t21 * t49 + t243) * t258 + (t39 * t57 + t173) * t259 + 0.2e1 * (-t31 * t59 - t32 * t60) * Ifges(6,4) + t59 * t140 + t7 * t216 + 0.2e1 * t275 * t13 + 0.2e1 * t276 * t36 + 0.2e1 * m(4) * (t120 * t186 + t173) + 0.2e1 * t18 * t14 + 0.2e1 * t19 * t15 + 0.2e1 * t2 * t41 + 0.2e1 * t3 * t42 + 0.2e1 * t49 * (mrSges(6,1) * t32 + mrSges(6,2) * t31) + 0.2e1 * t60 * t31 * Ifges(6,1) + 0.2e1 * t59 * Ifges(6,2) * t32 + 0.2e1 * t21 * (mrSges(6,1) * t59 + mrSges(6,2) * t60) + 0.2e1 * t57 * (mrSges(5,1) * t68 - mrSges(5,3) * t67) - t6 * t223 + 0.2e1 * t120 * (mrSges(4,1) * t68 + mrSges(4,2) * t67) + t32 * (Ifges(7,3) * t59 + (-Ifges(7,6) * t132 + t228) * t60) - 0.2e1 * t12 * t236 - 0.2e1 * t37 * t239; m(4) * (-t246 * t46 + t134 * t45 + (t134 * t69 + t246 * t70) * qJD(3)) * pkin(2) + t139 + t119 * t234 + Ifges(3,5) * t172 + m(5) * (t115 * t70 + t117 * t45 + t119 * t46 + t185 * t69) - Ifges(3,6) * t196 + t71 * t9 - t14 * t222 + m(6) * (t12 * t74 - t13 * t73 + t238) + m(7) * (t13 * t71 + t238) - t74 * t239 - t73 * t240 - mrSges(7,3) * t241 + t275 * t55 + (-mrSges(3,1) * t172 + mrSges(3,2) * t196) * pkin(7) + (t220 + m(7) * ((-t132 * t19 - t136 * t18) * qJD(6) + t158) + t272) * t72 + (m(6) * t37 + t260) * t54 + (t115 * t149 - t117 * t68 + t162) * mrSges(5,2) + (t149 * t160 - t187 * t67 - t245 * t68 + t162) * mrSges(4,3); 0.2e1 * t115 * mrSges(5,3) + t71 * t254 + t55 * t161 + t152 * t54 + 0.2e1 * t264 + (t166 * t72 + t55 * t71) * t257 + (t54 * t74 - t55 * t73) * t258 + (t115 * t117 + t119 * t185) * t259 + t146; (m(7) * t143 + t156 + t272) * t105 + t139 + (t75 * t41 - t251) * t136 + m(7) * (t104 * t13 - t18 * t221 + t219 * t75 + t237) - t42 * t221 + m(6) * (-t106 * t13 + t107 * t12 + t37 * t75 + t237) + m(5) * (-pkin(3) * t46 + qJ(4) * t45 + qJD(4) * t70) + (-t106 * t31 - t107 * t32 - t59 * t75 + t60 * t76) * mrSges(6,3) + (-pkin(3) * t67 - qJ(4) * t68 + qJD(4) * t149) * mrSges(5,2) + t76 * t40 + t104 * t9; (-t104 - t71) * t100 + (qJD(4) + t115) * mrSges(5,3) + t264 + m(7) * (t104 * t55 + t165 * t54 + t167 * t72 + t71 * t76) + m(6) * (-t106 * t55 + t107 * t54 - t73 * t76 + t74 * t75) + m(5) * (-pkin(3) * t185 + qJ(4) * t115 + qJD(4) * t117) + t146 + t209 * (t76 + t55) + t284 * (t75 + t54); t104 * t254 + t76 * t161 + 0.2e1 * (m(5) * qJ(4) + mrSges(5,3)) * qJD(4) + t152 * t75 + (t104 * t76 + t165 * t75) * t257 + (-t106 * t76 + t107 * t75) * t258 + t146; t234 + m(5) * t46 + (-m(7) * t13 + m(6) * (-t13 + t207) + t260 * qJD(5) - t276) * t137 + (qJD(5) * t40 + (-t132 * t41 - t136 * t42) * qJD(6) + (-t32 + t206) * mrSges(6,3) + m(7) * (t143 + t208) + m(6) * (t12 + t208) + t156) * t133; m(5) * t185 + m(7) * (-t215 + t54 * t164 + (t133 * t71 + t163 * t72) * qJD(5)) + m(6) * (t133 * t54 - t215 + (-t133 * t73 + t137 * t74) * qJD(5)) + t144; m(7) * (-t214 + t75 * t164 + (t104 * t133 + t105 * t163) * qJD(5)) + m(6) * (t133 * t75 - t214 + (-t106 * t133 + t107 * t137) * qJD(5)) + t144; (-0.1e1 + t197) * t133 * t194 * t257; -pkin(5) * t9 + t263 * t13 + (t31 * t249 + (-t210 / 0.2e1 - t27 / 0.2e1) * qJD(6) + (m(7) * (-t3 - t204) - qJD(6) * t41 - t14) * pkin(10) - t281) * t132 + (t112 * t253 + t251 + t211 / 0.2e1 + t6 / 0.2e1 + (t60 * t249 + t28 / 0.2e1) * qJD(6) + (m(7) * (t2 - t205) + t15 - qJD(6) * t42) * pkin(10)) * t136 + t286; -t54 * mrSges(6,2) + (t71 + pkin(5)) * t100 + t150 + t266 * t166 + t263 * t55; -t75 * mrSges(6,2) + (pkin(5) + t104) * t100 + t150 + t266 * t167 + t263 * t76; -t137 * t100 + (t133 * t109 + m(7) * (-pkin(5) * t133 + pkin(10) * t163) + mrSges(7,3) * t163) * qJD(5) + t198; pkin(5) * t254 + (t199 - t202) * qJD(6) + t271; mrSges(7,1) * t3 - mrSges(7,2) * t2 + t140; t123 - t157 * t54 + (mrSges(7,2) * t222 + (-mrSges(7,1) * t72 - Ifges(7,5)) * t136) * qJD(6); t123 - t157 * t75 + (t105 * t109 - t228) * qJD(6); (t133 * t193 - t136 * t194) * mrSges(7,2) + (-t132 * t194 - t133 * t192) * mrSges(7,1); -t123 + (pkin(10) * t109 + t228) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
