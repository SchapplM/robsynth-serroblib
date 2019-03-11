% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% Datum: 2019-03-08 21:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRP5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:45:33
% EndTime: 2019-03-08 21:45:45
% DurationCPUTime: 5.87s
% Computational Cost: add. (3218->492), mult. (7913->646), div. (0->0), fcn. (4745->8), ass. (0->235)
t296 = Ifges(6,1) + Ifges(7,1);
t292 = Ifges(7,4) + Ifges(6,5);
t151 = sin(qJ(5));
t154 = cos(qJ(5));
t155 = cos(qJ(3));
t227 = qJD(2) * t155;
t108 = qJD(3) * t151 + t154 * t227;
t205 = t151 * t227;
t109 = qJD(3) * t154 - t205;
t148 = qJD(3) * qJ(4);
t153 = sin(qJ(2));
t149 = sin(pkin(6));
t232 = qJD(1) * t149;
t209 = t153 * t232;
t117 = qJD(2) * pkin(8) + t209;
t150 = cos(pkin(6));
t152 = sin(qJ(3));
t231 = qJD(1) * t152;
t133 = t150 * t231;
t85 = t155 * t117 + t133;
t65 = pkin(4) * t227 + t85;
t52 = t148 + t65;
t16 = pkin(5) * t108 - qJ(6) * t109 + t52;
t157 = -pkin(3) - pkin(9);
t241 = t150 * t155;
t134 = qJD(1) * t241;
t200 = pkin(4) * qJD(2) + t117;
t174 = t200 * t152;
t64 = t134 - t174;
t44 = qJD(3) * t157 + qJD(4) - t64;
t201 = -qJ(4) * t152 - pkin(2);
t103 = t155 * t157 + t201;
t156 = cos(qJ(2));
t208 = t156 * t232;
t66 = qJD(2) * t103 - t208;
t14 = -t151 * t66 + t154 * t44;
t15 = t151 * t44 + t154 * t66;
t177 = t14 * t151 - t15 * t154;
t247 = Ifges(7,5) * t151;
t181 = -Ifges(7,3) * t154 + t247;
t250 = Ifges(6,4) * t151;
t185 = Ifges(6,2) * t154 + t250;
t192 = mrSges(7,1) * t154 + mrSges(7,3) * t151;
t194 = mrSges(6,1) * t154 - mrSges(6,2) * t151;
t229 = qJD(2) * t152;
t139 = qJD(5) + t229;
t277 = qJD(6) - t14;
t7 = -pkin(5) * t139 + t277;
t8 = qJ(6) * t139 + t15;
t195 = t151 * t7 + t154 * t8;
t259 = t154 / 0.2e1;
t260 = -t154 / 0.2e1;
t262 = -t151 / 0.2e1;
t264 = -t139 / 0.2e1;
t266 = -t109 / 0.2e1;
t267 = t108 / 0.2e1;
t268 = -t108 / 0.2e1;
t105 = Ifges(6,4) * t108;
t248 = Ifges(7,5) * t108;
t278 = t109 * t296 + t292 * t139 - t105 + t248;
t279 = Ifges(6,6) - Ifges(7,6);
t246 = Ifges(7,5) * t154;
t249 = Ifges(6,4) * t154;
t288 = t151 * t296 - t246 + t249;
t104 = Ifges(7,5) * t109;
t35 = t139 * Ifges(7,6) + t108 * Ifges(7,3) + t104;
t251 = Ifges(6,4) * t109;
t38 = -t108 * Ifges(6,2) + t139 * Ifges(6,6) + t251;
t297 = t259 * t35 + t260 * t38 + t262 * t278 + t16 * t192 + t181 * t268 + t185 * t267 + t52 * t194 + (t151 * t292 + t154 * t279) * t264 + t288 * t266 - t195 * mrSges(7,2) + t177 * mrSges(6,3);
t245 = qJD(2) * pkin(2);
t118 = -t208 - t245;
t282 = qJD(3) / 0.2e1;
t283 = -qJD(3) / 0.2e1;
t284 = qJD(2) / 0.2e1;
t68 = -t148 - t85;
t120 = -pkin(3) * t155 + t201;
t86 = qJD(2) * t120 - t208;
t295 = t85 * mrSges(4,3) + t86 * mrSges(5,2) + Ifges(4,6) * t282 + (Ifges(4,4) * t152 + t155 * Ifges(4,2)) * t284 + Ifges(5,5) * t283 - (-Ifges(5,6) * t152 - t155 * Ifges(5,3)) * qJD(2) / 0.2e1 - t118 * mrSges(4,1) - t68 * mrSges(5,1) + t297;
t269 = pkin(4) + pkin(8);
t294 = -t227 / 0.2e1;
t281 = mrSges(5,1) + mrSges(4,3);
t293 = mrSges(5,2) - mrSges(4,1);
t221 = qJD(2) * qJD(3);
t202 = t155 * t221;
t220 = qJD(3) * qJD(5);
t222 = qJD(5) * t154;
t226 = qJD(3) * t152;
t75 = -t151 * t220 + (t151 * t226 - t155 * t222) * qJD(2);
t203 = t152 * t221;
t76 = -qJD(5) * t205 + (-t203 + t220) * t154;
t291 = (-Ifges(6,4) + Ifges(7,5)) * t76 + t296 * t75 + t292 * t202;
t84 = t117 * t152 - t134;
t290 = -qJD(4) - t84;
t289 = t154 * t296 + t247 - t250;
t287 = -Ifges(4,1) / 0.2e1;
t285 = Ifges(4,4) * t294;
t230 = qJD(2) * t149;
t204 = qJD(1) * t230;
t199 = t156 * t204;
t31 = t152 * t199 + (t155 * t200 + t133) * qJD(3);
t178 = pkin(9) * t152 - qJ(4) * t155;
t224 = qJD(4) * t152;
t162 = qJD(3) * t178 - t224;
t233 = pkin(3) * t203 + t153 * t204;
t43 = qJD(2) * t162 + t233;
t4 = -qJD(5) * t15 - t151 * t43 + t154 * t31;
t225 = qJD(3) * t155;
t115 = t269 * t225;
t127 = t269 * t152;
t243 = t154 * t103 + t151 * t127;
t144 = pkin(3) * t226;
t87 = t144 + t162;
t11 = -qJD(5) * t243 + t115 * t154 - t151 * t87;
t212 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t214 = Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t217 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1;
t67 = -qJD(3) * pkin(3) - t290;
t275 = -t212 * t108 - t217 * t109 - t214 * t139 + t15 * mrSges(6,2) + t7 * mrSges(7,1) + t86 * mrSges(5,3) + Ifges(6,6) * t267 + Ifges(7,6) * t268 + t229 * t287 + Ifges(4,5) * t283 + t285 + Ifges(5,4) * t282 + (-t152 * Ifges(5,2) - Ifges(5,6) * t155) * t284 - t118 * mrSges(4,2) - t14 * mrSges(6,1) - t67 * mrSges(5,1) - t8 * mrSges(7,3) - t84 * mrSges(4,3) + t292 * t266 + (Ifges(6,3) + Ifges(7,2)) * t264;
t273 = -m(5) / 0.2e1;
t271 = -t76 / 0.2e1;
t270 = t76 / 0.2e1;
t265 = t109 / 0.2e1;
t261 = t151 / 0.2e1;
t206 = t156 * t230;
t166 = qJD(3) * t150 + t206;
t42 = t117 * t225 + t166 * t231;
t242 = t149 * t153;
t211 = t152 * t242;
t93 = t211 - t241;
t258 = t42 * t93;
t47 = -mrSges(7,2) * t76 + mrSges(7,3) * t202;
t50 = -mrSges(6,2) * t202 - mrSges(6,3) * t76;
t257 = t47 + t50;
t48 = mrSges(6,1) * t202 - mrSges(6,3) * t75;
t49 = -mrSges(7,1) * t202 + t75 * mrSges(7,2);
t256 = t48 - t49;
t253 = mrSges(6,3) * t108;
t79 = -mrSges(6,2) * t139 - t253;
t80 = -mrSges(7,2) * t108 + mrSges(7,3) * t139;
t255 = t79 + t80;
t252 = mrSges(6,3) * t109;
t81 = mrSges(6,1) * t139 - t252;
t82 = -mrSges(7,1) * t139 + mrSges(7,2) * t109;
t254 = -t82 + t81;
t142 = pkin(3) * t229;
t90 = qJD(2) * t178 + t142;
t24 = t151 * t65 + t154 * t90;
t125 = -mrSges(5,1) * t227 - qJD(3) * mrSges(5,3);
t59 = mrSges(6,1) * t108 + mrSges(6,2) * t109;
t244 = -t125 + t59;
t240 = t151 * t156;
t239 = t151 * t157;
t238 = t154 * t157;
t111 = (mrSges(5,2) * t155 - mrSges(5,3) * t152) * qJD(2);
t237 = t111 + (-mrSges(4,1) * t155 + mrSges(4,2) * t152) * qJD(2);
t236 = qJD(3) * t134 + t155 * t199;
t235 = -qJD(3) * t293 - t229 * t281;
t124 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t227;
t234 = t124 - t125;
t128 = t269 * t155;
t228 = qJD(2) * t153;
t223 = qJD(5) * t151;
t219 = pkin(8) * t42 / 0.2e1;
t58 = mrSges(7,1) * t108 - mrSges(7,3) * t109;
t218 = -t58 - t244;
t216 = Ifges(4,5) / 0.2e1 - Ifges(5,4) / 0.2e1;
t215 = Ifges(5,5) / 0.2e1 - Ifges(4,6) / 0.2e1;
t213 = -0.3e1 / 0.2e1 * Ifges(5,6) - 0.3e1 / 0.2e1 * Ifges(4,4);
t210 = t149 * t154 * t156;
t207 = t149 * t228;
t198 = -t208 / 0.2e1;
t3 = t151 * t31 + t154 * t43 + t44 * t222 - t223 * t66;
t1 = qJ(6) * t202 + qJD(6) * t139 + t3;
t2 = -pkin(5) * t202 - t4;
t197 = -t1 * t151 + t154 * t2;
t196 = -t151 * t3 - t154 * t4;
t193 = mrSges(6,1) * t151 + mrSges(6,2) * t154;
t191 = mrSges(7,1) * t151 - mrSges(7,3) * t154;
t186 = -Ifges(6,2) * t151 + t249;
t182 = Ifges(7,3) * t151 + t246;
t180 = pkin(5) * t154 + qJ(6) * t151;
t179 = -pkin(5) * t151 + qJ(6) * t154;
t23 = -t151 * t90 + t154 * t65;
t54 = -t103 * t151 + t127 * t154;
t171 = -pkin(4) - t180;
t41 = -t117 * t226 + t236;
t62 = t149 * t240 + t154 * t93;
t94 = t150 * t152 + t155 * t242;
t168 = -m(4) * t85 + m(5) * t68 - t234;
t167 = -qJ(4) * t225 - t224;
t10 = -t103 * t223 + t151 * t115 + t127 * t222 + t154 * t87;
t165 = t151 * t257 + t154 * t256;
t164 = -t151 * t254 + t154 * t255;
t163 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t147 = qJD(3) * qJD(4);
t28 = -qJD(3) * t174 + t147 + t236;
t158 = qJD(2) ^ 2;
t137 = Ifges(7,2) * t202;
t136 = Ifges(6,3) * t202;
t119 = qJ(4) - t179;
t114 = t269 * t226;
t113 = -qJ(4) * t227 + t142;
t101 = (mrSges(4,1) * t152 + mrSges(4,2) * t155) * t221;
t100 = (-mrSges(5,2) * t152 - mrSges(5,3) * t155) * t221;
t92 = t144 + t167;
t89 = qJD(5) * t180 - qJD(6) * t154 + qJD(4);
t83 = t155 * t180 + t128;
t78 = (t152 * t240 + t153 * t154) * t232;
t77 = -t152 * t154 * t208 + t151 * t209;
t73 = Ifges(7,4) * t75;
t72 = Ifges(6,5) * t75;
t71 = Ifges(6,6) * t76;
t70 = Ifges(7,6) * t76;
t63 = t151 * t93 - t210;
t61 = -qJD(3) * t211 + t155 * t166;
t60 = qJD(3) * t94 + t152 * t206;
t57 = pkin(5) * t109 + qJ(6) * t108;
t56 = qJD(2) * t167 + t233;
t46 = -pkin(5) * t152 - t54;
t45 = qJ(6) * t152 + t243;
t32 = -t147 - t41;
t30 = t134 + (qJD(2) * t171 - t117) * t152;
t27 = mrSges(6,1) * t76 + mrSges(6,2) * t75;
t26 = mrSges(7,1) * t76 - mrSges(7,3) * t75;
t25 = (qJD(5) * t179 + qJD(6) * t151) * t155 + (-pkin(8) + t171) * t226;
t22 = -pkin(5) * t227 - t23;
t21 = qJ(6) * t227 + t24;
t18 = t75 * Ifges(6,4) - t76 * Ifges(6,2) + Ifges(6,6) * t202;
t17 = t75 * Ifges(7,5) + Ifges(7,6) * t202 + t76 * Ifges(7,3);
t13 = qJD(5) * t62 + t151 * t60 + t154 * t207;
t12 = -qJD(5) * t210 - t60 * t154 + (qJD(5) * t93 + t207) * t151;
t9 = -pkin(5) * t225 - t11;
t6 = qJ(6) * t225 + qJD(6) * t152 + t10;
t5 = pkin(5) * t76 - qJ(6) * t75 - qJD(6) * t109 + t28;
t19 = [(t26 + t27) * t94 + t257 * t63 + t256 * t62 - t235 * t60 + t255 * t13 - t254 * t12 + (t124 - t218) * t61 + (-t158 * t153 * mrSges(3,1) + (-mrSges(3,2) * t158 - t100 - t101) * t156) * t149 + (t237 * t242 + t281 * qJD(3) * (-t152 * t94 + t155 * t93)) * qJD(2) + m(4) * (t41 * t94 + t258 + t60 * t84 + t61 * t85 + (t118 - t208) * t207) + m(5) * (-t32 * t94 + t258 + t60 * t67 - t61 * t68 + (-t156 * t56 + t228 * t86) * t149) + m(6) * (-t12 * t14 + t13 * t15 + t28 * t94 + t3 * t63 + t4 * t62 + t52 * t61) + m(7) * (t1 * t63 + t12 * t7 + t13 * t8 + t16 * t61 - t2 * t62 + t5 * t94); t243 * t50 + m(6) * (t10 * t15 + t11 * t14 - t114 * t52 + t128 * t28 + t243 * t3 + t4 * t54) - t114 * t59 + t120 * t100 + t128 * t27 + t92 * t111 - pkin(2) * t101 + t9 * t82 + t83 * t26 + t10 * t79 + t6 * t80 + t11 * t81 + t54 * t48 + t25 * t58 + (t136 / 0.2e1 + t137 / 0.2e1 - t71 / 0.2e1 + t72 / 0.2e1 + t73 / 0.2e1 + t70 / 0.2e1 - t56 * mrSges(5,3) + t212 * t76 + t217 * t75 + t281 * t42 + (mrSges(4,2) * t228 + t156 * t235) * t232 + 0.2e1 * (t198 * t67 + t219) * m(5) + 0.2e1 * (t198 * t84 + t219) * m(4) + (pkin(8) * t168 + t215 * qJD(3) + t213 * t229 - t295) * qJD(3) + t163) * t152 - m(7) * (t7 * t77 + t78 * t8) - m(6) * (-t14 * t77 + t15 * t78) + t45 * t47 + t46 * t49 + m(7) * (t1 * t45 + t16 * t25 + t2 * t46 + t5 * t83 + t6 * t8 + t7 * t9) + (t56 * mrSges(5,2) + t41 * mrSges(4,3) - t32 * mrSges(5,1) + t17 * t259 + t18 * t260 + t181 * t271 + t185 * t270 + t28 * t194 + t5 * t192 + (t4 * t151 - t3 * t154) * mrSges(6,3) + (-t1 * t154 - t2 * t151) * mrSges(7,2) + (-mrSges(4,1) * t228 + (-m(6) * t52 - m(7) * t16 + t168 - t58 - t59) * t156) * t232 + (t38 * t261 + t186 * t267 + t182 * t268 - t52 * t193 - t16 * t191 + (t14 * t154 + t15 * t151) * mrSges(6,3) + (t8 * t151 - t7 * t154) * mrSges(7,2) + t289 * t266 + (t151 * t279 - t154 * t292) * t139 / 0.2e1 + t278 * t260) * qJD(5) + (((-t151 * t217 + t154 * t212 - t213) * t155 + (-0.3e1 / 0.2e1 * Ifges(4,2) - 0.3e1 / 0.2e1 * Ifges(5,3) + 0.3e1 / 0.2e1 * Ifges(5,2) + 0.3e1 / 0.2e1 * Ifges(4,1) + t214) * t152) * qJD(2) + t216 * qJD(3) - t275) * qJD(3) - t288 * t75 / 0.2e1 + (qJD(5) * t35 + t291) * t262 + (0.2e1 * t32 * t273 + m(4) * t41 + (m(4) * t84 + m(5) * t67 - t235) * qJD(3)) * pkin(8)) * t155 + 0.2e1 * (t86 * t273 + (-t245 / 0.2e1 - t118 / 0.2e1) * m(4)) * t209 + t254 * t77 - t255 * t78 - t237 * t209 + m(5) * (t120 * t56 + t86 * t92); ((-m(6) * t177 + m(7) * t195 + t164) * t157 + t297) * qJD(5) + t119 * t26 - t113 * t111 - t22 * t82 - t24 * t79 - t21 * t80 - t23 * t81 - t64 * t59 + t28 * t193 + t196 * mrSges(6,3) + t197 * mrSges(7,2) + t5 * t191 + ((t285 + (-pkin(3) * mrSges(5,1) + t151 * t212 + t154 * t217 + t216) * qJD(3) + Ifges(5,6) * t294 + t275) * t155 + ((Ifges(4,4) / 0.2e1 + Ifges(5,6) / 0.2e1) * t229 + (-qJ(4) * mrSges(5,1) + t215) * qJD(3) + (Ifges(4,2) / 0.2e1 + t287 - Ifges(5,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t227 + t295) * t152) * qJD(2) - t41 * mrSges(4,2) - t32 * mrSges(5,3) + qJ(4) * t27 - m(6) * (t14 * t23 + t15 * t24 + t52 * t64) - m(7) * (t16 * t30 + t21 * t8 + t22 * t7) + t289 * t75 / 0.2e1 + (-pkin(3) * t42 - qJ(4) * t32 - t113 * t86 + t290 * t68 - t67 * t85) * m(5) + (t89 - t30) * t58 + t165 * t157 + t291 * t259 + t293 * t42 + t182 * t270 + t186 * t271 + t17 * t261 + t18 * t262 + t234 * t84 + t235 * t85 + m(6) * (qJ(4) * t28 + qJD(4) * t52 + t238 * t4 + t239 * t3) + m(7) * (t1 * t239 + t119 * t5 + t16 * t89 - t2 * t238) + t244 * qJD(4); t218 * qJD(3) + t164 * qJD(5) + (mrSges(5,1) * t225 + (t111 + t164) * t152) * qJD(2) + t165 + (-qJD(3) * t16 + t139 * t195 - t197) * m(7) + (-qJD(3) * t52 - t139 * t177 - t196) * m(6) + (qJD(3) * t68 + t229 * t86 + t42) * m(5); t163 + t136 + t137 - t16 * (mrSges(7,1) * t109 + mrSges(7,3) * t108) - t52 * (mrSges(6,1) * t109 - mrSges(6,2) * t108) + qJD(6) * t80 - t57 * t58 - t71 + t72 + t73 + t70 + qJ(6) * t47 - pkin(5) * t49 + (t108 * t7 + t109 * t8) * mrSges(7,2) + t38 * t265 + (Ifges(7,3) * t109 - t248) * t268 + (t252 + t254) * t15 + (-t253 - t255) * t14 + (-t108 * t292 - t109 * t279) * t264 + (-pkin(5) * t2 + qJ(6) * t1 - t15 * t7 - t16 * t57 + t277 * t8) * m(7) + (-Ifges(6,2) * t109 - t105 + t278) * t267 + (-t108 * t296 + t104 - t251 + t35) * t266; t109 * t58 - t139 * t80 + 0.2e1 * (t2 / 0.2e1 + t16 * t265 + t8 * t264) * m(7) + t49;];
tauc  = t19(:);
