% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRR13_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR13_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR13_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR13_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR13_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:30:24
% EndTime: 2024-09-27 17:30:29
% DurationCPUTime: 3.92s
% Computational Cost: add. (6646->333), mult. (12050->487), div. (0->0), fcn. (7659->10), ass. (0->209)
t189 = cos(qJ(3));
t190 = cos(qJ(2));
t185 = sin(qJ(3));
t186 = sin(qJ(2));
t246 = t185 * t186;
t202 = t189 * t190 - t246;
t260 = pkin(1) * qJD(1);
t147 = t202 * t260;
t241 = qJD(3) * t189;
t300 = -pkin(2) * t241 + t147;
t180 = qJD(1) + qJD(2);
t160 = t180 * pkin(2) + t190 * t260;
t232 = t186 * t260;
t134 = t189 * t160 - t185 * t232;
t135 = t160 * t185 + t189 * t232;
t188 = cos(qJ(4));
t182 = cos(pkin(5));
t247 = t182 * t188;
t224 = qJD(4) * t247;
t184 = sin(qJ(4));
t248 = t182 * t184;
t299 = pkin(3) * t224 - t134 * t188 + t135 * t248;
t245 = t186 * t189;
t203 = -t185 * t190 - t245;
t146 = t203 * t260;
t298 = -t146 * t247 + t147 * t184 + (-t184 * t189 - t185 * t247) * qJD(3) * pkin(2);
t174 = pkin(2) * t189 + pkin(3);
t297 = -t146 * t248 + t174 * t224 - t300 * t188;
t242 = qJD(3) * t185;
t231 = pkin(2) * t242;
t217 = t182 * t231;
t181 = sin(pkin(5));
t177 = t181 * pkin(9);
t162 = pkin(2) * t185 + t177;
t270 = pkin(10) * t181;
t220 = -t162 - t270;
t296 = (t220 * qJD(4) - t217) * t184 + t297;
t158 = t174 * t248;
t295 = (t220 * t188 - t158) * qJD(4) + t298;
t229 = t181 * (-pkin(9) - pkin(10));
t216 = t184 * t229;
t294 = qJD(4) * t216 + t299;
t171 = pkin(3) * t248;
t69 = -t134 * t184 - t135 * t247;
t293 = (t188 * t229 - t171) * qJD(4) - t69;
t176 = qJD(3) + t180;
t183 = sin(qJ(5));
t187 = cos(qJ(5));
t138 = (-t183 * t184 + t187 * t188) * t181;
t278 = qJD(4) + qJD(5);
t96 = t278 * t138;
t64 = t176 * t96;
t292 = Ifges(6,5) * t64;
t139 = (t183 * t188 + t184 * t187) * t181;
t97 = t278 * t139;
t65 = t176 * t97;
t291 = Ifges(6,6) * t65;
t210 = t146 + t231;
t251 = t176 * t181;
t110 = pkin(9) * t251 + t135;
t211 = pkin(10) * t251 + t110;
t201 = t211 * t184;
t116 = pkin(3) * t176 + t134;
t193 = (t202 * qJD(2) - t186 * t242) * pkin(1);
t88 = qJD(1) * t193 + t160 * t241;
t192 = (t203 * qJD(2) - t186 * t241) * pkin(1);
t89 = qJD(1) * t192 - t160 * t242;
t234 = t116 * t224 + t188 * t88 + t89 * t248;
t12 = -qJD(4) * t201 + t234;
t219 = -t184 * t88 + t89 * t247;
t228 = t116 * t248;
t53 = t211 * t188 + t228;
t13 = -t53 * qJD(4) + t219;
t256 = t183 * t53;
t161 = t182 * t176 + qJD(4);
t103 = t116 * t247;
t52 = t103 - t201;
t47 = pkin(4) * t161 + t52;
t14 = t187 * t47 - t256;
t4 = t14 * qJD(5) + t12 * t187 + t13 * t183;
t255 = t187 * t53;
t15 = t183 * t47 + t255;
t5 = -t15 * qJD(5) - t12 * t183 + t13 * t187;
t290 = t5 * mrSges(6,1) - t4 * mrSges(6,2);
t239 = qJD(4) * t184;
t20 = -t110 * t239 + t234;
t60 = t110 * t188 + t228;
t21 = -t60 * qJD(4) + t219;
t289 = t21 * mrSges(5,1) - t20 * mrSges(5,2);
t159 = t174 * t247;
t178 = t182 * pkin(4);
t102 = t220 * t184 + t159 + t178;
t126 = t188 * t162 + t158;
t249 = t181 * t188;
t169 = pkin(10) * t249;
t112 = t169 + t126;
t54 = t102 * t187 - t112 * t183;
t287 = t54 * qJD(5) + t295 * t183 + t296 * t187;
t55 = t102 * t183 + t112 * t187;
t286 = -t55 * qJD(5) - t296 * t183 + t295 * t187;
t172 = pkin(3) * t247;
t122 = t172 + t178 + t216;
t150 = pkin(9) * t249 + t171;
t136 = t169 + t150;
t74 = t122 * t183 + t136 * t187;
t285 = -t74 * qJD(5) - t294 * t183 + t293 * t187;
t73 = t122 * t187 - t136 * t183;
t284 = t73 * qJD(5) + t293 * t183 + t294 * t187;
t283 = (-qJD(4) * t162 - t217) * t184 + t297;
t282 = -t126 * qJD(4) + t298;
t281 = -t150 * qJD(4) - t69;
t225 = t181 * t239;
t280 = -pkin(9) * t225 + t299;
t175 = pkin(1) * t190 + pkin(2);
t212 = -pkin(1) * t246 + t189 * t175;
t144 = pkin(3) + t212;
t132 = t144 * t248;
t243 = pkin(1) * t245 + t185 * t175;
t137 = t177 + t243;
t79 = t188 * t137 + t132;
t279 = mrSges(3,1) * t186 + mrSges(3,2) * t190;
t179 = t181 ^ 2;
t120 = t176 * t139;
t276 = t120 / 0.2e1;
t274 = t184 / 0.2e1;
t272 = pkin(4) * t184;
t271 = pkin(4) * t188;
t269 = t14 * t96;
t268 = -t291 + t292;
t266 = mrSges(4,1) * t176;
t119 = t176 * t138;
t264 = mrSges(6,3) * t119;
t263 = Ifges(5,4) * t184;
t262 = Ifges(5,4) * t188;
t261 = Ifges(6,4) * t120;
t259 = t120 * mrSges(6,3);
t258 = t179 * t89;
t87 = (-t176 * t271 - t116) * t181;
t257 = t181 * t87;
t254 = t188 * Ifges(5,2);
t253 = mrSges(5,3) * qJD(4);
t252 = t116 * t179;
t250 = t181 * t184;
t240 = qJD(4) * t181;
t238 = qJD(5) * t183;
t237 = qJD(5) * t187;
t104 = t175 * t241 + t193;
t105 = -t175 * t242 + t192;
t233 = t188 * t104 + t105 * t248 + t144 * t224;
t227 = t176 * t250;
t226 = t176 * t249;
t223 = t176 * t239;
t222 = t240 / 0.2e1;
t221 = -t137 - t270;
t218 = -t104 * t184 + t105 * t247;
t165 = pkin(4) * t225;
t214 = -t225 / 0.2e1;
t213 = t188 * t222;
t209 = t221 * t184;
t208 = mrSges(5,1) * t184 + mrSges(5,2) * t188;
t133 = t144 * t247;
t68 = t133 + t178 + t209;
t75 = t169 + t79;
t31 = -t183 * t75 + t187 * t68;
t32 = t183 * t68 + t187 * t75;
t59 = -t110 * t184 + t103;
t207 = t60 * t184 + t59 * t188;
t140 = t208 * t240;
t148 = (-mrSges(5,1) * t188 + mrSges(5,2) * t184) * t181;
t206 = -t116 * t140 - t89 * t148;
t200 = t161 * (Ifges(5,5) * t188 - Ifges(5,6) * t184);
t199 = t184 * (Ifges(5,1) * t188 - t263);
t198 = (t254 + t263) * t181;
t117 = Ifges(6,4) * t119;
t157 = qJD(5) + t161;
t57 = Ifges(6,2) * t119 + Ifges(6,6) * t157 + t261;
t58 = Ifges(6,1) * t120 + Ifges(6,5) * t157 + t117;
t194 = t14 * t264 + t57 * t276 - t87 * (mrSges(6,1) * t120 + mrSges(6,2) * t119) + t268 - t120 * (Ifges(6,1) * t119 - t261) / 0.2e1 - t157 * (Ifges(6,5) * t119 - Ifges(6,6) * t120) / 0.2e1 - (-Ifges(6,2) * t120 + t117 + t58) * t119 / 0.2e1 + t290;
t100 = Ifges(5,6) * t161 + t176 * t198;
t156 = Ifges(5,4) * t226;
t101 = Ifges(5,1) * t227 + Ifges(5,5) * t161 + t156;
t152 = Ifges(5,5) * qJD(4) * t226;
t61 = (pkin(4) * t223 - t89) * t181;
t191 = t200 * t222 + t89 * mrSges(4,1) + t96 * t58 / 0.2e1 - t97 * t57 / 0.2e1 + t87 * (mrSges(6,1) * t97 + mrSges(6,2) * t96) + t119 * (Ifges(6,4) * t96 - Ifges(6,2) * t97) / 0.2e1 + t157 * (Ifges(6,5) * t96 - Ifges(6,6) * t97) / 0.2e1 - t15 * t97 * mrSges(6,3) + (Ifges(6,1) * t96 - Ifges(6,4) * t97) * t276 + t101 * t213 + t100 * t214 + (-Ifges(5,6) * t181 * t223 + t152 + t268) * t182 / 0.2e1 + (t20 * t249 - t21 * t250) * mrSges(5,3) + (-t291 / 0.2e1 + t292 / 0.2e1 + t289 + t290) * t182 + (t61 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,1) * t64 - Ifges(6,4) * t65) * t139 + (-t61 * mrSges(6,1) + t4 * mrSges(6,3) + Ifges(6,4) * t64 - Ifges(6,2) * t65) * t138 + ((t188 * (-Ifges(5,2) * t184 + t262) + t199) * qJD(4) * t179 + (Ifges(5,5) * t182 + (t184 * Ifges(5,1) + t262) * t181) * t213 + (Ifges(5,6) * t182 + t198) * t214) * t176;
t155 = (-pkin(3) - t271) * t181;
t149 = -pkin(9) * t250 + t172;
t145 = (-t174 - t271) * t181;
t141 = t181 * t231 + t165;
t129 = t176 * t148;
t128 = -mrSges(5,2) * t161 + mrSges(5,3) * t226;
t127 = mrSges(5,1) * t161 - mrSges(5,3) * t227;
t125 = -t162 * t184 + t159;
t121 = t176 * t140;
t118 = (-t144 - t271) * t181;
t93 = mrSges(6,1) * t157 - t259;
t92 = -mrSges(6,2) * t157 + t264;
t81 = -t105 * t181 + t165;
t78 = -t137 * t184 + t133;
t67 = -mrSges(6,1) * t119 + mrSges(6,2) * t120;
t30 = mrSges(6,1) * t65 + mrSges(6,2) * t64;
t28 = -t79 * qJD(4) + t218;
t27 = -t137 * t239 + t233;
t25 = (t221 * t188 - t132) * qJD(4) + t218;
t24 = qJD(4) * t209 + t233;
t19 = t187 * t52 - t256;
t18 = -t183 * t52 - t255;
t7 = -t32 * qJD(5) - t183 * t24 + t187 * t25;
t6 = t31 * qJD(5) + t183 * t25 + t187 * t24;
t1 = [t105 * t266 + t81 * t67 + t6 * t92 + t7 * t93 + t118 * t30 + t28 * t127 + t27 * t128 + t191 + m(4) * (t135 * t104 + t134 * t105 + t89 * t212 + t88 * t243) + m(6) * (t118 * t61 + t14 * t7 + t15 * t6 + t31 * t5 + t32 * t4 + t81 * t87) + (-t105 * t129 - t144 * t121 + ((-t184 * t79 - t188 * t78) * t176 - t207) * t253 + t206) * t181 + (-t31 * t64 - t32 * t65 - t269) * mrSges(6,3) + (-t104 * t176 - t88) * mrSges(4,2) + m(5) * (t20 * t79 + t21 * t78 + t27 * t60 + t28 * t59 + (t105 * t116 + t144 * t89) * t179) + t279 * pkin(1) * qJD(2) * (-qJD(1) - t180); -t88 * mrSges(4,2) + t141 * t67 + t145 * t30 + (-t174 * t121 + t146 * t67 + t210 * t129 + ((-t125 * t188 - t126 * t184) * t176 - t207) * t253 + t206) * t181 + t191 + (-t210 * mrSges(4,1) + t300 * mrSges(4,2)) * t176 + t286 * t93 + t287 * t92 + t283 * t128 + t282 * t127 + (-t54 * t64 - t55 * t65 - t269) * mrSges(6,3) + t279 * t260 * (-qJD(2) + t180) + (t286 * t14 + t141 * t87 + t145 * t61 + t146 * t257 + t287 * t15 + t4 * t55 + t5 * t54) * m(6) + (-t134 * t146 - t135 * t147 + (-t134 * t242 + t135 * t241 + t185 * t88 + t189 * t89) * pkin(2)) * m(4) + (t125 * t21 + t126 * t20 + t174 * t258 - t210 * t252 + t282 * t59 + t283 * t60) * m(5); t135 * t266 + (-pkin(3) * t121 + (-t129 - t67) * t135 + (t67 * t272 + ((-t149 * t188 - t150 * t184) * t176 - t207) * mrSges(5,3)) * qJD(4) + t206) * t181 + t155 * t30 + t191 + t280 * t128 + t281 * t127 + (-t64 * t73 - t65 * t74 - t269) * mrSges(6,3) + (t134 * t176 - t88) * mrSges(4,2) + t285 * t93 + t284 * t92 + (-t135 * t257 + t285 * t14 + t284 * t15 + t155 * t61 + t87 * t165 + t4 * t74 + t5 * t73) * m(6) + (pkin(3) * t258 + t135 * t252 + t149 * t21 + t150 * t20 + t280 * t60 + t281 * t59) * m(5); -t19 * t92 - t18 * t93 + t60 * t127 - t59 * t128 + t15 * t259 + t194 + (-t200 / 0.2e1 + t100 * t274 - Ifges(5,6) * t239 + (t116 * t208 + (t254 * t274 - t199 / 0.2e1) * t176) * t181 + (-m(6) * t87 - t67) * t272 + t207 * mrSges(5,3) - (t156 + t101) * t188 / 0.2e1) * t251 - m(6) * (t14 * t18 + t15 * t19) + (m(6) * (-t14 * t238 + t15 * t237 + t183 * t4 + t187 * t5) + t92 * t237 - t93 * t238 + (-t183 * t65 - t187 * t64) * mrSges(6,3)) * pkin(4) + t152 + t289; t194 - t14 * t92 + (t93 + t259) * t15;];
tauc = t1(:);
