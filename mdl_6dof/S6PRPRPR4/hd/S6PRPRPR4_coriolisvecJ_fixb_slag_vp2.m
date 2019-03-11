% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-03-08 19:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRPR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:39:05
% EndTime: 2019-03-08 19:39:18
% DurationCPUTime: 8.00s
% Computational Cost: add. (6287->467), mult. (16696->664), div. (0->0), fcn. (13176->12), ass. (0->221)
t190 = cos(pkin(11));
t268 = cos(qJ(4));
t223 = t268 * t190;
t187 = sin(pkin(11));
t193 = sin(qJ(4));
t235 = t193 * t187;
t203 = t223 - t235;
t188 = sin(pkin(6));
t196 = cos(qJ(2));
t237 = t188 * t196;
t200 = t203 * t237;
t263 = pkin(8) + qJ(3);
t173 = t263 * t187;
t175 = t263 * t190;
t287 = -t268 * t173 - t193 * t175;
t310 = qJD(1) * t200 - t203 * qJD(3) - qJD(4) * t287;
t160 = t203 * qJD(4);
t168 = t187 * t268 + t193 * t190;
t161 = t168 * qJD(4);
t194 = sin(qJ(2));
t229 = qJD(1) * t188;
t221 = t194 * t229;
t309 = pkin(4) * t161 - qJ(5) * t160 - qJD(5) * t168 - t221;
t186 = sin(pkin(12));
t189 = cos(pkin(12));
t301 = t310 * t186 + t309 * t189;
t300 = t309 * t186 - t310 * t189;
t266 = pkin(9) * t189;
t308 = pkin(5) * t161 - t160 * t266 + t301;
t240 = t160 * t186;
t307 = pkin(9) * t240 - t300;
t306 = Ifges(5,1) / 0.2e1;
t156 = t203 * qJD(2);
t305 = t156 / 0.2e1;
t192 = sin(qJ(6));
t195 = cos(qJ(6));
t181 = -pkin(3) * t190 - pkin(2);
t123 = -pkin(4) * t203 - qJ(5) * t168 + t181;
t132 = -t193 * t173 + t175 * t268;
t69 = t189 * t123 - t132 * t186;
t50 = -pkin(5) * t203 - t168 * t266 + t69;
t239 = t168 * t186;
t70 = t186 * t123 + t189 * t132;
t54 = -pkin(9) * t239 + t70;
t17 = t192 * t50 + t195 * t54;
t304 = -qJD(6) * t17 + t307 * t192 + t308 * t195;
t16 = -t192 * t54 + t195 * t50;
t303 = qJD(6) * t16 + t308 * t192 - t307 * t195;
t302 = Ifges(5,2) + Ifges(6,3);
t100 = qJD(3) * t168 + qJD(4) * t132;
t201 = t168 * t237;
t133 = qJD(1) * t201;
t298 = t100 - t133;
t157 = t168 * qJD(2);
t137 = qJD(4) * t186 + t157 * t189;
t217 = t189 * qJD(4) - t157 * t186;
t297 = -t137 * t192 + t195 * t217;
t84 = t137 * t195 + t192 * t217;
t149 = qJD(2) * t160;
t205 = t186 * t192 - t189 * t195;
t44 = qJD(6) * t297 - t149 * t205;
t281 = t44 / 0.2e1;
t167 = t186 * t195 + t189 * t192;
t45 = -qJD(6) * t84 - t149 * t167;
t280 = t45 / 0.2e1;
t150 = qJD(2) * t161;
t273 = t150 / 0.2e1;
t296 = Ifges(5,4) * t305;
t295 = -t156 / 0.2e1;
t294 = t157 * t306;
t262 = pkin(9) + qJ(5);
t172 = t262 * t186;
t174 = t262 * t189;
t131 = -t172 * t192 + t174 * t195;
t121 = pkin(4) * t157 - qJ(5) * t156;
t171 = qJD(2) * qJ(3) + t221;
t191 = cos(pkin(6));
t228 = qJD(1) * t191;
t140 = t190 * t171 + t187 * t228;
t256 = pkin(8) * qJD(2);
t128 = t190 * t256 + t140;
t125 = t193 * t128;
t177 = t190 * t228;
t127 = t177 + (-t171 - t256) * t187;
t224 = t268 * t127;
t78 = -t125 + t224;
t42 = t189 * t121 - t186 * t78;
t22 = pkin(5) * t157 - t156 * t266 + t42;
t241 = t156 * t186;
t43 = t186 * t121 + t189 * t78;
t32 = -pkin(9) * t241 + t43;
t293 = -qJD(5) * t167 - qJD(6) * t131 + t192 * t32 - t195 * t22;
t14 = -t45 * mrSges(7,1) + t44 * mrSges(7,2);
t242 = t149 * t189;
t243 = t149 * t186;
t96 = mrSges(6,1) * t243 + mrSges(6,2) * t242;
t292 = t14 + t96;
t129 = -t172 * t195 - t174 * t192;
t291 = -qJD(5) * t205 + qJD(6) * t129 - t192 * t22 - t195 * t32;
t230 = t187 ^ 2 + t190 ^ 2;
t290 = mrSges(4,3) * t230;
t105 = t205 * t156;
t158 = t205 * qJD(6);
t232 = -t158 + t105;
t104 = t167 * t156;
t159 = t167 * qJD(6);
t231 = -t159 + t104;
t227 = qJD(1) * t196;
t220 = t188 * t227;
t209 = qJD(3) - t220;
t164 = (qJD(3) + t220) * qJD(2);
t233 = qJD(4) * t224 + t164 * t223;
t41 = -t164 * t235 + (qJD(5) - t125) * qJD(4) + t233;
t238 = t188 * t194;
t219 = qJD(2) * t238;
t216 = qJD(1) * t219;
t75 = pkin(4) * t150 - qJ(5) * t149 - qJD(5) * t157 + t216;
t18 = -t186 * t41 + t189 * t75;
t13 = pkin(5) * t150 - pkin(9) * t242 + t18;
t19 = t186 * t75 + t189 * t41;
t15 = -pkin(9) * t243 + t19;
t79 = t193 * t127 + t128 * t268;
t72 = qJD(4) * qJ(5) + t79;
t151 = qJD(2) * t181 + t209;
t91 = -pkin(4) * t156 - qJ(5) * t157 + t151;
t29 = -t186 * t72 + t189 * t91;
t20 = -pkin(5) * t156 - pkin(9) * t137 + t29;
t30 = t186 * t91 + t189 * t72;
t21 = pkin(9) * t217 + t30;
t5 = -t192 * t21 + t195 * t20;
t1 = qJD(6) * t5 + t13 * t192 + t15 * t195;
t6 = t192 * t20 + t195 * t21;
t2 = -qJD(6) * t6 + t13 * t195 - t15 * t192;
t286 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t44 + Ifges(7,6) * t45;
t247 = t186 * Ifges(6,2);
t259 = Ifges(6,4) * t189;
t211 = -t247 + t259;
t260 = Ifges(6,4) * t186;
t212 = t189 * Ifges(6,1) - t260;
t213 = mrSges(6,1) * t186 + mrSges(6,2) * t189;
t269 = t189 / 0.2e1;
t270 = -t186 / 0.2e1;
t68 = -qJD(4) * pkin(4) + qJD(5) - t78;
t285 = t211 * t217 / 0.2e1 + t212 * t137 / 0.2e1 + t68 * t213 + t151 * mrSges(5,2) + (-t186 * t30 - t189 * t29) * mrSges(6,3) + t296 + Ifges(5,5) * qJD(4) + t294 + (t137 * Ifges(6,4) + Ifges(6,2) * t217 - Ifges(6,6) * t156) * t270 + (t137 * Ifges(6,1) + Ifges(6,4) * t217 - Ifges(6,5) * t156) * t269;
t284 = m(4) / 0.2e1;
t283 = Ifges(7,4) * t281 + Ifges(7,2) * t280 + Ifges(7,6) * t273;
t282 = Ifges(7,1) * t281 + Ifges(7,4) * t280 + Ifges(7,5) * t273;
t279 = -t297 / 0.2e1;
t278 = t297 / 0.2e1;
t277 = -t84 / 0.2e1;
t276 = t84 / 0.2e1;
t153 = qJD(6) - t156;
t272 = -t153 / 0.2e1;
t271 = t153 / 0.2e1;
t267 = Ifges(7,4) * t84;
t261 = Ifges(5,4) * t157;
t258 = Ifges(6,5) * t189;
t257 = Ifges(6,6) * t186;
t154 = -t187 * t238 + t190 * t191;
t155 = t187 * t191 + t190 * t238;
t204 = t154 * t268 - t193 * t155;
t47 = qJD(4) * t79 + t164 * t168;
t255 = t204 * t47;
t254 = t287 * t47;
t251 = t150 * Ifges(6,5);
t250 = t150 * Ifges(6,6);
t246 = qJD(4) * mrSges(5,1) + mrSges(6,1) * t217 - mrSges(6,2) * t137 - mrSges(5,3) * t157;
t214 = -mrSges(4,1) * t190 + mrSges(4,2) * t187;
t234 = -mrSges(5,1) * t156 + mrSges(5,2) * t157 + qJD(2) * t214;
t34 = -mrSges(7,1) * t297 + mrSges(7,2) * t84;
t226 = -t34 + t246;
t225 = Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t222 = t188 ^ 2 * t227;
t110 = t150 * mrSges(5,1) + t149 * mrSges(5,2);
t218 = t230 * t164;
t210 = -t257 + t258;
t208 = t18 * t189 + t186 * t19;
t207 = t186 * t29 - t189 * t30;
t112 = t193 * t154 + t155 * t268;
t92 = -t112 * t186 - t189 * t237;
t93 = t112 * t189 - t186 * t237;
t36 = -t192 * t93 + t195 * t92;
t37 = t192 * t92 + t195 * t93;
t206 = -(-t171 * t187 + t177) * t187 + t140 * t190;
t202 = t206 * t196;
t198 = t151 * mrSges(5,1) + t29 * mrSges(6,1) + t5 * mrSges(7,1) - Ifges(5,6) * qJD(4) - t261 / 0.2e1 + t153 * Ifges(7,3) + t84 * Ifges(7,5) + t297 * Ifges(7,6) + t137 * Ifges(6,5) + t217 * Ifges(6,6) - t30 * mrSges(6,2) - t6 * mrSges(7,2) + t302 * t295;
t197 = qJD(2) ^ 2;
t180 = -pkin(5) * t189 - pkin(4);
t170 = -qJD(2) * pkin(2) + t209;
t146 = Ifges(7,3) * t150;
t143 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t156;
t116 = t205 * t168;
t115 = t167 * t168;
t107 = mrSges(6,1) * t150 - mrSges(6,3) * t242;
t106 = -mrSges(6,2) * t150 - mrSges(6,3) * t243;
t103 = pkin(5) * t239 - t287;
t102 = -mrSges(6,1) * t156 - mrSges(6,3) * t137;
t101 = mrSges(6,2) * t156 + mrSges(6,3) * t217;
t81 = Ifges(7,4) * t297;
t80 = pkin(5) * t240 + t100;
t77 = t149 * t212 + t251;
t76 = t149 * t211 + t250;
t74 = qJD(2) * t201 + qJD(4) * t112;
t73 = qJD(2) * t200 + qJD(4) * t204;
t61 = mrSges(7,1) * t153 - mrSges(7,3) * t84;
t60 = -mrSges(7,2) * t153 + mrSges(7,3) * t297;
t59 = t186 * t219 + t189 * t73;
t58 = -t186 * t73 + t189 * t219;
t57 = t158 * t168 - t160 * t167;
t56 = -t159 * t168 - t160 * t205;
t55 = pkin(5) * t241 + t79;
t51 = -pkin(5) * t217 + t68;
t46 = (-qJD(4) * t128 - t164 * t187) * t193 + t233;
t31 = pkin(5) * t243 + t47;
t28 = -mrSges(7,2) * t150 + mrSges(7,3) * t45;
t27 = mrSges(7,1) * t150 - mrSges(7,3) * t44;
t26 = Ifges(7,1) * t84 + Ifges(7,5) * t153 + t81;
t25 = Ifges(7,2) * t297 + Ifges(7,6) * t153 + t267;
t8 = -qJD(6) * t37 - t192 * t59 + t195 * t58;
t7 = qJD(6) * t36 + t192 * t58 + t195 * t59;
t3 = [-t112 * t150 * mrSges(5,3) + t59 * t101 + t58 * t102 + t93 * t106 + t92 * t107 + t73 * t143 + t36 * t27 + t37 * t28 + t7 * t60 + t8 * t61 - t226 * t74 - (mrSges(5,3) * t149 + t292) * t204 + ((-mrSges(3,1) * t197 + qJD(2) * t234) * t194 + (-t110 + (-mrSges(3,2) + t290) * t197) * t196) * t188 + m(7) * (t1 * t37 + t2 * t36 - t204 * t31 + t5 * t8 + t51 * t74 + t6 * t7) + m(6) * (t18 * t92 + t19 * t93 + t29 * t58 + t30 * t59 + t68 * t74 - t255) + m(5) * (t112 * t46 + t73 * t79 - t74 * t78 - t255) + m(4) * (-t154 * t187 + t155 * t190) * t164 + 0.2e1 * (t188 * t202 * t284 + (m(5) * (t151 * t188 - t222) / 0.2e1 + (t170 * t188 - t222) * t284) * t194) * qJD(2); (t132 * t46 - t151 * t221 + t181 * t216 - t298 * t78 - t310 * t79 - t254) * m(5) - t310 * t143 + ((-mrSges(5,1) * t203 + mrSges(5,2) * t168 + t214) * qJD(2) - t234) * t221 + (t210 * t295 + t285 + t294) * t160 + (t149 * t203 - t150 * t168 + t160 * t305 - t157 * t161 / 0.2e1) * Ifges(5,4) + t303 * t60 + (t1 * t17 + t103 * t31 + t16 * t2 + t303 * t6 + (-t133 + t80) * t51 + t304 * t5) * m(7) + t304 * t61 - (t18 * mrSges(6,1) - t19 * mrSges(6,2) + t146 / 0.2e1 + t210 * t149 + (Ifges(7,3) / 0.2e1 + t302) * t150 + t286) * t203 + t300 * t101 + (t18 * t69 + t19 * t70 + t301 * t29 + t298 * t68 + t300 * t30 - t254) * m(6) + t301 * t102 + (qJD(2) * t209 * t230 + t218) * mrSges(4,3) + (-t156 * t225 + t198) * t161 + (-t132 * t150 - t149 * t287 - t160 * t78 - t161 * t79 + t168 * t47 + t203 * t46) * mrSges(5,3) - t287 * t96 + (-t1 * t115 + t116 * t2 - t5 * t56 + t57 * t6) * mrSges(7,3) + (-Ifges(7,5) * t116 - Ifges(7,6) * t115) * t273 + (-Ifges(7,4) * t116 - Ifges(7,2) * t115) * t280 + (-Ifges(7,1) * t116 - Ifges(7,4) * t115) * t281 + t31 * (mrSges(7,1) * t115 - mrSges(7,2) * t116) + (-(t170 * t194 + t202) * t229 - pkin(2) * t216 + qJ(3) * t218 + qJD(3) * t206) * m(4) - t115 * t283 + (Ifges(7,5) * t56 + Ifges(7,6) * t57) * t271 + (Ifges(7,1) * t56 + Ifges(7,4) * t57) * t276 + (Ifges(7,4) * t56 + Ifges(7,2) * t57) * t278 - t116 * t282 + t181 * t110 + t70 * t106 + t69 * t107 + t103 * t14 + t80 * t34 + t51 * (-mrSges(7,1) * t57 + mrSges(7,2) * t56) + t57 * t25 / 0.2e1 + t56 * t26 / 0.2e1 + t16 * t27 + t17 * t28 + t226 * t133 - t246 * t100 + (t210 * t273 + t47 * t213 + t76 * t270 + t77 * t269 - t208 * mrSges(6,3) + (Ifges(5,1) + t189 ^ 2 * Ifges(6,1) / 0.2e1 + (-t259 + t247 / 0.2e1) * t186) * t149) * t168; t186 * t106 + t189 * t107 - t205 * t27 + t167 * t28 + t231 * t61 + t232 * t60 + (m(4) + m(5)) * t216 - t197 * t290 + t226 * t157 - (t101 * t189 - t102 * t186 + t143) * t156 - m(5) * (t156 * t79 - t157 * t78) - m(4) * t206 * qJD(2) + t110 + (t1 * t167 - t157 * t51 - t2 * t205 + t231 * t5 + t232 * t6) * m(7) + (t156 * t207 - t157 * t68 + t208) * m(6); (-t159 / 0.2e1 + t104 / 0.2e1) * t25 + (-Ifges(7,5) * t105 - Ifges(7,6) * t104) * t272 + (-Ifges(7,1) * t105 - Ifges(7,4) * t104) * t277 + (-Ifges(7,4) * t105 - Ifges(7,2) * t104) * t279 + (-t158 / 0.2e1 + t105 / 0.2e1) * t26 + (Ifges(7,5) * t167 - Ifges(7,6) * t205) * t273 + (Ifges(7,4) * t167 - Ifges(7,2) * t205) * t280 + (Ifges(7,1) * t167 - Ifges(7,4) * t205) * t281 + t31 * (mrSges(7,1) * t205 + mrSges(7,2) * t167) + (-t1 * t205 - t167 * t2 + t231 * t6 - t232 * t5) * mrSges(7,3) - t205 * t283 + (t1 * t131 + t129 * t2 + t180 * t31 + t291 * t6 + t293 * t5 - t51 * t55) * m(7) + t293 * t61 + t291 * t60 - (t296 - t78 * mrSges(5,3) - (t258 / 0.2e1 - t257 / 0.2e1) * t156 + (t306 - t225) * t157 + t285) * t156 + (-pkin(4) * t47 - t207 * qJD(5) + (-t18 * t186 + t189 * t19) * qJ(5) - t29 * t42 - t30 * t43 - t68 * t79) * m(6) + (-Ifges(7,5) * t158 - Ifges(7,6) * t159) * t271 + (-Ifges(7,1) * t158 - Ifges(7,4) * t159) * t276 + (-Ifges(7,4) * t158 - Ifges(7,2) * t159) * t278 + t167 * t282 + t180 * t14 - Ifges(5,6) * t150 - t78 * t143 + t129 * t27 + t131 * t28 - pkin(4) * t96 - t43 * t101 - t42 * t102 - t55 * t34 - t46 * mrSges(5,2) - t47 * mrSges(5,1) + (-mrSges(7,1) * t231 + mrSges(7,2) * t232) * t51 + t246 * t79 + (t250 / 0.2e1 - t47 * mrSges(6,1) + t76 / 0.2e1 + qJD(5) * t101 + qJ(5) * t106 + t19 * mrSges(6,3)) * t189 + (t251 / 0.2e1 + t47 * mrSges(6,2) + t77 / 0.2e1 - qJD(5) * t102 - qJ(5) * t107 - t18 * mrSges(6,3)) * t186 + (t79 * mrSges(5,3) + t261 / 0.2e1 - t198) * t157 + (Ifges(5,5) + (Ifges(6,1) * t186 + t259) * t269 + (Ifges(6,2) * t189 + t260) * t270) * t149; -t217 * t101 + t137 * t102 - t297 * t60 + t84 * t61 + (-t297 * t6 + t5 * t84 + t31) * m(7) + (t137 * t29 - t217 * t30 + t47) * m(6) + t292; t146 - t51 * (mrSges(7,1) * t84 + mrSges(7,2) * t297) + (Ifges(7,1) * t297 - t267) * t277 + t25 * t276 + (Ifges(7,5) * t297 - Ifges(7,6) * t84) * t272 - t5 * t60 + t6 * t61 + (t297 * t5 + t6 * t84) * mrSges(7,3) + (-Ifges(7,2) * t84 + t26 + t81) * t279 + t286;];
tauc  = t3(:);
