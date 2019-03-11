% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRRPRP4
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
% Datum: 2019-03-08 21:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRP4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:40:17
% EndTime: 2019-03-08 21:40:29
% DurationCPUTime: 6.75s
% Computational Cost: add. (3208->478), mult. (7947->637), div. (0->0), fcn. (4802->8), ass. (0->226)
t311 = Ifges(6,4) + Ifges(7,4);
t312 = Ifges(6,1) + Ifges(7,1);
t299 = Ifges(6,5) + Ifges(7,5);
t310 = Ifges(6,2) + Ifges(7,2);
t298 = Ifges(6,6) + Ifges(7,6);
t155 = sin(qJ(5));
t158 = cos(qJ(5));
t161 = -pkin(3) - pkin(9);
t154 = cos(pkin(6));
t159 = cos(qJ(3));
t244 = t154 * t159;
t137 = qJD(1) * t244;
t156 = sin(qJ(3));
t157 = sin(qJ(2));
t153 = sin(pkin(6));
t234 = qJD(1) * t153;
t208 = t157 * t234;
t121 = qJD(2) * pkin(8) + t208;
t200 = pkin(4) * qJD(2) + t121;
t176 = t200 * t156;
t65 = t137 - t176;
t48 = qJD(3) * t161 + qJD(4) - t65;
t201 = -qJ(4) * t156 - pkin(2);
t107 = t159 * t161 + t201;
t160 = cos(qJ(2));
t207 = t160 * t234;
t67 = qJD(2) * t107 - t207;
t15 = -t155 * t67 + t158 * t48;
t16 = t155 * t48 + t158 * t67;
t178 = t15 * t155 - t158 * t16;
t191 = mrSges(7,1) * t158 - mrSges(7,2) * t155;
t193 = mrSges(6,1) * t158 - mrSges(6,2) * t155;
t231 = qJD(2) * t156;
t142 = qJD(5) + t231;
t227 = qJD(3) * t158;
t229 = qJD(2) * t159;
t112 = -t155 * t229 + t227;
t8 = -qJ(6) * t112 + t15;
t7 = pkin(5) * t142 + t8;
t111 = -qJD(3) * t155 - t158 * t229;
t9 = qJ(6) * t111 + t16;
t194 = t155 * t7 - t158 * t9;
t262 = -t158 / 0.2e1;
t264 = -t155 / 0.2e1;
t266 = -t142 / 0.2e1;
t268 = -t112 / 0.2e1;
t269 = -t111 / 0.2e1;
t307 = t311 * t111;
t278 = t112 * t312 + t299 * t142 + t307;
t305 = t311 * t112;
t279 = t111 * t310 + t142 * t298 + t305;
t304 = t311 * t158;
t289 = t155 * t312 + t304;
t303 = t311 * t155;
t291 = t158 * t310 + t303;
t152 = qJD(3) * qJ(4);
t233 = qJD(1) * t156;
t136 = t154 * t233;
t86 = t159 * t121 + t136;
t66 = pkin(4) * t229 + t86;
t53 = t152 + t66;
t30 = -pkin(5) * t111 + qJD(6) + t53;
t309 = t262 * t279 + t264 * t278 + t30 * t191 + t53 * t193 + (t155 * t299 + t158 * t298) * t266 + t291 * t269 + t289 * t268 + t178 * mrSges(6,3) + t194 * mrSges(7,3);
t308 = -t229 / 0.2e1;
t249 = qJD(2) * pkin(2);
t122 = -t207 - t249;
t281 = qJD(3) / 0.2e1;
t282 = -qJD(3) / 0.2e1;
t283 = -qJD(2) / 0.2e1;
t70 = -t152 - t86;
t125 = -pkin(3) * t159 + t201;
t87 = qJD(2) * t125 - t207;
t306 = t86 * mrSges(4,3) + t87 * mrSges(5,2) + Ifges(5,5) * t282 + (-Ifges(5,6) * t156 - t159 * Ifges(5,3)) * t283 + Ifges(4,6) * t281 + (Ifges(4,4) * t156 + t159 * Ifges(4,2)) * qJD(2) / 0.2e1 - t122 * mrSges(4,1) - t70 * mrSges(5,1) + t309;
t302 = -Ifges(4,1) / 0.2e1;
t270 = pkin(4) + pkin(8);
t301 = Ifges(4,4) * t308;
t280 = mrSges(5,1) + mrSges(4,3);
t300 = mrSges(5,2) - mrSges(4,1);
t220 = qJD(2) * qJD(3);
t202 = t159 * t220;
t219 = qJD(3) * qJD(5);
t222 = qJD(5) * t159;
t228 = qJD(3) * t156;
t77 = -t155 * t219 + (t155 * t228 - t158 * t222) * qJD(2);
t204 = t155 * t222;
t169 = t156 * t227 + t204;
t78 = qJD(2) * t169 - t158 * t219;
t297 = t202 * t298 + t310 * t78 + t311 * t77;
t296 = t299 * t202 + t311 * t78 + t312 * t77;
t224 = qJD(5) * t155;
t146 = pkin(3) * t231;
t179 = pkin(9) * t156 - qJ(4) * t159;
t92 = qJD(2) * t179 + t146;
t23 = -t155 * t92 + t158 * t66;
t240 = qJ(6) - t161;
t295 = -qJD(6) * t158 + t224 * t240 - (-qJ(6) * t155 * t156 + pkin(5) * t159) * qJD(2) - t23;
t124 = t240 * t158;
t24 = t155 * t66 + t158 * t92;
t294 = -qJ(6) * t158 * t231 - qJD(5) * t124 - qJD(6) * t155 - t24;
t209 = -pkin(5) * t158 - pkin(4);
t223 = qJD(5) * t158;
t293 = pkin(5) * t223 + qJD(4) - t137 - (qJD(2) * t209 - t121) * t156;
t85 = t121 * t156 - t137;
t292 = -qJD(4) - t85;
t290 = -t155 * t310 + t304;
t288 = t158 * t312 - t303;
t267 = t112 / 0.2e1;
t284 = t142 / 0.2e1;
t232 = qJD(2) * t153;
t203 = qJD(1) * t232;
t198 = t160 * t203;
t31 = t156 * t198 + (t159 * t200 + t136) * qJD(3);
t225 = qJD(4) * t156;
t166 = qJD(3) * t179 - t225;
t235 = qJD(3) * t146 + t157 * t203;
t46 = qJD(2) * t166 + t235;
t3 = t155 * t31 + t158 * t46 + t48 * t223 - t224 * t67;
t4 = -qJD(5) * t16 - t155 * t46 + t158 * t31;
t276 = t155 * t3 + t158 * t4;
t211 = Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1;
t212 = Ifges(6,6) / 0.2e1 + Ifges(7,6) / 0.2e1;
t215 = Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1;
t69 = -qJD(3) * pkin(3) - t292;
t275 = -t212 * t111 - t215 * t112 - t211 * t142 - t122 * mrSges(4,2) - t15 * mrSges(6,1) - t69 * mrSges(5,1) - t7 * mrSges(7,1) - t85 * mrSges(4,3) - Ifges(5,4) * t282 - (-t156 * Ifges(5,2) - Ifges(5,6) * t159) * t283 + t231 * t302 - Ifges(4,5) * t281 + t301 + t16 * mrSges(6,2) + t87 * mrSges(5,3) + t9 * mrSges(7,2) + t298 * t269 - (Ifges(7,3) + Ifges(6,3)) * t284 - t299 * t267;
t273 = -m(5) / 0.2e1;
t205 = t160 * t232;
t170 = qJD(3) * t154 + t205;
t226 = qJD(3) * t159;
t45 = t121 * t226 + t170 * t233;
t245 = t153 * t157;
t210 = t156 * t245;
t96 = t210 - t244;
t258 = t45 * t96;
t49 = mrSges(7,1) * t202 - mrSges(7,3) * t77;
t50 = mrSges(6,1) * t202 - mrSges(6,3) * t77;
t257 = t49 + t50;
t51 = -mrSges(7,2) * t202 + mrSges(7,3) * t78;
t52 = -mrSges(6,2) * t202 + mrSges(6,3) * t78;
t256 = t51 + t52;
t81 = -mrSges(7,2) * t142 + mrSges(7,3) * t111;
t82 = -mrSges(6,2) * t142 + mrSges(6,3) * t111;
t255 = t81 + t82;
t83 = mrSges(7,1) * t142 - mrSges(7,3) * t112;
t84 = mrSges(6,1) * t142 - mrSges(6,3) * t112;
t254 = t83 + t84;
t129 = -mrSges(5,1) * t229 - qJD(3) * mrSges(5,3);
t60 = -mrSges(6,1) * t111 + mrSges(6,2) * t112;
t246 = -t129 + t60;
t131 = t270 * t156;
t57 = t158 * t107 + t155 * t131;
t243 = t155 * t160;
t242 = t158 * t159;
t241 = t158 * t160;
t115 = (mrSges(5,2) * t159 - mrSges(5,3) * t156) * qJD(2);
t239 = t115 + (-mrSges(4,1) * t159 + mrSges(4,2) * t156) * qJD(2);
t238 = qJD(3) * t137 + t159 * t198;
t237 = -qJD(3) * t300 - t231 * t280;
t128 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t229;
t236 = t128 - t129;
t132 = t270 * t159;
t230 = qJD(2) * t157;
t221 = qJD(6) * t159;
t218 = pkin(8) * t45 / 0.2e1;
t59 = -mrSges(7,1) * t111 + mrSges(7,2) * t112;
t217 = -t59 - t246;
t216 = Ifges(4,5) / 0.2e1 - Ifges(5,4) / 0.2e1;
t214 = -Ifges(4,6) / 0.2e1 + Ifges(5,5) / 0.2e1;
t213 = -0.3e1 / 0.2e1 * Ifges(5,6) - 0.3e1 / 0.2e1 * Ifges(4,4);
t206 = t153 * t230;
t25 = -t78 * mrSges(7,1) + t77 * mrSges(7,2);
t199 = qJ(6) * t159 - t107;
t197 = -t207 / 0.2e1;
t1 = pkin(5) * t202 - qJ(6) * t77 - qJD(6) * t112 + t4;
t2 = qJ(6) * t78 + qJD(6) * t111 + t3;
t196 = -t1 * t158 - t155 * t2;
t192 = mrSges(6,1) * t155 + mrSges(6,2) * t158;
t190 = mrSges(7,1) * t155 + mrSges(7,2) * t158;
t44 = -t121 * t228 + t238;
t174 = t153 * t241 - t155 * t96;
t63 = t153 * t243 + t158 * t96;
t97 = t154 * t156 + t159 * t245;
t172 = -m(4) * t86 + m(5) * t70 - t236;
t171 = -qJ(4) * t226 - t225;
t119 = t270 * t226;
t148 = pkin(3) * t228;
t88 = t148 + t166;
t10 = -t107 * t224 + t155 * t119 + t131 * t223 + t158 * t88;
t168 = -t155 * t254 + t158 * t255;
t167 = t4 * mrSges(6,1) + t1 * mrSges(7,1) - t3 * mrSges(6,2) - t2 * mrSges(7,2);
t151 = qJD(3) * qJD(4);
t27 = -qJD(3) * t176 + t151 + t238;
t162 = qJD(2) ^ 2;
t144 = pkin(5) * t155 + qJ(4);
t140 = Ifges(6,3) * t202;
t139 = Ifges(7,3) * t202;
t123 = t240 * t155;
t118 = t270 * t228;
t117 = -qJ(4) * t229 + t146;
t114 = t158 * t131;
t105 = (mrSges(4,1) * t156 + mrSges(4,2) * t159) * t220;
t104 = (-mrSges(5,2) * t156 - mrSges(5,3) * t159) * t220;
t103 = t158 * t119;
t95 = pkin(5) * t242 + t132;
t94 = t148 + t171;
t80 = (t156 * t243 + t157 * t158) * t234;
t79 = (-t155 * t157 + t156 * t241) * t234;
t75 = Ifges(6,5) * t77;
t74 = Ifges(7,5) * t77;
t73 = Ifges(6,6) * t78;
t72 = Ifges(7,6) * t78;
t68 = -pkin(5) * t204 + (-pkin(8) + t209) * t228;
t62 = -qJD(3) * t210 + t159 * t170;
t61 = qJD(3) * t97 + t156 * t205;
t58 = qJD(2) * t171 + t235;
t56 = -t107 * t155 + t114;
t43 = -qJ(6) * t242 + t57;
t33 = -t151 - t44;
t32 = pkin(5) * t156 + t155 * t199 + t114;
t26 = -mrSges(6,1) * t78 + mrSges(6,2) * t77;
t14 = -pkin(5) * t78 + t27;
t13 = qJD(5) * t63 + t155 * t61 + t158 * t206;
t12 = qJD(5) * t174 - t155 * t206 + t158 * t61;
t11 = -qJD(5) * t57 - t155 * t88 + t103;
t6 = qJ(6) * t169 - t158 * t221 + t10;
t5 = pkin(5) * t226 + t103 + t199 * t223 + (-qJ(6) * t228 - qJD(5) * t131 + t221 - t88) * t155;
t17 = [(t25 + t26) * t97 - t256 * t174 + t257 * t63 - t237 * t61 + t255 * t13 + t254 * t12 + (t128 - t217) * t62 + (-t162 * t157 * mrSges(3,1) + (-mrSges(3,2) * t162 - t104 - t105) * t160) * t153 + (t239 * t245 + t280 * qJD(3) * (-t156 * t97 + t159 * t96)) * qJD(2) + m(4) * (t44 * t97 + t258 + t61 * t85 + t62 * t86 + (t122 - t207) * t206) + m(5) * (-t33 * t97 + t258 + t61 * t69 - t62 * t70 + (-t160 * t58 + t230 * t87) * t153) + m(6) * (t12 * t15 + t13 * t16 - t174 * t3 + t27 * t97 + t4 * t63 + t53 * t62) + m(7) * (t1 * t63 + t12 * t7 + t13 * t9 + t14 * t97 - t174 * t2 + t30 * t62); -t118 * t60 + t125 * t104 + t132 * t26 + t94 * t115 + t95 * t25 - pkin(2) * t105 + t6 * t81 + t10 * t82 + t5 * t83 + t11 * t84 + t68 * t59 + t56 * t50 + t57 * t52 + t32 * t49 + t43 * t51 + (t139 / 0.2e1 + t140 / 0.2e1 + t74 / 0.2e1 + t75 / 0.2e1 + t72 / 0.2e1 + t73 / 0.2e1 - t58 * mrSges(5,3) + t212 * t78 + t215 * t77 + t280 * t45 + (mrSges(4,2) * t230 + t160 * t237) * t234 + 0.2e1 * (t197 * t69 + t218) * m(5) + 0.2e1 * (t197 * t85 + t218) * m(4) + (pkin(8) * t172 + t214 * qJD(3) + t213 * t231 - t306) * qJD(3) + t167) * t156 + m(6) * (t10 * t16 + t11 * t15 - t118 * t53 + t132 * t27 + t3 * t57 + t4 * t56) + m(7) * (t1 * t32 + t14 * t95 + t2 * t43 + t30 * t68 + t5 * t7 + t6 * t9) + (t27 * t193 + t14 * t191 + t58 * mrSges(5,2) - t33 * mrSges(5,1) + t44 * mrSges(4,3) + (t1 * t155 - t158 * t2) * mrSges(7,3) + (t155 * t4 - t158 * t3) * mrSges(6,3) + (-mrSges(4,1) * t230 + (-m(6) * t53 - m(7) * t30 + t172 - t59 - t60) * t160) * t234 + (-t53 * t192 - t30 * t190 + (t155 * t9 + t158 * t7) * mrSges(7,3) + (t15 * t158 + t155 * t16) * mrSges(6,3) + t290 * t269 + t288 * t268 + (t155 * t298 - t158 * t299) * t284 + t279 * t155 / 0.2e1) * qJD(5) + (((-t155 * t215 - t158 * t212 - t213) * t159 + (-0.3e1 / 0.2e1 * Ifges(5,3) - 0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(5,2) + t211) * t156) * qJD(2) + t216 * qJD(3) - t275) * qJD(3) - t289 * t77 / 0.2e1 - t291 * t78 / 0.2e1 + t296 * t264 + (qJD(5) * t278 + t297) * t262 + (m(4) * t44 + 0.2e1 * t33 * t273 + (m(4) * t85 + m(5) * t69 - t237) * qJD(3)) * pkin(8)) * t159 + m(5) * (t125 * t58 + t87 * t94) - m(6) * (t15 * t79 + t16 * t80) - m(7) * (t7 * t79 + t80 * t9) + 0.2e1 * (t87 * t273 + (-t122 / 0.2e1 - t249 / 0.2e1) * m(4)) * t208 - t255 * t80 - t254 * t79 - t239 * t208; t144 * t25 - t123 * t51 - t124 * t49 + ((-m(6) * t178 - t155 * t84 + t158 * t82) * t161 + t309) * qJD(5) - m(6) * (t15 * t23 + t16 * t24 + t53 * t65) - t117 * t115 - t24 * t82 - t23 * t84 - t65 * t60 - t44 * mrSges(4,2) - t33 * mrSges(5,3) + qJ(4) * t26 + t14 * t190 + t27 * t192 + t294 * t81 + (-t1 * t124 - t123 * t2 + t14 * t144 + t293 * t30 + t294 * t9 + t295 * t7) * m(7) + t295 * t83 + t288 * t77 / 0.2e1 + t290 * t78 / 0.2e1 + (-pkin(3) * t45 - qJ(4) * t33 - t117 * t87 + t292 * t70 - t69 * t86) * m(5) + m(6) * (qJ(4) * t27 + qJD(4) * t53 + t161 * t276) + ((t301 + (-pkin(3) * mrSges(5,1) - t155 * t212 + t158 * t215 + t216) * qJD(3) + Ifges(5,6) * t308 + t275) * t159 + ((-qJ(4) * mrSges(5,1) + t214) * qJD(3) + (-Ifges(5,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + t302 + Ifges(4,2) / 0.2e1) * t229 + (Ifges(5,6) / 0.2e1 + Ifges(4,4) / 0.2e1) * t231 + t306) * t156) * qJD(2) - t276 * mrSges(6,3) + (t155 * t52 + t158 * t50) * t161 + t293 * t59 + t296 * t158 / 0.2e1 + t297 * t264 + t300 * t45 + t196 * mrSges(7,3) + t236 * t85 + t237 * t86 + t246 * qJD(4); t257 * t158 + t256 * t155 + t217 * qJD(3) + t168 * qJD(5) + (mrSges(5,1) * t226 + (t115 + t168) * t156) * qJD(2) + (-qJD(3) * t30 - t142 * t194 - t196) * m(7) + (-qJD(3) * t53 - t142 * t178 + t276) * m(6) + (qJD(3) * t70 + t231 * t87 + t45) * m(5); t139 + t140 + t74 + t75 + t72 + t73 + t167 - t53 * (mrSges(6,1) * t112 + mrSges(6,2) * t111) - t30 * (mrSges(7,1) * t112 + mrSges(7,2) * t111) - t8 * t81 - t15 * t82 + t9 * t83 + t16 * t84 + (-t112 * t59 + t49) * pkin(5) + (t111 * t7 + t112 * t9) * mrSges(7,3) + (t111 * t15 + t112 * t16) * mrSges(6,3) + (-(-t7 + t8) * t9 + (-t112 * t30 + t1) * pkin(5)) * m(7) + (t111 * t312 - t305) * t268 + t279 * t267 + (t111 * t299 - t112 * t298) * t266 + (-t112 * t310 + t278 + t307) * t269; -t111 * t81 + t112 * t83 + 0.2e1 * (t14 / 0.2e1 + t9 * t269 + t7 * t267) * m(7) + t25;];
tauc  = t17(:);
