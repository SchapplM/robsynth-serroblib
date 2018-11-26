% Calculate vector of centrifugal and coriolis load on the joints for
% S6PRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:04
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PRPRRR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:04:43
% EndTime: 2018-11-23 15:04:49
% DurationCPUTime: 6.44s
% Computational Cost: add. (9456->441), mult. (24819->621), div. (0->0), fcn. (20019->12), ass. (0->212)
t181 = sin(pkin(12));
t274 = pkin(8) + qJ(3);
t168 = t274 * t181;
t183 = cos(pkin(12));
t169 = t274 * t183;
t187 = sin(qJ(4));
t191 = cos(qJ(4));
t132 = -t187 * t168 + t191 * t169;
t164 = t181 * t191 + t183 * t187;
t182 = sin(pkin(6));
t192 = cos(qJ(2));
t251 = t182 * t192;
t204 = t164 * t251;
t300 = qJD(1) * t204 - t164 * qJD(3) - qJD(4) * t132;
t250 = t183 * t191;
t163 = -t181 * t187 + t250;
t203 = t163 * t251;
t244 = qJD(4) * t191;
t299 = -qJD(1) * t203 - t168 * t244 + qJD(3) * t250 + (-qJD(3) * t181 - qJD(4) * t169) * t187;
t156 = t164 * qJD(4);
t312 = -pkin(9) * t156 + t299;
t155 = t163 * qJD(4);
t311 = -pkin(9) * t155 + t300;
t188 = sin(qJ(2));
t248 = qJD(1) * t182;
t236 = t188 * t248;
t279 = pkin(4) * t156;
t310 = -t236 + t279;
t186 = sin(qJ(5));
t190 = cos(qJ(5));
t131 = -t191 * t168 - t169 * t187;
t108 = -pkin(9) * t164 + t131;
t109 = pkin(9) * t163 + t132;
t215 = t190 * t108 - t109 * t186;
t304 = qJD(5) * t215 + t311 * t186 + t312 * t190;
t211 = t190 * t163 - t164 * t186;
t85 = qJD(5) * t211 + t155 * t190 - t156 * t186;
t125 = t163 * t186 + t164 * t190;
t86 = qJD(5) * t125 + t155 * t186 + t190 * t156;
t309 = pkin(5) * t86 - pkin(10) * t85 + t310;
t180 = qJD(4) + qJD(5);
t185 = sin(qJ(6));
t189 = cos(qJ(6));
t153 = t163 * qJD(2);
t154 = t164 * qJD(2);
t212 = t153 * t186 + t190 * t154;
t100 = t180 * t189 - t185 * t212;
t146 = qJD(2) * t155;
t147 = qJD(2) * t156;
t231 = t190 * t153 - t154 * t186;
t77 = qJD(5) * t231 + t146 * t190 - t147 * t186;
t50 = qJD(6) * t100 + t189 * t77;
t101 = t180 * t185 + t189 * t212;
t51 = -qJD(6) * t101 - t185 * t77;
t18 = -mrSges(7,1) * t51 + mrSges(7,2) * t50;
t246 = qJD(1) * t192;
t235 = t182 * t246;
t162 = (qJD(3) + t235) * qJD(2);
t167 = qJD(2) * qJ(3) + t236;
t184 = cos(pkin(6));
t247 = qJD(1) * t184;
t172 = t183 * t247;
t267 = pkin(8) * qJD(2);
t129 = t172 + (-t167 - t267) * t181;
t136 = t183 * t167 + t181 * t247;
t130 = t183 * t267 + t136;
t91 = t129 * t187 + t130 * t191;
t56 = -qJD(4) * t91 - t162 * t164;
t200 = -pkin(9) * t146 + t56;
t80 = pkin(9) * t153 + t91;
t255 = t190 * t80;
t90 = t191 * t129 - t130 * t187;
t79 = -pkin(9) * t154 + t90;
t76 = qJD(4) * pkin(4) + t79;
t36 = t186 * t76 + t255;
t55 = t129 * t244 + t162 * t250 + (-qJD(4) * t130 - t162 * t181) * t187;
t53 = -pkin(9) * t147 + t55;
t9 = qJD(5) * t36 + t186 * t53 - t190 * t200;
t308 = m(7) * t9 + t18;
t112 = qJD(6) - t231;
t221 = Ifges(7,5) * t189 - Ifges(7,6) * t185;
t268 = Ifges(7,4) * t189;
t223 = -Ifges(7,2) * t185 + t268;
t269 = Ifges(7,4) * t185;
t225 = Ifges(7,1) * t189 - t269;
t226 = mrSges(7,1) * t185 + mrSges(7,2) * t189;
t280 = t189 / 0.2e1;
t281 = -t185 / 0.2e1;
t287 = t101 / 0.2e1;
t256 = t186 * t80;
t35 = t190 * t76 - t256;
t33 = -pkin(5) * t180 - t35;
t270 = Ifges(7,4) * t101;
t46 = Ifges(7,2) * t100 + Ifges(7,6) * t112 + t270;
t99 = Ifges(7,4) * t100;
t47 = t101 * Ifges(7,1) + t112 * Ifges(7,5) + t99;
t307 = t33 * t226 + t46 * t281 + t47 * t280 + t100 * t223 / 0.2e1 + t225 * t287 + t112 * t221 / 0.2e1;
t65 = t108 * t186 + t109 * t190;
t174 = -pkin(3) * t183 - pkin(2);
t140 = -pkin(4) * t163 + t174;
t66 = -pkin(5) * t211 - pkin(10) * t125 + t140;
t31 = -t185 * t65 + t189 * t66;
t306 = qJD(6) * t31 + t309 * t185 + t189 * t304;
t32 = t185 * t66 + t189 * t65;
t305 = -qJD(6) * t32 - t185 * t304 + t309 * t189;
t303 = qJD(5) * t65 + t312 * t186 - t311 * t190;
t34 = pkin(10) * t180 + t36;
t219 = qJD(3) - t235;
t148 = qJD(2) * t174 + t219;
t121 = -pkin(4) * t153 + t148;
t54 = -pkin(5) * t231 - pkin(10) * t212 + t121;
t17 = t185 * t54 + t189 * t34;
t16 = -t185 * t34 + t189 * t54;
t260 = t16 * t189;
t218 = t17 * t185 + t260;
t302 = t218 * mrSges(7,3);
t249 = t181 ^ 2 + t183 ^ 2;
t233 = t249 * mrSges(4,3);
t110 = t147 * mrSges(5,1) + t146 * mrSges(5,2);
t78 = qJD(5) * t212 + t146 * t186 + t190 * t147;
t39 = t78 * mrSges(6,1) + t77 * mrSges(6,2);
t301 = -t110 - t39;
t254 = mrSges(6,1) * t180 + mrSges(7,1) * t100 - mrSges(7,2) * t101 - mrSges(6,3) * t212;
t298 = -t16 * t185 + t17 * t189;
t245 = qJD(2) * t188;
t234 = t182 * t245;
t170 = qJD(1) * t234;
t128 = pkin(4) * t147 + t170;
t29 = pkin(5) * t78 - pkin(10) * t77 + t128;
t8 = qJD(5) * t35 + t186 * t200 + t190 * t53;
t2 = qJD(6) * t16 + t185 * t29 + t189 * t8;
t3 = -qJD(6) * t17 - t185 * t8 + t189 * t29;
t297 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t50 + Ifges(7,6) * t51;
t111 = Ifges(6,4) * t231;
t296 = t111 / 0.2e1 + t212 * Ifges(6,1) / 0.2e1;
t82 = pkin(5) * t212 - pkin(10) * t231;
t294 = t50 / 0.2e1;
t293 = t51 / 0.2e1;
t292 = t78 / 0.2e1;
t291 = t215 * t9;
t252 = t182 * t188;
t150 = -t181 * t252 + t183 * t184;
t151 = t181 * t184 + t183 * t252;
t113 = t150 * t191 - t151 * t187;
t114 = t150 * t187 + t151 * t191;
t214 = t190 * t113 - t114 * t186;
t290 = t214 * t9;
t289 = -t100 / 0.2e1;
t288 = -t101 / 0.2e1;
t286 = -t112 / 0.2e1;
t284 = -t154 / 0.2e1;
t283 = t155 / 0.2e1;
t282 = -t156 / 0.2e1;
t278 = t189 * t2;
t277 = t3 * t185;
t276 = t35 * mrSges(6,3);
t275 = t36 * mrSges(6,3);
t272 = Ifges(5,4) * t154;
t263 = t231 * Ifges(6,2);
t243 = qJD(6) * t185;
t242 = qJD(6) * t189;
t228 = -mrSges(4,1) * t183 + mrSges(4,2) * t181;
t81 = -mrSges(6,1) * t231 + mrSges(6,2) * t212;
t240 = -mrSges(5,1) * t153 + mrSges(5,2) * t154 + qJD(2) * t228 + t81;
t237 = t182 ^ 2 * t246;
t232 = t249 * t162;
t229 = -t2 * t185 - t3 * t189;
t227 = mrSges(7,1) * t189 - mrSges(7,2) * t185;
t224 = Ifges(7,1) * t185 + t268;
t222 = Ifges(7,2) * t189 + t269;
t220 = Ifges(7,5) * t185 + Ifges(7,6) * t189;
t62 = -mrSges(7,2) * t112 + mrSges(7,3) * t100;
t63 = mrSges(7,1) * t112 - mrSges(7,3) * t101;
t216 = -t185 * t63 + t189 * t62;
t70 = t113 * t186 + t114 * t190;
t213 = -(-t167 * t181 + t172) * t181 + t136 * t183;
t102 = -mrSges(6,2) * t180 + mrSges(6,3) * t231;
t210 = -t102 - t216;
t58 = -t185 * t70 - t189 * t251;
t209 = t185 * t251 - t189 * t70;
t202 = ((-qJD(2) * pkin(2) + t219) * t188 + t192 * t213) * t182;
t12 = t50 * Ifges(7,4) + t51 * Ifges(7,2) + t78 * Ifges(7,6);
t13 = t50 * Ifges(7,1) + t51 * Ifges(7,4) + t78 * Ifges(7,5);
t201 = -t8 * mrSges(6,2) + mrSges(7,3) * t278 + t185 * t13 / 0.2e1 + t12 * t280 + t224 * t294 + t222 * t293 + t220 * t292 - Ifges(6,6) * t78 + Ifges(6,5) * t77 + (-mrSges(6,1) - t227) * t9 + t307 * qJD(6);
t198 = t17 * mrSges(7,2) - t112 * Ifges(7,3) - t101 * Ifges(7,5) - t100 * Ifges(7,6) + t180 * Ifges(6,6) + t263 / 0.2e1 + Ifges(6,4) * t212 - t121 * mrSges(6,1) - t16 * mrSges(7,1);
t25 = mrSges(7,1) * t78 - mrSges(7,3) * t50;
t26 = -mrSges(7,2) * t78 + mrSges(7,3) * t51;
t197 = -t63 * t242 - t62 * t243 - t185 * t25 + t189 * t26 + m(7) * (-t16 * t242 - t17 * t243 - t277 + t278);
t196 = -t263 / 0.2e1 - t198;
t195 = t121 * mrSges(6,2) + t180 * Ifges(6,5) + t296 + t307;
t194 = t195 + t296;
t193 = qJD(2) ^ 2;
t149 = Ifges(5,4) * t153;
t139 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t154;
t138 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t153;
t116 = t154 * Ifges(5,1) + Ifges(5,5) * qJD(4) + t149;
t115 = t153 * Ifges(5,2) + Ifges(5,6) * qJD(4) + t272;
t89 = -qJD(2) * t204 - qJD(4) * t114;
t88 = qJD(2) * t203 + qJD(4) * t113;
t73 = Ifges(7,3) * t78;
t60 = pkin(4) * t154 + t82;
t38 = t190 * t79 - t256;
t37 = t186 * t79 + t255;
t24 = qJD(5) * t70 + t186 * t88 - t190 * t89;
t23 = qJD(5) * t214 + t186 * t89 + t190 * t88;
t22 = t185 * t82 + t189 * t35;
t21 = -t185 * t35 + t189 * t82;
t20 = t185 * t60 + t189 * t38;
t19 = -t185 * t38 + t189 * t60;
t15 = qJD(6) * t209 - t185 * t23 + t189 * t234;
t14 = qJD(6) * t58 + t185 * t234 + t189 * t23;
t1 = [t23 * t102 + t88 * t138 + t89 * t139 + t14 * t62 + t15 * t63 - t214 * t18 + t58 * t25 - t209 * t26 - t254 * t24 + (-t214 * t77 - t70 * t78) * mrSges(6,3) + (-t113 * t146 - t114 * t147) * mrSges(5,3) + ((-mrSges(3,1) * t193 + qJD(2) * t240) * t188 + ((-mrSges(3,2) + t233) * t193 + t301) * t192) * t182 + m(7) * (t14 * t17 + t15 * t16 - t2 * t209 + t24 * t33 + t3 * t58 - t290) + m(6) * (t23 * t36 - t24 * t35 - t290 + t70 * t8 + (t121 * t245 - t128 * t192) * t182) + m(5) * (t113 * t56 + t114 * t55 + t88 * t91 + t89 * t90 + (t148 * t182 - t237) * t245) + m(4) * ((-t150 * t181 + t151 * t183) * t162 + (-t188 * t237 + t202) * qJD(2)); (t225 * t294 + t223 * t293 + t221 * t292 + Ifges(6,1) * t77 - Ifges(6,4) * t78 + t128 * mrSges(6,2) + t12 * t281 + t13 * t280 + (mrSges(6,3) + t226) * t9 + t229 * mrSges(7,3) + (-t189 * t46 / 0.2e1 + t47 * t281 + t220 * t286 + t222 * t289 + t224 * t288 + t33 * t227 - t298 * mrSges(7,3)) * qJD(6)) * t125 + (qJD(2) * qJD(3) * t249 + t232) * mrSges(4,3) + t81 * t279 + (t146 * t164 + t154 * t283) * Ifges(5,1) + t305 * t63 + (t305 * t16 + t17 * t306 + t2 * t32 + t3 * t31 + t303 * t33 - t291) * m(7) + t306 * t62 - t254 * t303 + t304 * t102 + (-pkin(2) * t170 + qJ(3) * t232 - qJD(1) * t202 + qJD(3) * t213) * m(4) + (-t215 * t77 - t35 * t85 - t36 * t86 - t65 * t78) * mrSges(6,3) - t215 * t18 - (t73 / 0.2e1 - Ifges(6,4) * t77 + t128 * mrSges(6,1) - t8 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2)) * t78 + t297) * t211 + (-t163 * t147 + t153 * t282) * Ifges(5,2) + (t163 * t146 - t147 * t164 + t153 * t283 + t154 * t282) * Ifges(5,4) + (-t131 * t146 - t132 * t147 - t155 * t90 - t156 * t91 + t163 * t55 - t164 * t56) * mrSges(5,3) + (t121 * t310 + t128 * t140 - t303 * t35 + t304 * t36 + t65 * t8 - t291) * m(6) + t115 * t282 + t116 * t283 + (t194 - t302) * t85 + t299 * t138 + t300 * t139 + (t131 * t56 + t132 * t55 - t148 * t236 + t170 * t174 + t299 * t91 + t300 * t90) * m(5) + qJD(4) * (Ifges(5,5) * t155 - Ifges(5,6) * t156) / 0.2e1 + t148 * (mrSges(5,1) * t156 + mrSges(5,2) * t155) + (-t192 * qJD(2) * t233 + ((-mrSges(5,1) * t163 + mrSges(5,2) * t164 + t228) * qJD(2) - t240) * t188) * t248 + t196 * t86 + t174 * t110 + t140 * t39 + t31 * t25 + t32 * t26; -t153 * t138 + t154 * t139 + t185 * t26 + t189 * t25 + t254 * t212 + t216 * qJD(6) + (m(4) + m(5)) * t170 - t193 * t233 + t210 * t231 - m(5) * (t153 * t91 - t154 * t90) - m(4) * t213 * qJD(2) + (t112 * t298 - t212 * t33 - t229) * m(7) + (t212 * t35 - t231 * t36 + t128) * m(6) - t301; (-t112 * t260 + (-t112 * t17 - t3) * t185) * mrSges(7,3) + t254 * t37 + t201 - (-Ifges(5,2) * t154 + t116 + t149) * t153 / 0.2e1 + (Ifges(5,1) * t153 - t272) * t284 + (-t154 * t81 + (-t186 * t78 - t190 * t77) * mrSges(6,3) + (-t254 * t186 - t210 * t190 + m(7) * (t186 * t33 + t190 * t298)) * qJD(5) + (0.2e1 * t121 * t284 + t186 * t8 - t190 * t9 + (-t186 * t35 + t190 * t36) * qJD(5)) * m(6)) * pkin(4) - m(7) * (t16 * t19 + t17 * t20 + t33 * t37) + (-t196 + t275) * t212 + (-t194 + t276) * t231 + t197 * (pkin(4) * t186 + pkin(10)) - m(6) * (-t35 * t37 + t36 * t38) - qJD(4) * (Ifges(5,5) * t153 - Ifges(5,6) * t154) / 0.2e1 - t148 * (mrSges(5,1) * t154 + mrSges(5,2) * t153) + t154 * t115 / 0.2e1 - Ifges(5,6) * t147 + Ifges(5,5) * t146 - t90 * t138 + t91 * t139 - t38 * t102 + t308 * (-pkin(4) * t190 - pkin(5)) + (t153 * t90 + t154 * t91) * mrSges(5,3) - t55 * mrSges(5,2) + t56 * mrSges(5,1) - t20 * t62 - t19 * t63; (-qJD(6) * t218 - t277) * mrSges(7,3) + t254 * t36 + t201 + t197 * pkin(10) + (t198 + t275) * t212 + (-t111 / 0.2e1 + t276 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t212 + t302 - t195) * t231 - t35 * t102 - m(7) * (t16 * t21 + t17 * t22 + t33 * t36) - t22 * t62 - t21 * t63 - t308 * pkin(5); t73 - t33 * (mrSges(7,1) * t101 + mrSges(7,2) * t100) + (Ifges(7,1) * t100 - t270) * t288 + t46 * t287 + (Ifges(7,5) * t100 - Ifges(7,6) * t101) * t286 - t16 * t62 + t17 * t63 + (t100 * t16 + t101 * t17) * mrSges(7,3) + (-Ifges(7,2) * t101 + t47 + t99) * t289 + t297;];
tauc  = t1(:);
