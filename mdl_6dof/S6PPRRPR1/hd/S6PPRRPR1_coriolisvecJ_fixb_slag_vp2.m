% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
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
% Datum: 2019-03-08 18:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PPRRPR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRPR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRPR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:45:24
% EndTime: 2019-03-08 18:45:35
% DurationCPUTime: 5.26s
% Computational Cost: add. (5276->475), mult. (14380->704), div. (0->0), fcn. (11910->14), ass. (0->229)
t182 = sin(qJ(4));
t185 = cos(qJ(4));
t203 = pkin(4) * t182 - qJ(5) * t185;
t135 = qJD(4) * t203 - qJD(5) * t182;
t173 = sin(pkin(13));
t177 = cos(pkin(13));
t224 = qJD(4) * t182;
t219 = pkin(9) * t224;
t108 = t177 * t135 + t173 * t219;
t239 = t173 * t185;
t180 = cos(pkin(6));
t166 = qJD(1) * t180 + qJD(2);
t175 = sin(pkin(7));
t179 = cos(pkin(7));
t186 = cos(qJ(3));
t178 = cos(pkin(12));
t176 = sin(pkin(6));
t229 = qJD(1) * t176;
t215 = t178 * t229;
t174 = sin(pkin(12));
t183 = sin(qJ(3));
t238 = t174 * t183;
t91 = t186 * (t166 * t175 + t179 * t215) - t229 * t238;
t233 = t178 * t179;
t193 = (t174 * t186 + t183 * t233) * t176;
t237 = t175 * t183;
t93 = qJD(1) * t193 + t166 * t237;
t44 = t177 * t93 - t239 * t91;
t291 = t108 - t44;
t234 = t177 * t185;
t200 = pkin(5) * t182 - pkin(10) * t234;
t196 = t200 * qJD(4);
t290 = t196 + t291;
t128 = t173 * t135;
t220 = pkin(10) * t239;
t235 = t177 * t182;
t45 = t173 * t93 + t234 * t91;
t289 = t45 - t128 - (-pkin(9) * t235 - t220) * qJD(4);
t263 = pkin(10) + qJ(5);
t159 = t263 * t173;
t160 = t263 * t177;
t181 = sin(qJ(6));
t184 = cos(qJ(6));
t118 = -t159 * t181 + t160 * t184;
t153 = t173 * t184 + t177 * t181;
t155 = t203 * qJD(3);
t120 = t166 * t179 - t175 * t215;
t90 = qJD(3) * pkin(9) + t93;
t81 = t182 * t90;
t55 = t120 * t185 - t81;
t41 = t177 * t155 - t173 * t55;
t35 = qJD(3) * t200 + t41;
t42 = t173 * t155 + t177 * t55;
t36 = -qJD(3) * t220 + t42;
t288 = -qJD(5) * t153 - qJD(6) * t118 + t181 * t36 - t184 * t35;
t117 = -t159 * t184 - t160 * t181;
t201 = t173 * t181 - t177 * t184;
t287 = -qJD(5) * t201 + qJD(6) * t117 - t181 * t35 - t184 * t36;
t286 = qJD(4) / 0.2e1;
t194 = (t186 * t233 - t238) * t176;
t236 = t175 * t186;
t285 = t180 * t236 + t194;
t228 = qJD(3) * t182;
t284 = t228 / 0.2e1;
t212 = -Ifges(5,6) * qJD(4) / 0.2e1;
t283 = Ifges(5,5) * t286;
t158 = -pkin(4) * t185 - qJ(5) * t182 - pkin(3);
t149 = t177 * t158;
t104 = -pkin(10) * t235 + t149 + (-pkin(9) * t173 - pkin(5)) * t185;
t122 = pkin(9) * t234 + t173 * t158;
t114 = -pkin(10) * t173 * t182 + t122;
t60 = t104 * t181 + t114 * t184;
t282 = -qJD(6) * t60 + t289 * t181 + t290 * t184;
t59 = t104 * t184 - t114 * t181;
t281 = qJD(6) * t59 + t290 * t181 - t289 * t184;
t150 = qJD(4) * t177 - t173 * t228;
t151 = t173 * qJD(4) + t177 * t228;
t209 = t184 * t150 - t151 * t181;
t138 = t201 * qJD(6);
t223 = qJD(4) * t185;
t85 = (qJD(1) * t194 + t166 * t236) * qJD(3);
t245 = t120 * t223 + t185 * t85;
t32 = (qJD(5) - t81) * qJD(4) + t245;
t65 = (t93 + t135) * qJD(3);
t14 = t173 * t65 + t177 * t32;
t222 = qJD(3) * qJD(4);
t210 = t185 * t222;
t208 = t173 * t210;
t10 = -pkin(10) * t208 + t14;
t226 = qJD(3) * t185;
t241 = t120 * t182;
t56 = t185 * t90 + t241;
t53 = qJD(4) * qJ(5) + t56;
t75 = qJD(3) * t158 - t91;
t25 = -t173 * t53 + t177 * t75;
t21 = -pkin(5) * t226 - pkin(10) * t151 + t25;
t26 = t173 * t75 + t177 * t53;
t22 = pkin(10) * t150 + t26;
t5 = -t181 * t22 + t184 * t21;
t13 = -t173 * t32 + t177 * t65;
t9 = qJD(3) * t196 + t13;
t1 = qJD(6) * t5 + t10 * t184 + t181 * t9;
t6 = t181 * t21 + t184 * t22;
t2 = -qJD(6) * t6 - t10 * t181 + t184 * t9;
t197 = t201 * t185;
t191 = qJD(4) * t197;
t69 = -qJD(3) * t191 + qJD(6) * t209;
t103 = t150 * t181 + t151 * t184;
t198 = t153 * t185;
t192 = qJD(4) * t198;
t70 = -qJD(3) * t192 - qJD(6) * t103;
t280 = -t2 * mrSges(7,1) + t1 * mrSges(7,2) - Ifges(7,5) * t69 - Ifges(7,6) * t70;
t172 = Ifges(5,4) * t226;
t257 = Ifges(6,2) * t173;
t260 = Ifges(6,4) * t177;
t204 = -t257 + t260;
t261 = Ifges(6,4) * t173;
t205 = Ifges(6,1) * t177 - t261;
t262 = mrSges(6,2) * t177;
t206 = mrSges(6,1) * t173 + t262;
t267 = t177 / 0.2e1;
t268 = -t173 / 0.2e1;
t51 = -qJD(4) * pkin(4) + qJD(5) - t55;
t89 = -qJD(3) * pkin(3) - t91;
t279 = -(t173 * t26 + t177 * t25) * mrSges(6,3) + t51 * t206 + t89 * mrSges(5,2) + Ifges(5,1) * t284 + t172 / 0.2e1 + t283 + t150 * t204 / 0.2e1 + t151 * t205 / 0.2e1 + (t151 * Ifges(6,4) + t150 * Ifges(6,2) - Ifges(6,6) * t226) * t268 + (t151 * Ifges(6,1) + t150 * Ifges(6,4) - Ifges(6,5) * t226) * t267 - t55 * mrSges(5,3);
t278 = t69 / 0.2e1;
t277 = t70 / 0.2e1;
t276 = -t209 / 0.2e1;
t275 = t209 / 0.2e1;
t274 = -t103 / 0.2e1;
t273 = t103 / 0.2e1;
t131 = t153 * t182;
t272 = -t131 / 0.2e1;
t132 = t201 * t182;
t271 = -t132 / 0.2e1;
t169 = qJD(6) - t226;
t270 = -t169 / 0.2e1;
t269 = t169 / 0.2e1;
t266 = pkin(5) * t173;
t248 = t182 * t85;
t34 = qJD(4) * t56 + t248;
t106 = t180 * t237 + t193;
t136 = -t175 * t176 * t178 + t179 * t180;
t73 = t106 * t182 - t136 * t185;
t265 = t34 * t73;
t259 = Ifges(7,4) * t103;
t258 = Ifges(6,5) * t177;
t256 = Ifges(6,6) * t173;
t86 = t93 * qJD(3);
t253 = t285 * t86;
t140 = -t185 * t179 + t182 * t237;
t252 = t140 * t34;
t247 = t182 * t91;
t246 = t186 * t86;
t127 = mrSges(6,1) * t208 + t210 * t262;
t31 = -t70 * mrSges(7,1) + t69 * mrSges(7,2);
t244 = t127 + t31;
t125 = qJD(3) * t198;
t139 = t153 * qJD(6);
t232 = -t125 + t139;
t126 = qJD(3) * t197;
t231 = -t126 + t138;
t230 = qJD(4) * mrSges(5,1) + mrSges(6,1) * t150 - mrSges(6,2) * t151 - mrSges(5,3) * t228;
t227 = qJD(3) * t183;
t225 = qJD(3) * t186;
t221 = pkin(9) * t182 * t34;
t54 = -mrSges(7,1) * t209 + mrSges(7,2) * t103;
t218 = t54 - t230;
t216 = pkin(9) + t266;
t214 = t175 * t227;
t213 = t175 * t225;
t211 = t182 * t222;
t74 = t106 * t185 + t136 * t182;
t39 = -t173 * t74 - t177 * t285;
t40 = -t173 * t285 + t177 * t74;
t11 = -t181 * t40 + t184 * t39;
t12 = t181 * t39 + t184 * t40;
t141 = t179 * t182 + t185 * t237;
t110 = -t141 * t173 - t177 * t236;
t111 = t141 * t177 - t173 * t236;
t63 = t110 * t184 - t111 * t181;
t64 = t110 * t181 + t111 * t184;
t49 = t241 + (qJD(3) * t266 + t90) * t185;
t189 = t25 * mrSges(6,1) + t5 * mrSges(7,1) + t89 * mrSges(5,1) + t212 - (Ifges(5,4) * t182 + t185 * Ifges(5,2)) * qJD(3) / 0.2e1 + t169 * Ifges(7,3) + t103 * Ifges(7,5) + t209 * Ifges(7,6) - Ifges(6,3) * t226 / 0.2e1 + t151 * Ifges(6,5) + t150 * Ifges(6,6) - t26 * mrSges(6,2) - t56 * mrSges(5,3) - t6 * mrSges(7,2);
t187 = qJD(3) ^ 2;
t171 = -pkin(5) * t177 - pkin(4);
t167 = Ifges(7,3) * t211;
t164 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t226;
t156 = t216 * t182;
t154 = (-mrSges(5,1) * t185 + mrSges(5,2) * t182) * qJD(3);
t147 = t216 * t223;
t146 = (mrSges(5,1) * t182 + mrSges(5,2) * t185) * t222;
t134 = (mrSges(6,1) * t182 - mrSges(6,3) * t234) * t222;
t133 = (-mrSges(6,2) * t182 - mrSges(6,3) * t239) * t222;
t124 = -mrSges(6,1) * t226 - mrSges(6,3) * t151;
t123 = mrSges(6,2) * t226 + mrSges(6,3) * t150;
t121 = -pkin(9) * t239 + t149;
t116 = (Ifges(6,5) * t182 + t185 * t205) * t222;
t115 = (Ifges(6,6) * t182 + t185 * t204) * t222;
t113 = qJD(4) * t141 + t182 * t213;
t112 = -qJD(4) * t140 + t185 * t213;
t109 = -t177 * t219 + t128;
t99 = t106 * qJD(3);
t98 = t285 * qJD(3);
t97 = Ifges(7,4) * t209;
t88 = t138 * t182 - t192;
t87 = -t139 * t182 - t191;
t84 = t112 * t177 + t173 * t214;
t83 = -t112 * t173 + t177 * t214;
t80 = mrSges(7,1) * t169 - mrSges(7,3) * t103;
t79 = -mrSges(7,2) * t169 + mrSges(7,3) * t209;
t58 = -mrSges(7,2) * t211 + mrSges(7,3) * t70;
t57 = mrSges(7,1) * t211 - mrSges(7,3) * t69;
t48 = Ifges(7,1) * t103 + Ifges(7,5) * t169 + t97;
t47 = Ifges(7,2) * t209 + Ifges(7,6) * t169 + t259;
t43 = -pkin(5) * t150 + t51;
t38 = -qJD(4) * t73 + t185 * t98;
t37 = qJD(4) * t74 + t182 * t98;
t33 = -t224 * t90 + t245;
t29 = qJD(4) * t49 + t248;
t28 = t69 * Ifges(7,1) + t70 * Ifges(7,4) + Ifges(7,5) * t211;
t27 = t69 * Ifges(7,4) + t70 * Ifges(7,2) + Ifges(7,6) * t211;
t24 = t173 * t99 + t177 * t38;
t23 = -t173 * t38 + t177 * t99;
t18 = -qJD(6) * t64 - t181 * t84 + t184 * t83;
t17 = qJD(6) * t63 + t181 * t83 + t184 * t84;
t4 = -qJD(6) * t12 - t181 * t24 + t184 * t23;
t3 = qJD(6) * t11 + t181 * t23 + t184 * t24;
t7 = [-t285 * t146 + t11 * t57 + t12 * t58 + t24 * t123 + t23 * t124 + t40 * t133 + t39 * t134 + t99 * t154 + t38 * t164 + t3 * t79 + t4 * t80 + t244 * t73 + t218 * t37 + (-t99 * mrSges(4,1) - t98 * mrSges(4,2) + (-t182 * t74 + t185 * t73) * qJD(4) * mrSges(5,3)) * qJD(3) + m(7) * (t1 * t12 + t11 * t2 + t29 * t73 + t3 * t6 + t37 * t43 + t4 * t5) + m(6) * (t13 * t39 + t14 * t40 + t23 * t25 + t24 * t26 + t37 * t51 + t265) + m(4) * (t106 * t85 - t91 * t99 + t93 * t98 - t253) + m(5) * (t33 * t74 - t37 * t55 + t38 * t56 + t89 * t99 - t253 + t265); -t141 * mrSges(5,3) * t211 + t110 * t134 + t111 * t133 + t112 * t164 + t84 * t123 + t83 * t124 + t17 * t79 + t18 * t80 + t63 * t57 + t64 * t58 + (mrSges(5,3) * t210 + t244) * t140 + t218 * t113 + m(7) * (t1 * t64 + t113 * t43 + t140 * t29 + t17 * t6 + t18 * t5 + t2 * t63) + m(6) * (t110 * t13 + t111 * t14 + t113 * t51 + t25 * t83 + t26 * t84 + t252) + m(5) * (t112 * t56 - t113 * t55 + t141 * t33 + t252) + ((-mrSges(4,2) * t187 - t146) * t186 + (-mrSges(4,1) * t187 + qJD(3) * t154) * t183 + m(4) * (t183 * t85 + t225 * t93 - t227 * t91 - t246) + m(5) * (t227 * t89 - t246)) * t175; ((t189 + t212 + (-m(5) * t56 - t164) * pkin(9)) * t182 + (t283 + (-m(5) * t55 + m(6) * t51 - t230) * pkin(9) + t279) * t185) * qJD(4) + t281 * t79 + (t1 * t60 + t156 * t29 + t2 * t59 + t281 * t6 + t282 * t5 + (t147 - t247) * t43) * m(7) + t282 * t80 + (-t1 * t131 + t132 * t2 - t5 * t87 + t6 * t88) * mrSges(7,3) + (-t91 * t164 + t33 * mrSges(5,3) - t13 * mrSges(6,1) + t14 * mrSges(6,2) - t86 * mrSges(5,1) - t167 / 0.2e1 + t280) * t185 - m(5) * (t89 * t93 + (-t182 * t55 + t185 * t56) * t91) + (-Ifges(7,4) * t132 - Ifges(7,2) * t131) * t277 + (-Ifges(7,1) * t132 - Ifges(7,4) * t131) * t278 + t29 * (mrSges(7,1) * t131 - mrSges(7,2) * t132) + (t109 - t45) * t123 + t28 * t271 + t27 * t272 + (Ifges(7,1) * t87 + Ifges(7,4) * t88) * t273 + (Ifges(7,4) * t87 + Ifges(7,2) * t88) * t275 + (Ifges(7,5) * t87 + Ifges(7,6) * t88) * t269 + t291 * t124 - t93 * t154 + t156 * t31 - pkin(3) * t146 + t147 * t54 + t122 * t133 + t121 * t134 + t87 * t48 / 0.2e1 + t88 * t47 / 0.2e1 + t43 * (-mrSges(7,1) * t88 + mrSges(7,2) * t87) + m(6) * (t108 * t25 + t109 * t26 + t121 * t13 + t122 * t14 + t221) + m(5) * (pkin(9) * t185 * t33 - pkin(3) * t86 + t221) - m(6) * (t247 * t51 + t25 * t44 + t26 * t45) + (t116 * t267 + t115 * t268 + pkin(9) * t127 + t86 * mrSges(5,2) + (mrSges(5,3) + t206) * t34 + (-t13 * t177 - t14 * t173) * mrSges(6,3) - t218 * t91) * t182 + (t93 * mrSges(4,1) + t91 * mrSges(4,2) + (Ifges(7,5) * t271 + Ifges(7,6) * t272 + (-0.3e1 / 0.2e1 * Ifges(5,4) + t258 / 0.2e1 - t256 / 0.2e1) * t182) * t224 + ((-0.3e1 / 0.2e1 * t258 + 0.3e1 / 0.2e1 * t256 + 0.3e1 / 0.2e1 * Ifges(5,4)) * t185 + (Ifges(6,1) * t177 ^ 2 / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(5,1) - 0.3e1 / 0.2e1 * Ifges(5,2) - 0.3e1 / 0.2e1 * Ifges(6,3) - Ifges(7,3) / 0.2e1 + (-t260 + t257 / 0.2e1) * t173) * t182) * t223) * qJD(3) + t59 * t57 + t60 * t58 - t85 * mrSges(4,2) - t86 * mrSges(4,1); (-t25 * t41 - t26 * t42 - t51 * t56 - pkin(4) * t34 + (-t173 * t25 + t177 * t26) * qJD(5) + (-t13 * t173 + t14 * t177) * qJ(5)) * m(6) + t287 * t79 + t288 * t80 + (t1 * t118 + t117 * t2 + t171 * t29 + t287 * t6 + t288 * t5 - t43 * t49) * m(7) + (qJ(5) * t133 + qJD(5) * t123 + t14 * mrSges(6,3) - t34 * mrSges(6,1) + t115 / 0.2e1) * t177 + ((-t172 / 0.2e1 + (-t256 + t258) * t226 / 0.2e1 + ((Ifges(6,1) * t173 + t260) * t267 + (Ifges(6,2) * t177 + t261) * t268 + Ifges(5,5) / 0.2e1) * qJD(4) - t279) * t185 + (Ifges(5,4) * t284 - t189 + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t226 + t212 + (Ifges(6,5) * t173 + Ifges(7,5) * t153 + Ifges(6,6) * t177 - Ifges(7,6) * t201) * t286) * t182) * qJD(3) + (-Ifges(7,1) * t138 - Ifges(7,4) * t139) * t273 + (-Ifges(7,4) * t138 - Ifges(7,2) * t139) * t275 + (-Ifges(7,5) * t138 - Ifges(7,6) * t139) * t269 + (-qJ(5) * t134 - qJD(5) * t124 - t13 * mrSges(6,3) + t34 * mrSges(6,2) + t116 / 0.2e1) * t173 + (Ifges(7,4) * t153 - Ifges(7,2) * t201) * t277 + (Ifges(7,1) * t153 - Ifges(7,4) * t201) * t278 + t29 * (mrSges(7,1) * t201 + mrSges(7,2) * t153) + (-t1 * t201 - t153 * t2 + t231 * t5 - t232 * t6) * mrSges(7,3) - t201 * t27 / 0.2e1 + (-t139 / 0.2e1 + t125 / 0.2e1) * t47 + (-t138 / 0.2e1 + t126 / 0.2e1) * t48 + (-Ifges(7,1) * t126 - Ifges(7,4) * t125) * t274 + (-Ifges(7,4) * t126 - Ifges(7,2) * t125) * t276 + (-Ifges(7,5) * t126 - Ifges(7,6) * t125) * t270 + t171 * t31 - t55 * t164 + t153 * t28 / 0.2e1 - t42 * t123 - t41 * t124 - pkin(4) * t127 + t117 * t57 + t118 * t58 + t230 * t56 + (mrSges(7,1) * t232 - mrSges(7,2) * t231) * t43 - t33 * mrSges(5,2) - t34 * mrSges(5,1) - t49 * t54; t103 * t80 - t209 * t79 - t150 * t123 + t151 * t124 + t244 + (t103 * t5 - t209 * t6 + t29) * m(7) + (-t150 * t26 + t151 * t25 + t34) * m(6); t167 - t43 * (mrSges(7,1) * t103 + mrSges(7,2) * t209) + (Ifges(7,1) * t209 - t259) * t274 + t47 * t273 + (Ifges(7,5) * t209 - Ifges(7,6) * t103) * t270 - t5 * t79 + t6 * t80 + (t103 * t6 + t209 * t5) * mrSges(7,3) + (-Ifges(7,2) * t103 + t48 + t97) * t276 - t280;];
tauc  = t7(:);
