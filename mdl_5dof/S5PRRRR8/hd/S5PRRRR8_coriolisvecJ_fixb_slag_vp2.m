% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% m_mdh [6x1]
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
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRR8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR8_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR8_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:14:49
% EndTime: 2019-12-05 17:15:15
% DurationCPUTime: 4.70s
% Computational Cost: add. (4073->370), mult. (10510->531), div. (0->0), fcn. (7521->10), ass. (0->195)
t155 = sin(qJ(4));
t156 = sin(qJ(3));
t159 = cos(qJ(4));
t160 = cos(qJ(3));
t130 = t155 * t156 - t159 * t160;
t151 = qJD(3) + qJD(4);
t100 = t151 * t130;
t131 = t155 * t160 + t156 * t159;
t101 = t151 * t131;
t157 = sin(qJ(2));
t152 = sin(pkin(5));
t212 = qJD(1) * t152;
t198 = t157 * t212;
t206 = qJD(3) * t156;
t202 = pkin(3) * t206;
t280 = pkin(4) * t101 + pkin(9) * t100 - t198 + t202;
t255 = -pkin(8) - pkin(7);
t199 = qJD(3) * t255;
t134 = t156 * t199;
t140 = t255 * t156;
t141 = t255 * t160;
t174 = t159 * t140 + t141 * t155;
t189 = t160 * t199;
t161 = cos(qJ(2));
t197 = t161 * t212;
t271 = qJD(4) * t174 + t130 * t197 + t159 * t134 + t155 * t189;
t279 = t151 * Ifges(5,6) / 0.2e1;
t149 = -pkin(3) * t160 - pkin(2);
t119 = qJD(2) * t149 - t197;
t228 = t151 * Ifges(5,5);
t278 = t119 * mrSges(5,2) + t228 / 0.2e1;
t126 = t130 * qJD(2);
t127 = t131 * qJD(2);
t229 = t127 * Ifges(5,4);
t277 = t279 + t229 / 0.2e1 - t126 * Ifges(5,2) / 0.2e1;
t154 = sin(qJ(5));
t158 = cos(qJ(5));
t112 = t140 * t155 - t141 * t159;
t95 = pkin(4) * t130 - pkin(9) * t131 + t149;
t53 = -t112 * t154 + t158 * t95;
t275 = qJD(5) * t53 + t280 * t154 + t271 * t158;
t54 = t112 * t158 + t154 * t95;
t274 = -qJD(5) * t54 - t271 * t154 + t280 * t158;
t136 = qJD(2) * pkin(7) + t198;
t190 = pkin(8) * qJD(2) + t136;
t153 = cos(pkin(5));
t211 = qJD(1) * t156;
t194 = t153 * t211;
t105 = t160 * t190 + t194;
t210 = qJD(2) * t152;
t191 = qJD(1) * t210;
t188 = t161 * t191;
t262 = -t105 * qJD(3) - t156 * t188;
t214 = t159 * t105;
t215 = t153 * t160;
t144 = qJD(1) * t215;
t176 = t190 * t156;
t104 = t144 - t176;
t98 = qJD(3) * pkin(3) + t104;
t52 = t155 * t98 + t214;
t213 = qJD(3) * t144 + t160 * t188;
t71 = -qJD(3) * t176 + t213;
t15 = t52 * qJD(4) + t155 * t71 - t159 * t262;
t109 = -t127 * t154 + t151 * t158;
t90 = t100 * qJD(2);
t49 = qJD(5) * t109 - t158 * t90;
t110 = t127 * t158 + t151 * t154;
t50 = -qJD(5) * t110 + t154 * t90;
t16 = -mrSges(6,1) * t50 + mrSges(6,2) * t49;
t272 = m(6) * t15 + t16;
t270 = qJD(4) * t112 - t131 * t197 + t134 * t155 - t159 * t189;
t93 = mrSges(5,1) * t126 + mrSges(5,2) * t127;
t269 = -m(5) * t119 - t93;
t221 = t105 * t155;
t51 = t159 * t98 - t221;
t14 = qJD(4) * t51 + t155 * t262 + t159 * t71;
t45 = pkin(9) * t151 + t52;
t66 = pkin(4) * t126 - pkin(9) * t127 + t119;
t17 = -t154 * t45 + t158 * t66;
t121 = qJD(2) * t202 + t157 * t191;
t91 = t101 * qJD(2);
t30 = pkin(4) * t91 + pkin(9) * t90 + t121;
t2 = qJD(5) * t17 + t14 * t158 + t154 * t30;
t18 = t154 * t66 + t158 * t45;
t3 = -qJD(5) * t18 - t14 * t154 + t158 * t30;
t268 = -t3 * t154 + t158 * t2;
t267 = -t119 * mrSges(5,1) - t17 * mrSges(6,1) + t18 * mrSges(6,2) + t277;
t120 = qJD(5) + t126;
t185 = mrSges(6,1) * t154 + mrSges(6,2) * t158;
t44 = -pkin(4) * t151 - t51;
t171 = t44 * t185;
t180 = Ifges(6,5) * t158 - Ifges(6,6) * t154;
t237 = Ifges(6,4) * t158;
t182 = -Ifges(6,2) * t154 + t237;
t238 = Ifges(6,4) * t154;
t184 = Ifges(6,1) * t158 - t238;
t246 = t158 / 0.2e1;
t247 = -t154 / 0.2e1;
t251 = t110 / 0.2e1;
t239 = Ifges(6,4) * t110;
t42 = Ifges(6,2) * t109 + Ifges(6,6) * t120 + t239;
t106 = Ifges(6,4) * t109;
t43 = t110 * Ifges(6,1) + t120 * Ifges(6,5) + t106;
t266 = (-t154 * t18 - t158 * t17) * mrSges(6,3) + t120 * t180 / 0.2e1 + t184 * t251 + t109 * t182 / 0.2e1 + t171 + t43 * t246 + t42 * t247;
t265 = -Ifges(4,1) / 0.2e1;
t264 = -t42 / 0.2e1;
t207 = qJD(2) * t160;
t263 = -Ifges(4,4) * t207 / 0.2e1;
t261 = -t17 * t154 + t18 * t158;
t260 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t49 + Ifges(6,6) * t50;
t258 = t49 / 0.2e1;
t257 = t50 / 0.2e1;
t256 = t91 / 0.2e1;
t253 = -t109 / 0.2e1;
t252 = -t110 / 0.2e1;
t250 = -t120 / 0.2e1;
t249 = t126 / 0.2e1;
t248 = -t127 / 0.2e1;
t217 = t152 * t157;
t122 = -t156 * t217 + t215;
t123 = t153 * t156 + t160 * t217;
t175 = t159 * t122 - t123 * t155;
t244 = t15 * t175;
t241 = mrSges(5,3) * t126;
t240 = Ifges(4,4) * t156;
t118 = Ifges(5,4) * t126;
t236 = t109 * Ifges(6,6);
t235 = t110 * Ifges(6,5);
t234 = t174 * t15;
t233 = t120 * Ifges(6,3);
t231 = t127 * mrSges(5,3);
t230 = t127 * Ifges(5,1);
t224 = mrSges(5,1) * t151 + mrSges(6,1) * t109 - mrSges(6,2) * t110 - t231;
t223 = Ifges(4,5) * qJD(3);
t222 = Ifges(4,6) * qJD(3);
t220 = t126 * t154;
t219 = t126 * t158;
t218 = t136 * t160;
t216 = t152 * t161;
t209 = qJD(2) * t156;
t208 = qJD(2) * t157;
t205 = qJD(5) * t154;
t204 = qJD(5) * t158;
t203 = qJD(2) * qJD(3);
t196 = t152 * t208;
t195 = t161 * t210;
t193 = t223 / 0.2e1;
t192 = -t222 / 0.2e1;
t94 = pkin(4) * t127 + pkin(9) * t126;
t186 = mrSges(6,1) * t158 - mrSges(6,2) * t154;
t183 = Ifges(6,1) * t154 + t237;
t181 = Ifges(6,2) * t158 + t238;
t179 = Ifges(6,5) * t154 + Ifges(6,6) * t158;
t83 = -t136 * t206 + t213;
t84 = -qJD(3) * t218 + (-qJD(3) * t153 - t195) * t211;
t177 = -t84 * t156 + t83 * t160;
t82 = t122 * t155 + t123 * t159;
t64 = -t154 * t82 - t158 * t216;
t172 = t154 * t216 - t158 * t82;
t113 = -t136 * t156 + t144;
t137 = -qJD(2) * pkin(2) - t197;
t170 = t113 * mrSges(4,3) + t209 * t265 + t263 - t223 / 0.2e1 - t137 * mrSges(4,2);
t114 = t194 + t218;
t169 = t114 * mrSges(4,3) + t222 / 0.2e1 + (Ifges(4,2) * t160 + t240) * qJD(2) / 0.2e1 - t137 * mrSges(4,1);
t23 = mrSges(6,1) * t91 - mrSges(6,3) * t49;
t24 = -mrSges(6,2) * t91 + mrSges(6,3) * t50;
t68 = -mrSges(6,2) * t120 + mrSges(6,3) * t109;
t69 = mrSges(6,1) * t120 - mrSges(6,3) * t110;
t164 = -t154 * t23 + t158 * t24 + m(6) * (-t17 * t204 - t18 * t205 + t268) - t68 * t205 - t69 * t204;
t41 = t233 + t235 + t236;
t79 = -t118 + t228 + t230;
t8 = t49 * Ifges(6,4) + t50 * Ifges(6,2) + t91 * Ifges(6,6);
t9 = t49 * Ifges(6,1) + t50 * Ifges(6,4) + t91 * Ifges(6,5);
t163 = -t14 * mrSges(5,2) + t43 * t219 / 0.2e1 + t179 * t256 + t181 * t257 + t183 * t258 + t8 * t246 - t51 * t241 + t220 * t264 - Ifges(5,5) * t90 - Ifges(5,6) * t91 + t154 * t9 / 0.2e1 + (-t118 + t79) * t249 + (-t229 + t41) * t248 + (-mrSges(5,1) - t186) * t15 + (-Ifges(5,1) * t248 - t180 * t250 - t182 * t253 - t184 * t252 + t171 + t278) * t126 + (Ifges(6,5) * t252 - Ifges(5,2) * t249 + Ifges(6,6) * t253 + Ifges(6,3) * t250 + t267 + t279) * t127 + (-t17 * t219 - t18 * t220 + t268) * mrSges(6,3) + t266 * qJD(5);
t162 = qJD(2) ^ 2;
t139 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t207;
t138 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t209;
t129 = (mrSges(4,1) * t156 + mrSges(4,2) * t160) * t203;
t115 = -mrSges(5,2) * t151 - t241;
t103 = qJD(3) * t122 + t160 * t195;
t102 = -qJD(3) * t123 - t156 * t195;
t87 = Ifges(6,3) * t91;
t77 = pkin(3) * t209 + t94;
t56 = t104 * t159 - t221;
t55 = t104 * t155 + t214;
t40 = mrSges(5,1) * t91 - mrSges(5,2) * t90;
t26 = qJD(4) * t82 - t159 * t102 + t103 * t155;
t25 = qJD(4) * t175 + t102 * t155 + t103 * t159;
t22 = t154 * t94 + t158 * t51;
t21 = -t154 * t51 + t158 * t94;
t20 = t154 * t77 + t158 * t56;
t19 = -t154 * t56 + t158 * t77;
t12 = qJD(5) * t172 - t154 * t25 + t158 * t196;
t11 = qJD(5) * t64 + t154 * t196 + t158 * t25;
t1 = [t102 * t138 + t103 * t139 + t11 * t68 + t25 * t115 + t12 * t69 - t175 * t16 + t64 * t23 - t172 * t24 - t224 * t26 + (t175 * t90 - t82 * t91) * mrSges(5,3) + (-t122 * t160 - t123 * t156) * mrSges(4,3) * t203 + ((-mrSges(3,2) * t162 - t129 - t40) * t161 + (-mrSges(3,1) * t162 + (qJD(2) * (-mrSges(4,1) * t160 + mrSges(4,2) * t156) + t93) * qJD(2)) * t157) * t152 + m(4) * (t102 * t113 + t103 * t114 + t122 * t84 + t123 * t83 + (t137 - t197) * t196) + m(5) * (t14 * t82 - t244 + t25 * t52 - t26 * t51 + (t119 * t208 - t121 * t161) * t152) + m(6) * (t11 * t18 + t12 * t17 - t172 * t2 + t26 * t44 + t3 * t64 - t244); t274 * t69 + t271 * t115 + (t100 * t51 - t101 * t52 - t112 * t91 + t174 * t90) * mrSges(5,3) + (t236 / 0.2e1 + t235 / 0.2e1 + t233 / 0.2e1 + t41 / 0.2e1 - t267 - t277) * t101 + t53 * t23 + t54 * t24 + (t8 * t247 + t9 * t246 - Ifges(5,1) * t90 - Ifges(5,4) * t91 + t121 * mrSges(5,2) + t182 * t257 + t184 * t258 + t180 * t256 + (mrSges(5,3) + t185) * t15 + (-t2 * t154 - t3 * t158) * mrSges(6,3) + (-mrSges(6,3) * t261 + t158 * t264 + t179 * t250 + t181 * t253 + t183 * t252 + t44 * t186 + t43 * t247) * qJD(5)) * t131 + (-t14 * mrSges(5,3) + t87 / 0.2e1 + Ifges(5,4) * t90 + t121 * mrSges(5,1) + (Ifges(5,2) + Ifges(6,3) / 0.2e1) * t91 + t260) * t130 + (0.3e1 / 0.2e1 * t160 ^ 2 - 0.3e1 / 0.2e1 * t156 ^ 2) * Ifges(4,4) * t203 + t275 * t68 - (t79 / 0.2e1 - t118 / 0.2e1 + t230 / 0.2e1 + t266 + t278) * t100 + ((-pkin(7) * t138 - t170 + t193) * t160 + (pkin(3) * t93 - pkin(7) * t139 + t192 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2)) * t207 - t169) * t156) * qJD(3) + m(4) * ((-t113 * t160 - t114 * t156) * qJD(3) + t177) * pkin(7) + t177 * mrSges(4,3) - t174 * t16 - pkin(2) * t129 + t149 * t40 - t270 * t224 + (t17 * t274 + t18 * t275 + t2 * t54 + t270 * t44 + t3 * t53 - t234) * m(6) + (t112 * t14 + t119 * t202 + t121 * t149 - t270 * t51 + t271 * t52 - t234) * m(5) + (-(pkin(2) * t208 + (-t113 * t156 + t114 * t160) * t161) * m(4) + (t138 * t156 - t139 * t160) * t161 + (-m(4) * t137 + t269) * t157) * t212; -t20 * t68 - t19 * t69 - m(5) * (-t51 * t55 + t52 * t56) + t224 * t55 + t52 * t231 + ((t263 + t193 + t170) * t160 + (t192 + (t240 / 0.2e1 + (Ifges(4,2) / 0.2e1 + t265) * t160) * qJD(2) + t269 * pkin(3) + t169) * t156) * qJD(2) + (m(5) * (t14 * t155 - t15 * t159) + (-t155 * t91 + t159 * t90) * mrSges(5,3) + ((-m(5) * t51 + m(6) * t44 - t224) * t155 + (m(5) * t52 + m(6) * t261 - t154 * t69 + t158 * t68 + t115) * t159) * qJD(4)) * pkin(3) + t164 * (pkin(3) * t155 + pkin(9)) + t163 - m(6) * (t17 * t19 + t18 * t20 + t44 * t55) - t83 * mrSges(4,2) + t84 * mrSges(4,1) - t56 * t115 + t114 * t138 - t113 * t139 + t272 * (-pkin(3) * t159 - pkin(4)); -t22 * t68 - t21 * t69 + t164 * pkin(9) + (t224 + t231) * t52 + t163 - m(6) * (t17 * t21 + t18 * t22 + t44 * t52) - t51 * t115 - t272 * pkin(4); t87 - t44 * (mrSges(6,1) * t110 + mrSges(6,2) * t109) + (Ifges(6,1) * t109 - t239) * t252 + t42 * t251 + (Ifges(6,5) * t109 - Ifges(6,6) * t110) * t250 - t17 * t68 + t18 * t69 + (t109 * t17 + t110 * t18) * mrSges(6,3) + (-Ifges(6,2) * t110 + t106 + t43) * t253 + t260;];
tauc = t1(:);
