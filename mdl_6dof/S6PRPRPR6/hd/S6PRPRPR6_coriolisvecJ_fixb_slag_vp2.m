% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRPR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR6_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:46:52
% EndTime: 2019-03-08 19:46:59
% DurationCPUTime: 4.48s
% Computational Cost: add. (3669->459), mult. (8702->672), div. (0->0), fcn. (5885->10), ass. (0->219)
t276 = -qJD(4) / 0.2e1;
t159 = sin(qJ(6));
t162 = cos(qJ(6));
t155 = sin(pkin(11));
t157 = cos(pkin(11));
t177 = t155 * t159 - t157 * t162;
t163 = cos(qJ(4));
t160 = sin(qJ(4));
t222 = t157 * t160;
t205 = pkin(9) * t222;
t171 = (pkin(5) * t163 + t205) * qJD(2);
t182 = pkin(4) * t163 + qJ(5) * t160;
t137 = t182 * qJD(2);
t165 = -pkin(2) - pkin(8);
t164 = cos(qJ(2));
t156 = sin(pkin(6));
t215 = qJD(1) * t156;
t200 = t164 * t215;
t181 = qJD(3) - t200;
t124 = qJD(2) * t165 + t181;
t158 = cos(pkin(6));
t214 = qJD(1) * t158;
t151 = t160 * t214;
t90 = t124 * t163 - t151;
t43 = t157 * t137 - t155 * t90;
t35 = t171 + t43;
t225 = t155 * t160;
t206 = pkin(9) * t225;
t44 = t155 * t137 + t157 * t90;
t40 = qJD(2) * t206 + t44;
t242 = pkin(9) + qJ(5);
t145 = t242 * t155;
t146 = t242 * t157;
t88 = -t145 * t162 - t146 * t159;
t275 = -qJD(5) * t177 + qJD(6) * t88 - t159 * t35 - t162 * t40;
t135 = t155 * t162 + t157 * t159;
t89 = -t145 * t159 + t146 * t162;
t274 = -qJD(5) * t135 - qJD(6) * t89 + t159 * t40 - t162 * t35;
t111 = qJD(4) * t182 - qJD(5) * t163 + qJD(3);
t100 = t157 * t111;
t208 = qJD(4) * t163;
t196 = t165 * t208;
t161 = sin(qJ(2));
t221 = t160 * t161;
t94 = (-t155 * t221 + t157 * t164) * t215;
t273 = -t155 * t196 + t100 - t94;
t76 = t155 * t111 + t157 * t196;
t95 = (t155 * t164 + t157 * t221) * t215;
t272 = -t95 + t76;
t195 = Ifges(5,6) * t276;
t271 = Ifges(5,5) * t276;
t270 = -qJD(2) / 0.2e1;
t108 = t135 * t163;
t266 = qJD(6) * t177;
t269 = t177 * qJD(2) - qJD(4) * t108 + t160 * t266;
t110 = t177 * t163;
t117 = t135 * qJD(2);
t119 = t135 * qJD(6);
t268 = -qJD(4) * t110 - t119 * t160 - t117;
t267 = t177 * t160;
t211 = qJD(2) * t163;
t131 = qJD(4) * t157 - t155 * t211;
t132 = t155 * qJD(4) + t157 * t211;
t189 = t162 * t131 - t132 * t159;
t212 = qJD(2) * t161;
t198 = t156 * t212;
t187 = qJD(1) * t198;
t216 = t124 * t208 + t160 * t187;
t49 = (qJD(5) - t151) * qJD(4) + t216;
t74 = (t111 + t200) * qJD(2);
t16 = -t155 * t49 + t157 * t74;
t14 = qJD(4) * t171 + t16;
t17 = t155 * t74 + t157 * t49;
t188 = qJD(4) * t206;
t15 = qJD(2) * t188 + t17;
t213 = qJD(2) * t160;
t199 = t163 * t214;
t91 = t124 * t160 + t199;
t81 = qJD(4) * qJ(5) + t91;
t142 = pkin(4) * t160 - qJ(5) * t163 + qJ(3);
t201 = t161 * t215;
t98 = qJD(2) * t142 + t201;
t33 = -t155 * t81 + t157 * t98;
t18 = pkin(5) * t213 - pkin(9) * t132 + t33;
t34 = t155 * t98 + t157 * t81;
t22 = pkin(9) * t131 + t34;
t5 = -t159 * t22 + t162 * t18;
t1 = qJD(6) * t5 + t14 * t159 + t15 * t162;
t6 = t159 * t18 + t162 * t22;
t2 = -qJD(6) * t6 + t14 * t162 - t15 * t159;
t169 = qJD(4) * t267;
t36 = qJD(2) * t169 + qJD(6) * t189;
t209 = qJD(4) * t160;
t170 = t135 * t209;
t73 = t131 * t159 + t132 * t162;
t37 = qJD(2) * t170 - qJD(6) * t73;
t265 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t36 + Ifges(7,6) * t37;
t139 = qJD(2) * qJ(3) + t201;
t179 = t155 * t34 + t157 * t33;
t238 = Ifges(6,2) * t155;
t240 = Ifges(6,4) * t157;
t183 = t238 - t240;
t241 = Ifges(6,4) * t155;
t184 = -Ifges(6,1) * t157 + t241;
t185 = -mrSges(6,1) * t155 - mrSges(6,2) * t157;
t249 = -t157 / 0.2e1;
t250 = t155 / 0.2e1;
t77 = -qJD(4) * pkin(4) + qJD(5) - t90;
t264 = t179 * mrSges(6,3) + t77 * t185 + t90 * mrSges(5,3) + t271 + (t163 * Ifges(5,1) - Ifges(5,4) * t160) * t270 + t131 * t183 / 0.2e1 + t132 * t184 / 0.2e1 - t139 * mrSges(5,2) + (t132 * Ifges(6,4) + t131 * Ifges(6,2) + Ifges(6,6) * t213) * t250 + (t132 * Ifges(6,1) + t131 * Ifges(6,4) + Ifges(6,5) * t213) * t249;
t263 = -m(4) / 0.2e1;
t262 = t36 / 0.2e1;
t261 = t37 / 0.2e1;
t260 = -t189 / 0.2e1;
t259 = t189 / 0.2e1;
t258 = -t73 / 0.2e1;
t257 = t73 / 0.2e1;
t207 = qJD(2) * qJD(4);
t256 = (Ifges(6,5) * t163 + t160 * t184) * t207 / 0.2e1;
t255 = m(6) * t77;
t254 = -t108 / 0.2e1;
t253 = -t110 / 0.2e1;
t153 = qJD(6) + t213;
t252 = -t153 / 0.2e1;
t251 = t153 / 0.2e1;
t248 = Ifges(7,4) * t73;
t247 = pkin(5) * t155;
t246 = pkin(9) * t163;
t239 = Ifges(6,5) * t157;
t237 = Ifges(6,6) * t155;
t223 = t156 * t164;
t120 = t158 * t160 + t163 * t223;
t141 = t163 * t187;
t59 = qJD(4) * t91 - t141;
t236 = t120 * t59;
t232 = t163 * t59;
t193 = t160 * t207;
t105 = t185 * t193;
t13 = -t37 * mrSges(7,1) + t36 * mrSges(7,2);
t231 = -t105 - t13;
t230 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t131 + mrSges(6,2) * t132 + mrSges(5,3) * t211;
t226 = t139 * t164;
t224 = t156 * t161;
t220 = t160 * t165;
t219 = t161 * t163;
t103 = t160 * t117;
t218 = t103 + t119;
t104 = qJD(2) * t267;
t217 = t104 + t266;
t97 = t155 * t142 + t157 * t220;
t210 = qJD(2) * t164;
t25 = -mrSges(7,1) * t189 + mrSges(7,2) * t73;
t204 = t25 + t230;
t203 = t165 * t232;
t202 = t160 * t223;
t197 = t156 * t210;
t192 = t163 * t207;
t191 = -t155 * t165 + pkin(5);
t190 = -t165 + t247;
t180 = -t155 * t16 + t157 * t17;
t178 = -t155 * t33 + t157 * t34;
t129 = t157 * t142;
t68 = -t157 * t246 + t160 * t191 + t129;
t80 = -t155 * t246 + t97;
t23 = -t159 * t80 + t162 * t68;
t24 = t159 * t68 + t162 * t80;
t121 = t158 * t163 - t202;
t82 = -t121 * t155 + t157 * t224;
t83 = t121 * t157 + t155 * t224;
t28 = -t159 * t83 + t162 * t82;
t29 = t159 * t82 + t162 * t83;
t133 = (qJD(3) + t200) * qJD(2);
t176 = qJ(3) * t133 + qJD(3) * t139;
t175 = t239 / 0.2e1 - t237 / 0.2e1;
t174 = t133 * t161 + t139 * t210;
t66 = t199 + (-qJD(2) * t247 + t124) * t160;
t168 = t139 * mrSges(5,1) + t33 * mrSges(6,1) + t5 * mrSges(7,1) + t195 + (Ifges(5,4) * t163 - Ifges(5,2) * t160) * t270 + t153 * Ifges(7,3) + t73 * Ifges(7,5) + t189 * Ifges(7,6) + Ifges(6,3) * t213 / 0.2e1 + t132 * Ifges(6,5) + t131 * Ifges(6,6) - t34 * mrSges(6,2) - t6 * mrSges(7,2) - t91 * mrSges(5,3);
t166 = qJD(2) ^ 2;
t154 = -pkin(5) * t157 - pkin(4);
t152 = Ifges(7,3) * t192;
t148 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t213;
t138 = (mrSges(5,1) * t160 + mrSges(5,2) * t163) * qJD(2);
t136 = -qJD(2) * pkin(2) + t181;
t130 = t190 * t163;
t127 = (mrSges(5,1) * t163 - mrSges(5,2) * t160) * t207;
t115 = t190 * t209;
t113 = (mrSges(6,1) * t163 + mrSges(6,3) * t222) * t207;
t112 = (-mrSges(6,2) * t163 + mrSges(6,3) * t225) * t207;
t107 = t135 * t160;
t102 = mrSges(6,1) * t213 - mrSges(6,3) * t132;
t101 = -mrSges(6,2) * t213 + mrSges(6,3) * t131;
t96 = -t155 * t220 + t129;
t86 = (Ifges(6,6) * t163 + t160 * t183) * t207;
t85 = -qJD(4) * t202 + t158 * t208 - t163 * t198;
t84 = -qJD(4) * t120 + t160 * t198;
t67 = Ifges(7,4) * t189;
t60 = t188 + t76;
t58 = -qJD(4) * t151 + t216;
t57 = t163 * t266 + t170;
t55 = -t119 * t163 + t169;
t53 = t155 * t197 + t157 * t84;
t52 = -t155 * t84 + t157 * t197;
t51 = mrSges(7,1) * t153 - mrSges(7,3) * t73;
t50 = -mrSges(7,2) * t153 + mrSges(7,3) * t189;
t46 = -pkin(5) * t131 + t77;
t45 = t100 + (t163 * t191 + t205) * qJD(4);
t41 = qJD(4) * t66 - t141;
t39 = t159 * t94 + t162 * t95;
t38 = -t159 * t95 + t162 * t94;
t27 = -mrSges(7,2) * t192 + mrSges(7,3) * t37;
t26 = mrSges(7,1) * t192 - mrSges(7,3) * t36;
t21 = Ifges(7,1) * t73 + Ifges(7,5) * t153 + t67;
t20 = Ifges(7,2) * t189 + Ifges(7,6) * t153 + t248;
t10 = Ifges(7,1) * t36 + Ifges(7,4) * t37 + Ifges(7,5) * t192;
t9 = Ifges(7,4) * t36 + Ifges(7,2) * t37 + Ifges(7,6) * t192;
t8 = -qJD(6) * t29 - t159 * t53 + t162 * t52;
t7 = qJD(6) * t28 + t159 * t52 + t162 * t53;
t4 = -qJD(6) * t24 - t159 * t60 + t162 * t45;
t3 = qJD(6) * t23 + t159 * t45 + t162 * t60;
t11 = [-t121 * mrSges(5,3) * t192 + t53 * t101 + t52 * t102 + t83 * t112 + t82 * t113 + t84 * t148 + t28 * t26 + t29 * t27 + t7 * t50 + t8 * t51 + t204 * t85 + (-mrSges(5,3) * t193 - t231) * t120 + m(6) * (t16 * t82 + t17 * t83 + t33 * t52 + t34 * t53 + t77 * t85 + t236) + m(7) * (t1 * t29 + t120 * t41 + t2 * t28 + t46 * t85 + t5 * t8 + t6 * t7) + m(5) * (t121 * t58 + t84 * t91 - t85 * t90 + t236) + (t138 * t210 + t161 * t127 + ((-mrSges(3,2) + mrSges(4,3)) * t164 + (-mrSges(3,1) + mrSges(4,2)) * t161) * t166 + m(5) * t174 + m(4) * (t136 * t212 + t174) + 0.2e1 * t164 * t187 * t263) * t156; qJ(3) * t127 + t130 * t13 - t115 * t25 + t97 * t112 + t96 * t113 + t55 * t21 / 0.2e1 + t46 * (-mrSges(7,1) * t57 + mrSges(7,2) * t55) + t57 * t20 / 0.2e1 + t23 * t26 + t24 * t27 + ((t168 + (m(5) * t91 + t148) * t165 + t195) * t163 + (t271 + (-m(5) * t90 + t230 + t255) * t165 + t264) * t160) * qJD(4) + t272 * t101 + (t16 * t96 + t17 * t97 + t272 * t34 + t273 * t33 - t203) * m(6) + t273 * t102 + m(5) * (t58 * t220 + t176 - t203) + t181 * t138 + m(4) * t176 + 0.2e1 * ((t136 * t161 + t226) * t263 - m(5) * (t90 * t219 + t91 * t221 + t226) / 0.2e1 + (m(7) * t46 + t255) * t219 / 0.2e1) * t215 + (-t58 * mrSges(5,3) - t148 * t201 + t152 / 0.2e1 + t133 * mrSges(5,1) - t17 * mrSges(6,2) + t16 * mrSges(6,1) + t265) * t160 + (-t38 + t4) * t51 + (-t39 + t3) * t50 - m(7) * (t38 * t5 + t39 * t6) + t133 * mrSges(4,3) + (-m(4) * pkin(2) * t201 + t181 * mrSges(4,3) + (Ifges(7,5) * t253 + Ifges(7,6) * t254 + (-0.3e1 / 0.2e1 * Ifges(5,4) + t175) * t163) * t208 + ((0.3e1 / 0.2e1 * Ifges(5,4) - 0.3e1 / 0.2e1 * t239 + 0.3e1 / 0.2e1 * t237) * t160 + (-0.3e1 / 0.2e1 * Ifges(5,1) + 0.3e1 / 0.2e1 * Ifges(5,2) + 0.3e1 / 0.2e1 * Ifges(6,3) - Ifges(6,1) * t157 ^ 2 / 0.2e1 + Ifges(7,3) / 0.2e1 + (t240 - t238 / 0.2e1) * t155) * t163) * t209) * qJD(2) + (-t155 * t86 / 0.2e1 + t157 * t256 - t165 * t105 + t133 * mrSges(5,2) + (mrSges(5,3) - t185) * t59 + (-t155 * t17 - t157 * t16) * mrSges(6,3) + t204 * t201) * t163 + t41 * (mrSges(7,1) * t108 - mrSges(7,2) * t110) + (-t1 * t108 + t110 * t2 - t5 * t55 + t57 * t6) * mrSges(7,3) + (-Ifges(7,4) * t110 - Ifges(7,2) * t108) * t261 + (-Ifges(7,1) * t110 - Ifges(7,4) * t108) * t262 + m(7) * (t1 * t24 - t115 * t46 + t130 * t41 + t2 * t23 + t3 * t6 + t4 * t5) + (Ifges(7,5) * t55 + Ifges(7,6) * t57) * t251 + t10 * t253 + t9 * t254 + (Ifges(7,1) * t55 + Ifges(7,4) * t57) * t257 + (Ifges(7,4) * t55 + Ifges(7,2) * t57) * t259; -t166 * mrSges(4,3) - t107 * t26 - t267 * t27 + t269 * t51 + t268 * t50 + t231 * t163 + (t112 * t157 - t113 * t155) * t160 + ((t101 * t157 - t102 * t155 + t148) * t163 + t204 * t160) * qJD(4) + m(5) * (t160 * t58 - t232 + (-t160 * t90 + t163 * t91) * qJD(4)) + m(6) * (-t232 + t180 * t160 + (t160 * t77 + t163 * t178) * qJD(4)) + (-t155 * t101 - m(5) * t139 - m(6) * t179 - t157 * t102 - t138 + (-t139 + t201) * m(4)) * qJD(2) + (-t1 * t267 - t107 * t2 - t163 * t41 + t46 * t209 + t268 * t6 + t269 * t5) * m(7); ((-t168 + Ifges(5,4) * t211 / 0.2e1 + t195 + (Ifges(6,5) * t155 + Ifges(7,5) * t135 + Ifges(6,6) * t157 - Ifges(7,6) * t177) * qJD(4) / 0.2e1) * t163 + ((-Ifges(5,4) / 0.2e1 + t175) * t213 + (Ifges(5,1) / 0.2e1 - Ifges(5,2) / 0.2e1 - Ifges(6,3) / 0.2e1) * t211 + (-Ifges(5,5) / 0.2e1 + (Ifges(6,1) * t155 + t240) * t249 + (Ifges(6,2) * t157 + t241) * t250) * qJD(4) - t264) * t160) * qJD(2) - t177 * t9 / 0.2e1 + (-t1 * t177 - t135 * t2 + t217 * t5 - t218 * t6) * mrSges(7,3) + t41 * (mrSges(7,1) * t177 + mrSges(7,2) * t135) + (Ifges(7,4) * t135 - Ifges(7,2) * t177) * t261 + (Ifges(7,1) * t135 - Ifges(7,4) * t177) * t262 + (-t266 / 0.2e1 - t104 / 0.2e1) * t21 + (-Ifges(7,5) * t266 - Ifges(7,6) * t119) * t251 + (-Ifges(7,1) * t266 - Ifges(7,4) * t119) * t257 + (-Ifges(7,4) * t266 - Ifges(7,2) * t119) * t259 - t44 * t101 - t43 * t102 - pkin(4) * t105 + t88 * t26 + t89 * t27 - t66 * t25 - t58 * mrSges(5,2) - t59 * mrSges(5,1) + (qJ(5) * t112 + qJD(5) * t101 + t17 * mrSges(6,3) + t86 / 0.2e1 - t59 * mrSges(6,1)) * t157 + (-pkin(4) * t59 + qJ(5) * t180 + qJD(5) * t178 - t33 * t43 - t34 * t44 - t77 * t91) * m(6) + (mrSges(7,1) * t218 - mrSges(7,2) * t217) * t46 + t274 * t51 + (t1 * t89 + t154 * t41 + t2 * t88 + t274 * t5 + t275 * t6 - t46 * t66) * m(7) + t275 * t50 + t135 * t10 / 0.2e1 - t90 * t148 - t230 * t91 + t154 * t13 + (t59 * mrSges(6,2) - t16 * mrSges(6,3) - qJ(5) * t113 - qJD(5) * t102 + t256) * t155 + (-t119 / 0.2e1 - t103 / 0.2e1) * t20 + (Ifges(7,5) * t104 + Ifges(7,6) * t103) * t252 + (Ifges(7,1) * t104 + Ifges(7,4) * t103) * t258 + (Ifges(7,4) * t104 + Ifges(7,2) * t103) * t260; -t131 * t101 + t132 * t102 - t189 * t50 + t73 * t51 - t231 + (-t189 * t6 + t5 * t73 + t41) * m(7) + (-t131 * t34 + t132 * t33 + t59) * m(6); t152 - t46 * (mrSges(7,1) * t73 + mrSges(7,2) * t189) + (Ifges(7,1) * t189 - t248) * t258 + t20 * t257 + (Ifges(7,5) * t189 - Ifges(7,6) * t73) * t252 - t5 * t50 + t6 * t51 + (t189 * t5 + t6 * t73) * mrSges(7,3) + (-Ifges(7,2) * t73 + t21 + t67) * t260 + t265;];
tauc  = t11(:);
