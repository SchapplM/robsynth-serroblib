% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% Datum: 2018-11-23 15:48
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPPRRR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:48:04
% EndTime: 2018-11-23 15:48:09
% DurationCPUTime: 4.89s
% Computational Cost: add. (8451->452), mult. (20495->628), div. (0->0), fcn. (14975->10), ass. (0->205)
t169 = sin(qJ(5));
t255 = -pkin(9) - pkin(8);
t206 = qJD(5) * t255;
t173 = cos(qJ(4));
t166 = cos(pkin(11));
t210 = qJD(1) * t166;
t164 = sin(pkin(11));
t170 = sin(qJ(4));
t212 = t164 * t170;
t140 = -qJD(1) * t212 + t173 * t210;
t215 = t140 * t169;
t148 = t164 * t173 + t166 * t170;
t141 = t148 * qJD(1);
t113 = pkin(4) * t141 - pkin(8) * t140;
t172 = cos(qJ(5));
t157 = sin(pkin(10)) * pkin(1) + qJ(3);
t153 = t157 * qJD(1);
t161 = t166 * qJD(2);
t122 = t161 + (-pkin(7) * qJD(1) - t153) * t164;
t132 = t164 * qJD(2) + t166 * t153;
t123 = pkin(7) * t210 + t132;
t74 = t122 * t173 - t170 * t123;
t51 = t169 * t113 + t172 * t74;
t286 = pkin(9) * t215 + t169 * t206 - t51;
t242 = pkin(9) * t172;
t50 = t172 * t113 - t169 * t74;
t285 = -pkin(5) * t141 + t140 * t242 + t172 * t206 - t50;
t124 = qJD(4) * t172 - t141 * t169;
t138 = qJD(5) - t140;
t75 = t122 * t170 + t123 * t173;
t69 = qJD(4) * pkin(8) + t75;
t184 = -cos(pkin(10)) * pkin(1) - pkin(3) * t166 - pkin(2);
t139 = t184 * qJD(1) + qJD(3);
t86 = -pkin(4) * t140 - pkin(8) * t141 + t139;
t44 = -t169 * t69 + t172 * t86;
t45 = t169 * t86 + t172 * t69;
t191 = t169 * t45 + t172 * t44;
t195 = Ifges(6,5) * t172 - Ifges(6,6) * t169;
t236 = Ifges(6,4) * t172;
t197 = -Ifges(6,2) * t169 + t236;
t237 = Ifges(6,4) * t169;
t199 = Ifges(6,1) * t172 - t237;
t200 = mrSges(6,1) * t169 + mrSges(6,2) * t172;
t244 = t172 / 0.2e1;
t245 = -t169 / 0.2e1;
t125 = qJD(4) * t169 + t141 * t172;
t250 = t125 / 0.2e1;
t238 = Ifges(6,4) * t125;
t63 = t124 * Ifges(6,2) + t138 * Ifges(6,6) + t238;
t119 = Ifges(6,4) * t124;
t64 = t125 * Ifges(6,1) + t138 * Ifges(6,5) + t119;
t68 = -qJD(4) * pkin(4) - t74;
t284 = t63 * t245 + t64 * t244 + t68 * t200 + t124 * t197 / 0.2e1 + t199 * t250 + t138 * t195 / 0.2e1 - t191 * mrSges(6,3);
t171 = cos(qJ(6));
t168 = sin(qJ(6));
t33 = pkin(9) * t124 + t45;
t226 = t168 * t33;
t32 = -pkin(9) * t125 + t44;
t31 = pkin(5) * t138 + t32;
t10 = t171 * t31 - t226;
t224 = t171 * t33;
t11 = t168 * t31 + t224;
t143 = t148 * qJD(4);
t135 = qJD(1) * t143;
t129 = Ifges(7,3) * t135;
t203 = t171 * t124 - t125 * t168;
t79 = t124 * t168 + t125 * t171;
t243 = Ifges(7,4) * t79;
t133 = qJD(6) + t138;
t249 = -t133 / 0.2e1;
t259 = -t79 / 0.2e1;
t261 = -t203 / 0.2e1;
t73 = Ifges(7,4) * t203;
t37 = Ifges(7,1) * t79 + Ifges(7,5) * t133 + t73;
t53 = -pkin(5) * t124 + t68;
t283 = t129 + (Ifges(7,5) * t203 - Ifges(7,6) * t79) * t249 + (t10 * t203 + t11 * t79) * mrSges(7,3) + (-Ifges(7,2) * t79 + t37 + t73) * t261 - t53 * (mrSges(7,1) * t79 + mrSges(7,2) * t203) + (Ifges(7,1) * t203 - t243) * t259;
t202 = qJD(1) * (t164 ^ 2 + t166 ^ 2);
t282 = mrSges(4,3) * t202;
t137 = Ifges(5,4) * t140;
t267 = t137 / 0.2e1 + t141 * Ifges(5,1) / 0.2e1;
t281 = t139 * mrSges(5,2) + Ifges(5,5) * qJD(4) + t267 + t284;
t147 = -t173 * t166 + t212;
t142 = t147 * qJD(4);
t134 = qJD(1) * t142;
t100 = pkin(4) * t135 + pkin(8) * t134;
t208 = qJD(5) * t172;
t209 = qJD(5) * t169;
t178 = t147 * qJD(3);
t59 = -qJD(1) * t178 + t74 * qJD(4);
t17 = t169 * t100 + t172 * t59 + t86 * t208 - t69 * t209;
t90 = -t125 * qJD(5) + t134 * t169;
t12 = pkin(9) * t90 + t17;
t18 = -t45 * qJD(5) + t172 * t100 - t169 * t59;
t89 = t124 * qJD(5) - t134 * t172;
t9 = pkin(5) * t135 - pkin(9) * t89 + t18;
t2 = t10 * qJD(6) + t12 * t171 + t168 * t9;
t29 = t203 * qJD(6) + t168 * t90 + t171 * t89;
t3 = -t11 * qJD(6) - t12 * t168 + t171 * t9;
t30 = -t79 * qJD(6) - t168 * t89 + t171 * t90;
t280 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t29 + Ifges(7,6) * t30;
t263 = t29 / 0.2e1;
t262 = t30 / 0.2e1;
t36 = Ifges(7,2) * t203 + Ifges(7,6) * t133 + t243;
t278 = t36 / 0.2e1;
t247 = t135 / 0.2e1;
t154 = t255 * t169;
t155 = t255 * t172;
t120 = t154 * t171 + t155 * t168;
t277 = t120 * qJD(6) + t285 * t168 + t286 * t171;
t121 = t154 * t168 - t155 * t171;
t276 = -t121 * qJD(6) - t286 * t168 + t285 * t171;
t47 = -mrSges(7,1) * t203 + mrSges(7,2) * t79;
t271 = m(7) * t53 + t47;
t150 = t168 * t172 + t169 * t171;
t269 = qJD(5) + qJD(6);
t116 = t269 * t150;
t96 = t150 * t140;
t221 = t96 - t116;
t185 = t168 * t169 - t171 * t172;
t115 = t269 * t185;
t97 = t185 * t140;
t220 = t97 - t115;
t110 = t185 * t148;
t239 = pkin(7) + t157;
t144 = t239 * t164;
t145 = t239 * t166;
t270 = -t173 * t144 - t145 * t170;
t193 = -t18 * t169 + t17 * t172;
t266 = t18 * mrSges(6,1) - t17 * mrSges(6,2) + Ifges(6,5) * t89 + Ifges(6,6) * t90 + t280;
t265 = Ifges(7,4) * t263 + Ifges(7,2) * t262 + Ifges(7,6) * t247;
t264 = Ifges(7,1) * t263 + Ifges(7,4) * t262 + Ifges(7,5) * t247;
t260 = t203 / 0.2e1;
t258 = t79 / 0.2e1;
t257 = t89 / 0.2e1;
t256 = t90 / 0.2e1;
t252 = -t124 / 0.2e1;
t251 = -t125 / 0.2e1;
t248 = t133 / 0.2e1;
t246 = -t138 / 0.2e1;
t179 = t148 * qJD(3);
t60 = qJD(1) * t179 + t75 * qJD(4);
t235 = t270 * t60;
t230 = t140 * Ifges(5,2);
t227 = t147 * t60;
t112 = -t144 * t170 + t145 * t173;
t101 = t172 * t112;
t102 = pkin(4) * t147 - pkin(8) * t148 + t184;
t55 = t169 * t102 + t101;
t222 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t124 + mrSges(6,2) * t125 + t141 * mrSges(5,3);
t213 = t148 * t169;
t207 = -t47 - t222;
t205 = t135 * mrSges(5,1) - t134 * mrSges(5,2);
t114 = pkin(4) * t143 + pkin(8) * t142;
t81 = t270 * qJD(4) - t178;
t204 = t172 * t114 - t169 * t81;
t54 = t172 * t102 - t112 * t169;
t201 = mrSges(6,1) * t172 - mrSges(6,2) * t169;
t198 = Ifges(6,1) * t169 + t236;
t196 = Ifges(6,2) * t172 + t237;
t194 = Ifges(6,5) * t169 + Ifges(6,6) * t172;
t46 = pkin(5) * t147 - t148 * t242 + t54;
t49 = -pkin(9) * t213 + t55;
t20 = -t168 * t49 + t171 * t46;
t21 = t168 * t46 + t171 * t49;
t192 = -t169 * t17 - t172 * t18;
t190 = t169 * t44 - t172 * t45;
t65 = mrSges(6,1) * t135 - mrSges(6,3) * t89;
t66 = -mrSges(6,2) * t135 + mrSges(6,3) * t90;
t189 = -t169 * t65 + t172 * t66;
t91 = -mrSges(6,2) * t138 + mrSges(6,3) * t124;
t92 = mrSges(6,1) * t138 - mrSges(6,3) * t125;
t188 = -t169 * t92 + t172 * t91;
t187 = -t169 * t91 - t172 * t92;
t186 = -(-t153 * t164 + t161) * t164 + t132 * t166;
t126 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t140;
t183 = -t126 - t188;
t181 = -t142 * t169 + t148 * t208;
t25 = t102 * t208 - t112 * t209 + t169 * t114 + t172 * t81;
t82 = t112 * qJD(4) + t179;
t176 = t11 * mrSges(7,2) + t45 * mrSges(6,2) + Ifges(5,6) * qJD(4) + t141 * Ifges(5,4) + t230 / 0.2e1 - t133 * Ifges(7,3) - t79 * Ifges(7,5) - t203 * Ifges(7,6) - t138 * Ifges(6,3) - t125 * Ifges(6,5) - t124 * Ifges(6,6) - t10 * mrSges(7,1) - t139 * mrSges(5,1) - t44 * mrSges(6,1);
t159 = -pkin(5) * t172 - pkin(4);
t130 = Ifges(6,3) * t135;
t109 = t150 * t148;
t83 = pkin(5) * t213 - t270;
t61 = pkin(5) * t215 + t75;
t58 = mrSges(7,1) * t133 - mrSges(7,3) * t79;
t57 = -mrSges(7,2) * t133 + mrSges(7,3) * t203;
t52 = t181 * pkin(5) + t82;
t48 = -mrSges(6,1) * t90 + mrSges(6,2) * t89;
t42 = t89 * Ifges(6,1) + t90 * Ifges(6,4) + t135 * Ifges(6,5);
t41 = t89 * Ifges(6,4) + t90 * Ifges(6,2) + t135 * Ifges(6,6);
t40 = t269 * t110 + t150 * t142;
t39 = -t116 * t148 + t185 * t142;
t38 = -pkin(5) * t90 + t60;
t26 = -t55 * qJD(5) + t204;
t24 = -mrSges(7,2) * t135 + mrSges(7,3) * t30;
t23 = mrSges(7,1) * t135 - mrSges(7,3) * t29;
t22 = -t181 * pkin(9) + t25;
t19 = t142 * t242 + pkin(5) * t143 + (-t101 + (pkin(9) * t148 - t102) * t169) * qJD(5) + t204;
t14 = t171 * t32 - t226;
t13 = -t168 * t32 - t224;
t8 = -mrSges(7,1) * t30 + mrSges(7,2) * t29;
t5 = -t21 * qJD(6) - t168 * t22 + t171 * t19;
t4 = t20 * qJD(6) + t168 * t19 + t171 * t22;
t1 = [(-t112 * t135 + t134 * t270 + t142 * t74 - t143 * t75) * mrSges(5,3) - t270 * t48 + (t41 * t245 + t42 * t244 - Ifges(5,1) * t134 - Ifges(5,4) * t135 + t199 * t257 + t197 * t256 + t195 * t247 + (mrSges(5,3) + t200) * t60 + t192 * mrSges(6,3) + (-t172 * t63 / 0.2e1 + t64 * t245 + t68 * t201 + t196 * t252 + t198 * t251 + t194 * t246 + t190 * mrSges(6,3)) * qJD(5)) * t148 + (-t59 * mrSges(5,3) + t129 / 0.2e1 + t130 / 0.2e1 + Ifges(5,4) * t134 + (Ifges(7,3) / 0.2e1 + Ifges(5,2) + Ifges(6,3) / 0.2e1) * t135 + t266) * t147 + m(6) * (t17 * t55 + t18 * t54 + t25 * t45 + t26 * t44 + t68 * t82 - t235) + m(5) * (t112 * t59 - t74 * t82 + t75 * t81 - t235) + (-t230 / 0.2e1 - t176) * t143 + t222 * t82 + t184 * t205 + (m(4) * (t157 * t202 + t186) + 0.2e1 * t282) * qJD(3) - (t267 + t281) * t142 + (-Ifges(7,5) * t110 - Ifges(7,6) * t109) * t247 + t38 * (mrSges(7,1) * t109 - mrSges(7,2) * t110) + (-t10 * t39 - t109 * t2 + t11 * t40 + t110 * t3) * mrSges(7,3) + (-Ifges(7,4) * t110 - Ifges(7,2) * t109) * t262 + (-Ifges(7,1) * t110 - Ifges(7,4) * t109) * t263 + t40 * t278 + t55 * t66 + t4 * t57 + t5 * t58 + t54 * t65 + t52 * t47 + t53 * (-mrSges(7,1) * t40 + mrSges(7,2) * t39) + t39 * t37 / 0.2e1 + t21 * t24 + t20 * t23 + m(7) * (t10 * t5 + t11 * t4 + t2 * t21 + t20 * t3 + t38 * t83 + t52 * t53) + (Ifges(7,4) * t39 + Ifges(7,2) * t40) * t260 - t110 * t264 - t109 * t265 + (Ifges(7,5) * t39 + Ifges(7,6) * t40) * t248 + (Ifges(7,1) * t39 + Ifges(7,4) * t40) * t258 + t83 * t8 + t25 * t91 + t26 * t92 + t81 * t126; -t109 * t23 - t110 * t24 + t39 * t57 + t40 * t58 + (-mrSges(5,3) * t134 + t48 + t8) * t147 - t207 * t143 + t183 * t142 + (-t135 * mrSges(5,3) + t187 * qJD(5) + t189) * t148 + m(5) * (-t142 * t75 - t143 * t74 + t148 * t59 + t227) + m(6) * (t143 * t68 + t227 + t190 * t142 + (-t191 * qJD(5) + t193) * t148) + m(7) * (t10 * t40 - t109 * t3 + t11 * t39 - t110 * t2 + t143 * t53 + t147 * t38); -t185 * t23 + t150 * t24 + t169 * t66 + t172 * t65 + t221 * t58 + t220 * t57 + t188 * qJD(5) + t207 * t141 + t183 * t140 - m(5) * (t140 * t75 - t141 * t74) + t205 + (-m(4) * t186 - t282) * qJD(1) + (t221 * t10 + t220 * t11 - t141 * t53 + t150 * t2 - t185 * t3) * m(7) + (-t138 * t190 - t141 * t68 - t192) * m(6); (-Ifges(7,5) * t115 - Ifges(7,6) * t116) * t248 + (-Ifges(7,1) * t115 - Ifges(7,4) * t116) * t258 + (t96 / 0.2e1 - t116 / 0.2e1) * t36 + (Ifges(7,5) * t150 - Ifges(7,6) * t185 + t194) * t247 + (-t220 * t10 + t221 * t11 - t3 * t150 - t185 * t2) * mrSges(7,3) + (Ifges(7,4) * t150 - Ifges(7,2) * t185) * t262 + (Ifges(7,1) * t150 - Ifges(7,4) * t185) * t263 + t38 * (mrSges(7,1) * t185 + mrSges(7,2) * t150) - t185 * t265 + (t97 / 0.2e1 - t115 / 0.2e1) * t37 + (-Ifges(7,4) * t115 - Ifges(7,2) * t116) * t260 + (-t221 * mrSges(7,1) + t220 * mrSges(7,2)) * t53 - t222 * t75 + (-mrSges(5,1) - t201) * t60 + t193 * mrSges(6,3) + (-pkin(4) * t60 - t44 * t50 - t45 * t51 - t68 * t75) * m(6) + (t74 * mrSges(5,3) - t137 / 0.2e1 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t141 - t281) * t140 + (t271 * t169 * pkin(5) + t284) * qJD(5) + (-Ifges(7,4) * t97 - Ifges(7,2) * t96) * t261 + (-Ifges(7,5) * t97 - Ifges(7,6) * t96) * t249 + (-Ifges(7,1) * t97 - Ifges(7,4) * t96) * t259 - t59 * mrSges(5,2) - t61 * t47 - pkin(4) * t48 + (t75 * mrSges(5,3) + t176) * t141 + t277 * t57 + (t276 * t10 + t277 * t11 + t120 * t3 + t121 * t2 + t159 * t38 - t53 * t61) * m(7) + t276 * t58 + t150 * t264 + t41 * t244 + t196 * t256 + t198 * t257 + ((-m(6) * t191 + t187) * qJD(5) + m(6) * t193 + t189) * pkin(8) - t51 * t91 - t50 * t92 + t120 * t23 + t121 * t24 - t74 * t126 - Ifges(5,5) * t134 - Ifges(5,6) * t135 + t159 * t8 + t169 * t42 / 0.2e1; t79 * t278 + t266 + (-Ifges(6,2) * t125 + t119 + t64) * t252 + (t124 * t44 + t125 * t45) * mrSges(6,3) - m(7) * (t10 * t13 + t11 * t14) - t14 * t57 - t13 * t58 + (t168 * t24 + t171 * t23 + m(7) * (t168 * t2 + t171 * t3) - t271 * t125 + (-t168 * t58 + t171 * t57 + m(7) * (-t10 * t168 + t11 * t171)) * qJD(6)) * pkin(5) + t130 + (Ifges(6,5) * t124 - Ifges(6,6) * t125) * t246 + t63 * t250 + (Ifges(6,1) * t124 - t238) * t251 - t44 * t91 + t45 * t92 - t68 * (mrSges(6,1) * t125 + mrSges(6,2) * t124) + t283; -t10 * t57 + t11 * t58 + t36 * t258 + t280 + t283;];
tauc  = t1(:);
