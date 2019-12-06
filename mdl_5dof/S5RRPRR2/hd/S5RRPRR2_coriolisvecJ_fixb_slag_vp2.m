% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR2_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:27:04
% EndTime: 2019-12-05 18:27:23
% DurationCPUTime: 6.12s
% Computational Cost: add. (7049->387), mult. (18848->545), div. (0->0), fcn. (14135->8), ass. (0->186)
t189 = sin(qJ(5));
t192 = cos(qJ(5));
t190 = sin(qJ(4));
t193 = cos(qJ(4));
t206 = qJD(4) * t193;
t207 = qJD(4) * t190;
t187 = sin(pkin(9));
t188 = cos(pkin(9));
t191 = sin(qJ(2));
t194 = cos(qJ(2));
t164 = -t187 * t191 + t188 * t194;
t156 = t164 * qJD(2);
t144 = qJD(1) * t156;
t224 = -qJ(3) - pkin(6);
t200 = qJD(2) * t224;
t151 = qJD(3) * t194 + t191 * t200;
t133 = t151 * qJD(1);
t152 = -t191 * qJD(3) + t194 * t200;
t134 = t152 * qJD(1);
t91 = -t133 * t187 + t188 * t134;
t74 = -pkin(7) * t144 + t91;
t165 = t187 * t194 + t188 * t191;
t155 = t165 * qJD(2);
t143 = qJD(1) * t155;
t92 = t188 * t133 + t187 * t134;
t75 = -pkin(7) * t143 + t92;
t175 = t224 * t194;
t170 = qJD(1) * t175;
t157 = t187 * t170;
t174 = t224 * t191;
t169 = qJD(1) * t174;
t163 = qJD(2) * pkin(2) + t169;
t116 = t188 * t163 + t157;
t209 = qJD(1) * t194;
t210 = qJD(1) * t191;
t154 = -t187 * t209 - t188 * t210;
t225 = pkin(7) * t154;
t85 = qJD(2) * pkin(3) + t116 + t225;
t212 = t188 * t170;
t117 = t187 * t163 - t212;
t153 = t164 * qJD(1);
t226 = pkin(7) * t153;
t90 = t117 + t226;
t22 = t190 * t74 + t193 * t75 + t85 * t206 - t207 * t90;
t114 = t153 * t190 - t154 * t193;
t59 = -qJD(4) * t114 - t143 * t193 - t144 * t190;
t6 = pkin(8) * t59 + t22;
t41 = t190 * t85 + t193 * t90;
t23 = -qJD(4) * t41 - t190 * t75 + t193 * t74;
t197 = t193 * t153 + t154 * t190;
t58 = qJD(4) * t197 - t143 * t190 + t144 * t193;
t7 = -pkin(8) * t58 + t23;
t260 = pkin(8) * t197;
t35 = t41 + t260;
t220 = t189 * t35;
t186 = qJD(2) + qJD(4);
t267 = pkin(8) * t114;
t40 = -t190 * t90 + t193 * t85;
t34 = t40 - t267;
t33 = pkin(4) * t186 + t34;
t8 = t192 * t33 - t220;
t2 = qJD(5) * t8 + t189 * t7 + t192 * t6;
t264 = -t114 * t189 + t192 * t197;
t20 = qJD(5) * t264 + t189 * t59 + t192 * t58;
t65 = t114 * t192 + t189 * t197;
t21 = -qJD(5) * t65 - t189 * t58 + t192 * t59;
t241 = t65 / 0.2e1;
t219 = t192 * t35;
t9 = t189 * t33 + t219;
t3 = -qJD(5) * t9 - t189 * t6 + t192 * t7;
t185 = qJD(5) + t186;
t228 = Ifges(6,4) * t65;
t30 = Ifges(6,2) * t264 + Ifges(6,6) * t185 + t228;
t57 = Ifges(6,4) * t264;
t31 = Ifges(6,1) * t65 + Ifges(6,5) * t185 + t57;
t282 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t20 + Ifges(6,6) * t21 + t30 * t241 - (-Ifges(6,2) * t65 + t31 + t57) * t264 / 0.2e1;
t108 = Ifges(5,4) * t197;
t181 = -pkin(2) * t194 - pkin(1);
t211 = qJD(1) * t181;
t171 = qJD(3) + t211;
t122 = -t153 * pkin(3) + t171;
t221 = Ifges(5,4) * t114;
t53 = Ifges(5,1) * t114 + Ifges(5,5) * t186 + t108;
t280 = t23 * mrSges(5,1) - t22 * mrSges(5,2) + Ifges(5,5) * t58 + Ifges(5,6) * t59 - (Ifges(5,5) * t197 - Ifges(5,6) * t114) * t186 / 0.2e1 - (-Ifges(5,2) * t114 + t108 + t53) * t197 / 0.2e1 - t122 * (mrSges(5,1) * t114 + mrSges(5,2) * t197) - (Ifges(5,1) * t197 - t221) * t114 / 0.2e1 + t282;
t76 = -pkin(4) * t197 + t122;
t279 = -(Ifges(6,5) * t264 - Ifges(6,6) * t65) * t185 / 0.2e1 - t76 * (mrSges(6,1) * t65 + mrSges(6,2) * t264) - (Ifges(6,1) * t264 - t228) * t65 / 0.2e1;
t180 = pkin(2) * t188 + pkin(3);
t227 = pkin(2) * t187;
t149 = t193 * t180 - t190 * t227;
t120 = -t169 * t187 + t212;
t93 = t120 - t226;
t121 = t188 * t169 + t157;
t94 = t121 + t225;
t249 = t149 * qJD(4) - t190 * t93 - t193 * t94;
t150 = t180 * t190 + t193 * t227;
t248 = -t150 * qJD(4) + t190 * t94 - t193 * t93;
t277 = t264 * t8 + t65 * t9;
t271 = t267 + t249;
t270 = t260 + t248;
t269 = t277 * mrSges(6,3) + t279;
t265 = t114 * t41 + t197 * t40;
t52 = Ifges(5,2) * t197 + Ifges(5,6) * t186 + t221;
t261 = t52 / 0.2e1;
t147 = pkin(4) + t149;
t101 = t147 * t192 - t150 * t189;
t258 = qJD(5) * t101 + t270 * t189 + t271 * t192;
t102 = t147 * t189 + t150 * t192;
t257 = -qJD(5) * t102 - t271 * t189 + t270 * t192;
t243 = t264 / 0.2e1;
t240 = pkin(1) * mrSges(3,1);
t239 = pkin(1) * mrSges(3,2);
t237 = t197 / 0.2e1;
t235 = t114 / 0.2e1;
t233 = -t154 / 0.2e1;
t232 = -t155 / 0.2e1;
t231 = t156 / 0.2e1;
t223 = Ifges(3,4) * t191;
t222 = Ifges(4,4) * t154;
t124 = t188 * t174 + t175 * t187;
t106 = -pkin(7) * t165 + t124;
t125 = t187 * t174 - t188 * t175;
t107 = pkin(7) * t164 + t125;
t51 = t190 * t106 + t193 * t107;
t218 = Ifges(3,5) * qJD(2);
t217 = Ifges(3,6) * qJD(2);
t216 = qJD(2) * mrSges(3,1);
t215 = qJD(2) * mrSges(3,2);
t105 = t188 * t151 + t187 * t152;
t208 = qJD(2) * t191;
t183 = pkin(2) * t210;
t204 = t218 / 0.2e1;
t203 = -t217 / 0.2e1;
t202 = -t59 * mrSges(5,1) + t58 * mrSges(5,2);
t201 = -t21 * mrSges(6,1) + t20 * mrSges(6,2);
t179 = qJD(2) * t183;
t123 = pkin(3) * t143 + t179;
t129 = -pkin(3) * t154 + t183;
t130 = pkin(2) * t208 + pkin(3) * t155;
t199 = t143 * mrSges(4,1) + t144 * mrSges(4,2);
t50 = t193 * t106 - t107 * t190;
t104 = -t151 * t187 + t188 * t152;
t119 = t164 * t190 + t165 * t193;
t38 = -pkin(8) * t119 + t50;
t118 = t164 * t193 - t165 * t190;
t39 = pkin(8) * t118 + t51;
t24 = -t189 * t39 + t192 * t38;
t25 = t189 * t38 + t192 * t39;
t69 = t118 * t192 - t119 * t189;
t70 = t118 * t189 + t119 * t192;
t135 = -t164 * pkin(3) + t181;
t82 = -pkin(7) * t156 + t104;
t83 = -pkin(7) * t155 + t105;
t28 = t106 * t206 - t107 * t207 + t190 * t82 + t193 * t83;
t29 = -qJD(4) * t51 - t190 * t83 + t193 * t82;
t182 = Ifges(3,4) * t209;
t173 = mrSges(3,3) * t209 - t215;
t172 = -mrSges(3,3) * t210 + t216;
t162 = Ifges(3,1) * t210 + t182 + t218;
t161 = t217 + (t194 * Ifges(3,2) + t223) * qJD(1);
t148 = Ifges(4,4) * t153;
t132 = qJD(2) * mrSges(4,1) + mrSges(4,3) * t154;
t131 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t153;
t115 = -mrSges(4,1) * t153 - mrSges(4,2) * t154;
t110 = -Ifges(4,1) * t154 + Ifges(4,5) * qJD(2) + t148;
t109 = Ifges(4,2) * t153 + Ifges(4,6) * qJD(2) - t222;
t96 = mrSges(5,1) * t186 - mrSges(5,3) * t114;
t95 = -mrSges(5,2) * t186 + mrSges(5,3) * t197;
t86 = -t118 * pkin(4) + t135;
t77 = pkin(4) * t114 + t129;
t68 = -qJD(4) * t119 - t155 * t193 - t156 * t190;
t67 = qJD(4) * t118 - t155 * t190 + t156 * t193;
t66 = -mrSges(5,1) * t197 + mrSges(5,2) * t114;
t47 = mrSges(6,1) * t185 - mrSges(6,3) * t65;
t46 = -mrSges(6,2) * t185 + mrSges(6,3) * t264;
t45 = -pkin(4) * t68 + t130;
t42 = -pkin(4) * t59 + t123;
t32 = -mrSges(6,1) * t264 + mrSges(6,2) * t65;
t27 = -qJD(5) * t70 - t189 * t67 + t192 * t68;
t26 = qJD(5) * t69 + t189 * t68 + t192 * t67;
t15 = -pkin(8) * t67 + t29;
t14 = pkin(8) * t68 + t28;
t11 = t192 * t34 - t220;
t10 = -t189 * t34 - t219;
t5 = -qJD(5) * t25 - t14 * t189 + t15 * t192;
t4 = qJD(5) * t24 + t14 * t192 + t15 * t189;
t1 = [t185 * (Ifges(6,5) * t26 + Ifges(6,6) * t27) / 0.2e1 + t186 * (Ifges(5,5) * t67 + Ifges(5,6) * t68) / 0.2e1 + t122 * (-mrSges(5,1) * t68 + mrSges(5,2) * t67) + t123 * (-mrSges(5,1) * t118 + mrSges(5,2) * t119) + t130 * t66 + t105 * t131 + t104 * t132 + t181 * t199 + t86 * t201 + t135 * t202 + (t2 * t69 - t20 * t24 + t21 * t25 - t26 * t8 + t27 * t9 - t3 * t70) * mrSges(6,3) + t171 * (mrSges(4,1) * t155 + mrSges(4,2) * t156) + (t118 * t22 - t119 * t23 - t40 * t67 + t41 * t68 - t50 * t58 + t51 * t59) * mrSges(5,3) + t110 * t231 + t109 * t232 + (t144 * t165 + t156 * t233) * Ifges(4,1) + (t58 * t119 + t235 * t67) * Ifges(5,1) + (t118 * t59 + t237 * t68) * Ifges(5,2) + (t118 * t58 + t59 * t119 + t235 * t68 + t237 * t67) * Ifges(5,4) + (Ifges(4,5) * t231 + Ifges(4,6) * t232 + (t162 / 0.2e1 - pkin(6) * t172 + t204 + (-0.2e1 * t239 + 0.3e1 / 0.2e1 * Ifges(3,4) * t194) * qJD(1)) * t194) * qJD(2) + (-t161 / 0.2e1 - pkin(6) * t173 + t203 + (-0.2e1 * t240 - 0.3e1 / 0.2e1 * t223 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t194) * qJD(1) + (qJD(1) * (-mrSges(4,1) * t164 + mrSges(4,2) * t165) + m(4) * (t171 + t211) + t115) * pkin(2)) * t208 + (t20 * t70 + t241 * t26) * Ifges(6,1) + (t69 * t21 + t243 * t27) * Ifges(6,2) + (t69 * t20 + t21 * t70 + t241 * t27 + t243 * t26) * Ifges(6,4) + m(4) * (t104 * t116 + t105 * t117 + t124 * t91 + t125 * t92) + t28 * t95 + t29 * t96 + t76 * (-mrSges(6,1) * t27 + mrSges(6,2) * t26) + t67 * t53 / 0.2e1 + t42 * (-mrSges(6,1) * t69 + mrSges(6,2) * t70) + t45 * t32 + t4 * t46 + t5 * t47 + t26 * t31 / 0.2e1 + t27 * t30 / 0.2e1 + m(6) * (t2 * t25 + t24 * t3 + t4 * t9 + t42 * t86 + t45 * t76 + t5 * t8) + m(5) * (t122 * t130 + t123 * t135 + t22 * t51 + t23 * t50 + t28 * t41 + t29 * t40) + (-t116 * t156 - t117 * t155 - t124 * t144 - t125 * t143 + t164 * t92 - t165 * t91) * mrSges(4,3) + (-t143 * t165 + t164 * t144 + t153 * t231 - t155 * t233) * Ifges(4,4) + (-t164 * t143 + t153 * t232) * Ifges(4,2) + t68 * t261; t279 + (-t101 * t20 + t102 * t21 + t277) * mrSges(6,3) + m(4) * (t187 * t92 + t188 * t91) * pkin(2) - t171 * (-mrSges(4,1) * t154 + mrSges(4,2) * t153) - qJD(2) * (Ifges(4,5) * t153 + Ifges(4,6) * t154) / 0.2e1 - Ifges(4,6) * t143 + Ifges(4,5) * t144 - t129 * t66 - t121 * t131 - t120 * t132 - m(4) * (t116 * t120 + t117 * t121) + t280 - (Ifges(4,2) * t154 + t110 + t148) * t153 / 0.2e1 + t109 * t233 + t154 * (Ifges(4,1) * t153 + t222) / 0.2e1 + ((t204 - t162 / 0.2e1 - t182 / 0.2e1 + qJD(1) * t239 + (t172 - t216) * pkin(6)) * t194 + (t203 + t161 / 0.2e1 + (t240 + t223 / 0.2e1 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t194) * qJD(1) + (t173 + t215) * pkin(6) + (-m(4) * t171 - t115) * pkin(2)) * t191) * qJD(1) + t257 * t47 + (t101 * t3 + t102 * t2 + t257 * t8 + t258 * t9 - t76 * t77) * m(6) + t258 * t46 + t91 * mrSges(4,1) - t92 * mrSges(4,2) - t77 * t32 + (-t149 * t58 + t150 * t59 + t265) * mrSges(5,3) + t114 * t261 + (t116 * t153 - t117 * t154 + (-t143 * t187 - t144 * t188) * pkin(2)) * mrSges(4,3) + t248 * t96 + t249 * t95 + (-t122 * t129 + t149 * t23 + t150 * t22 + t248 * t40 + t249 * t41) * m(5); t114 * t96 - t197 * t95 - t153 * t131 - t154 * t132 - t264 * t46 + t65 * t47 + t199 + t201 + t202 + (-t264 * t9 + t65 * t8 + t42) * m(6) + (t114 * t40 - t197 * t41 + t123) * m(5) + (-t116 * t154 - t117 * t153 + t179) * m(4); -m(6) * (t10 * t8 + t11 * t9) + t52 * t235 - t40 * t95 + t41 * t96 - t11 * t46 - t10 * t47 + t265 * mrSges(5,3) + (-t114 * t32 + (-t189 * t47 + t192 * t46) * qJD(5) + (t189 * t21 - t192 * t20) * mrSges(6,3) + (-t114 * t76 + t189 * t2 + t192 * t3 + (-t189 * t8 + t192 * t9) * qJD(5)) * m(6)) * pkin(4) + t269 + t280; -t8 * t46 + t9 * t47 + t269 + t282;];
tauc = t1(:);
