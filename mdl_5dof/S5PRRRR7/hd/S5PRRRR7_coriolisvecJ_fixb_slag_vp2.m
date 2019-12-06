% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRR7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR7_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR7_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:11:48
% EndTime: 2019-12-05 17:11:58
% DurationCPUTime: 3.72s
% Computational Cost: add. (3685->333), mult. (9339->487), div. (0->0), fcn. (6479->8), ass. (0->170)
t157 = sin(qJ(3));
t214 = -pkin(7) - pkin(6);
t180 = qJD(3) * t214;
t137 = t157 * t180;
t161 = cos(qJ(3));
t138 = t161 * t180;
t141 = t214 * t157;
t142 = t214 * t161;
t156 = sin(qJ(4));
t160 = cos(qJ(4));
t166 = t156 * t157 - t160 * t161;
t186 = qJD(4) * t160;
t187 = qJD(4) * t156;
t162 = cos(qJ(2));
t193 = qJD(1) * t162;
t220 = t160 * t137 + t156 * t138 + t141 * t186 + t142 * t187 + t166 * t193;
t132 = t156 * t161 + t157 * t160;
t165 = t132 * t162;
t97 = t156 * t141 - t160 * t142;
t219 = qJD(1) * t165 - t97 * qJD(4) - t137 * t156 + t160 * t138;
t152 = qJD(3) + qJD(4);
t91 = t152 * t166;
t236 = pkin(8) * t91 + t219;
t92 = t152 * t132;
t235 = -pkin(8) * t92 + t220;
t151 = qJD(5) + t152;
t122 = t166 * qJD(2);
t123 = t132 * qJD(2);
t155 = sin(qJ(5));
t159 = cos(qJ(5));
t173 = -t159 * t122 - t123 * t155;
t158 = sin(qJ(2));
t194 = qJD(1) * t158;
t144 = qJD(2) * pkin(6) + t194;
t175 = pkin(7) * qJD(2) + t144;
t112 = t175 * t157;
t107 = qJD(3) * pkin(3) - t112;
t113 = t175 * t161;
t183 = qJD(1) * qJD(2);
t176 = t162 * t183;
t143 = t161 * t176;
t93 = -qJD(3) * t112 + t143;
t188 = qJD(3) * t161;
t179 = t144 * t188;
t94 = -t179 + (-pkin(7) * t188 - t157 * t193) * qJD(2);
t27 = t107 * t186 - t113 * t187 + t156 * t94 + t160 * t93;
t82 = t92 * qJD(2);
t10 = -pkin(8) * t82 + t27;
t106 = t160 * t113;
t61 = t107 * t156 + t106;
t28 = -qJD(4) * t61 - t156 * t93 + t160 * t94;
t81 = t91 * qJD(2);
t11 = pkin(8) * t81 + t28;
t208 = pkin(8) * t122;
t44 = t61 - t208;
t202 = t155 * t44;
t117 = t123 * pkin(8);
t104 = t156 * t113;
t60 = t160 * t107 - t104;
t43 = -t117 + t60;
t41 = pkin(4) * t152 + t43;
t12 = t159 * t41 - t202;
t2 = qJD(5) * t12 + t10 * t159 + t11 * t155;
t75 = -t122 * t155 + t123 * t159;
t209 = Ifges(6,4) * t75;
t24 = t173 * qJD(5) - t155 * t82 - t159 * t81;
t25 = -qJD(5) * t75 + t155 * t81 - t159 * t82;
t201 = t159 * t44;
t13 = t155 * t41 + t201;
t3 = -qJD(5) * t13 - t10 * t155 + t11 * t159;
t69 = Ifges(6,4) * t173;
t35 = Ifges(6,1) * t75 + Ifges(6,5) * t151 + t69;
t149 = -pkin(3) * t161 - pkin(2);
t128 = qJD(2) * t149 - t193;
t90 = pkin(4) * t122 + t128;
t234 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t24 + Ifges(6,6) * t25 - (Ifges(6,5) * t173 - Ifges(6,6) * t75) * t151 / 0.2e1 - (-Ifges(6,2) * t75 + t35 + t69) * t173 / 0.2e1 - t90 * (mrSges(6,1) * t75 + mrSges(6,2) * t173) - (Ifges(6,1) * t173 - t209) * t75 / 0.2e1;
t207 = t13 * t75;
t232 = -Ifges(4,1) / 0.2e1;
t34 = Ifges(6,2) * t173 + Ifges(6,6) * t151 + t209;
t231 = t34 / 0.2e1;
t190 = qJD(2) * t161;
t230 = -Ifges(4,4) * t190 / 0.2e1;
t96 = t160 * t141 + t142 * t156;
t70 = -pkin(8) * t132 + t96;
t71 = -pkin(8) * t166 + t97;
t37 = t155 * t70 + t159 * t71;
t229 = -qJD(5) * t37 - t235 * t155 + t236 * t159;
t36 = -t155 * t71 + t159 * t70;
t228 = qJD(5) * t36 + t236 * t155 + t235 * t159;
t148 = pkin(3) * t160 + pkin(4);
t184 = qJD(5) * t159;
t185 = qJD(5) * t155;
t195 = t156 * t159;
t62 = t112 * t156 - t106;
t47 = t62 + t208;
t63 = -t160 * t112 - t104;
t48 = -t117 + t63;
t227 = t155 * t48 - t159 * t47 - t148 * t185 + (-t156 * t184 + (-t155 * t160 - t195) * qJD(4)) * pkin(3);
t196 = t155 * t156;
t226 = -t155 * t47 - t159 * t48 + t148 * t184 + (-t156 * t185 + (t159 * t160 - t196) * qJD(4)) * pkin(3);
t39 = -mrSges(6,1) * t173 + mrSges(6,2) * t75;
t84 = mrSges(5,1) * t122 + mrSges(5,2) * t123;
t225 = -t84 - t39;
t221 = t12 * t173;
t115 = t166 * t158;
t153 = t157 ^ 2;
t154 = t161 ^ 2;
t217 = t173 / 0.2e1;
t215 = t75 / 0.2e1;
t212 = -t122 / 0.2e1;
t211 = t123 / 0.2e1;
t206 = mrSges(5,3) * t122;
t205 = Ifges(4,4) * t157;
t204 = t123 * mrSges(5,3);
t203 = t123 * Ifges(5,4);
t200 = Ifges(4,5) * qJD(3);
t199 = Ifges(4,6) * qJD(3);
t145 = -qJD(2) * pkin(2) - t193;
t197 = t145 * t158;
t147 = t158 * t183;
t189 = qJD(3) * t157;
t181 = pkin(3) * t189;
t129 = qJD(2) * t181 + t147;
t192 = qJD(2) * t157;
t191 = qJD(2) * t158;
t182 = qJD(2) * qJD(3);
t178 = t200 / 0.2e1;
t177 = -t199 / 0.2e1;
t174 = (t153 + t154) * t144;
t172 = t199 / 0.2e1 + (Ifges(4,2) * t161 + t205) * qJD(2) / 0.2e1 - t145 * mrSges(4,1);
t171 = t192 * t232 + t230 - t200 / 0.2e1 - t145 * mrSges(4,2);
t102 = -t144 * t189 + t143;
t103 = -t157 * t176 - t179;
t169 = t102 * t161 - t103 * t157;
t114 = t132 * t158;
t64 = -t114 * t159 + t115 * t155;
t65 = -t114 * t155 - t115 * t159;
t88 = -t132 * t155 - t159 * t166;
t89 = t132 * t159 - t155 * t166;
t139 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t192;
t140 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t190;
t168 = t161 * t139 + t157 * t140;
t167 = t157 * t139 - t161 * t140;
t116 = Ifges(5,4) * t122;
t67 = -t122 * Ifges(5,2) + t152 * Ifges(5,6) + t203;
t68 = t123 * Ifges(5,1) + t152 * Ifges(5,5) - t116;
t164 = mrSges(6,3) * t221 - t27 * mrSges(5,2) - t128 * (mrSges(5,1) * t123 - mrSges(5,2) * t122) + t75 * t231 + t28 * mrSges(5,1) - t60 * t206 + t67 * t211 - t123 * (-Ifges(5,1) * t122 - t203) / 0.2e1 - t152 * (-Ifges(5,5) * t122 - Ifges(5,6) * t123) / 0.2e1 - Ifges(5,6) * t82 - Ifges(5,5) * t81 + (-Ifges(5,2) * t123 - t116 + t68) * t122 / 0.2e1 + t234;
t163 = qJD(2) ^ 2;
t127 = (mrSges(4,1) * t157 + mrSges(4,2) * t161) * t182;
t119 = pkin(3) * t195 + t148 * t155;
t118 = -pkin(3) * t196 + t148 * t159;
t110 = pkin(4) * t166 + t149;
t101 = mrSges(5,1) * t152 - t204;
t100 = -mrSges(5,2) * t152 - t206;
t99 = pkin(3) * t192 + pkin(4) * t123;
t76 = pkin(4) * t92 + t181;
t59 = mrSges(6,1) * t151 - mrSges(6,3) * t75;
t58 = -mrSges(6,2) * t151 + mrSges(6,3) * t173;
t55 = pkin(4) * t82 + t129;
t46 = -qJD(2) * t165 + t152 * t115;
t45 = -t122 * t162 - t158 * t92;
t40 = mrSges(5,1) * t82 - mrSges(5,2) * t81;
t30 = -qJD(5) * t89 + t155 * t91 - t159 * t92;
t29 = qJD(5) * t88 - t155 * t92 - t159 * t91;
t15 = t159 * t43 - t202;
t14 = -t155 * t43 - t201;
t9 = -qJD(5) * t65 - t155 * t45 + t159 * t46;
t8 = qJD(5) * t64 + t155 * t46 + t159 * t45;
t6 = -mrSges(6,1) * t25 + mrSges(6,2) * t24;
t1 = [t45 * t100 + t46 * t101 + t8 * t58 + t9 * t59 + (-t24 * t64 + t25 * t65) * mrSges(6,3) + (-t114 * t81 + t115 * t82) * mrSges(5,3) + (-t163 * mrSges(3,2) - t167 * qJD(2) - t127 - t40 - t6) * t162 + (-t163 * mrSges(3,1) - t168 * qJD(3) + (qJD(2) * (-mrSges(4,1) * t161 + mrSges(4,2) * t157) - t225) * qJD(2)) * t158 + m(5) * (-t114 * t28 - t115 * t27 + t128 * t191 - t129 * t162 + t45 * t61 + t46 * t60) + m(6) * (t12 * t9 + t13 * t8 - t162 * t55 + t191 * t90 + t2 * t65 + t3 * t64) + m(4) * (t169 * t158 + (t197 + (t174 - t194) * t162) * qJD(2)); t228 * t58 + (t110 * t55 + t2 * t37 + t3 * t36 + (-t194 + t76) * t90 + t228 * t13 + t229 * t12) * m(6) + t229 * t59 + (t225 * t158 + t167 * t162) * qJD(1) + t219 * t101 + (t129 * t149 + t27 * t97 + t28 * t96 + t220 * t61 + t219 * t60 + (t181 - t194) * t128) * m(5) + t220 * t100 + t169 * mrSges(4,3) + t29 * t35 / 0.2e1 + (t166 * t82 - t212 * t92) * Ifges(5,2) + (-t132 * t81 - t211 * t91) * Ifges(5,1) + (-t132 * t82 + t166 * t81 - t211 * t92 - t212 * t91) * Ifges(5,4) + (-t132 * t28 - t166 * t27 + t60 * t91 - t61 * t92 + t81 * t96 - t82 * t97) * mrSges(5,3) + t30 * t231 + (-pkin(2) * t147 + t169 * pkin(6) - (t162 * t174 + t197) * qJD(1)) * m(4) + (-t12 * t29 + t13 * t30 + t2 * t88 - t24 * t36 + t25 * t37 - t3 * t89) * mrSges(6,3) + t129 * (mrSges(5,1) * t166 + t132 * mrSges(5,2)) + t128 * (t92 * mrSges(5,1) - t91 * mrSges(5,2)) + t152 * (-Ifges(5,5) * t91 - Ifges(5,6) * t92) / 0.2e1 + ((-pkin(6) * t139 - t171 + t178) * t161 + (-pkin(6) * t140 + pkin(3) * t84 + t177 + (-0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1)) * t190 - t172) * t157) * qJD(3) + (0.3e1 / 0.2e1 * t154 - 0.3e1 / 0.2e1 * t153) * Ifges(4,4) * t182 + (t215 * t29 + t89 * t24) * Ifges(6,1) + (t217 * t30 + t88 * t25) * Ifges(6,2) + (t215 * t30 + t217 * t29 + t88 * t24 + t89 * t25) * Ifges(6,4) + t76 * t39 + t55 * (-mrSges(6,1) * t88 + mrSges(6,2) * t89) + t90 * (-mrSges(6,1) * t30 + mrSges(6,2) * t29) - t91 * t68 / 0.2e1 - t92 * t67 / 0.2e1 + t110 * t6 - pkin(2) * t127 + t149 * t40 + t151 * (Ifges(6,5) * t29 + Ifges(6,6) * t30) / 0.2e1; t164 - m(5) * (t60 * t62 + t61 * t63) + ((t178 + t230 + t171) * t161 + (t177 + (t205 / 0.2e1 + (Ifges(4,2) / 0.2e1 + t232) * t161) * qJD(2) + (-m(5) * t128 - t84) * pkin(3) + t172) * t157) * qJD(2) + (t100 * t186 - t101 * t187 + m(5) * (t156 * t27 + t160 * t28 + t186 * t61 - t187 * t60) + (-t156 * t82 + t160 * t81) * mrSges(5,3)) * pkin(3) + t227 * t59 + t226 * t58 + t168 * t144 + (-t118 * t24 + t119 * t25 + t207) * mrSges(6,3) + t61 * t204 - t99 * t39 - t63 * t100 - t62 * t101 - t102 * mrSges(4,2) + t103 * mrSges(4,1) + (t3 * t118 + t2 * t119 + t227 * t12 + t226 * t13 - t90 * t99) * m(6); t164 - t15 * t58 - t14 * t59 + (-t123 * t39 + (-t155 * t59 + t159 * t58) * qJD(5) + (t155 * t25 - t159 * t24) * mrSges(6,3) + (-t12 * t185 - t123 * t90 + t13 * t184 + t155 * t2 + t159 * t3) * m(6)) * pkin(4) - m(6) * (t12 * t14 + t13 * t15) + mrSges(6,3) * t207 + (t101 + t204) * t61 - t60 * t100; t34 * t215 - t12 * t58 + t13 * t59 + (t221 + t207) * mrSges(6,3) + t234;];
tauc = t1(:);
