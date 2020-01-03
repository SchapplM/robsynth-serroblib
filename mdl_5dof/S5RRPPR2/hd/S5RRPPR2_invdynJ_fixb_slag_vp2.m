% Calculate vector of inverse dynamics joint torques for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:57:21
% EndTime: 2020-01-03 11:57:25
% DurationCPUTime: 1.98s
% Computational Cost: add. (2409->303), mult. (3834->425), div. (0->0), fcn. (2190->14), ass. (0->148)
t149 = sin(qJ(5));
t152 = cos(qJ(5));
t147 = cos(pkin(9));
t193 = qJD(4) * t147;
t199 = t147 * t149;
t146 = sin(pkin(8));
t119 = pkin(2) * t146 + qJ(4);
t198 = t147 * t152;
t145 = sin(pkin(9));
t169 = -pkin(4) * t147 - pkin(7) * t145 - pkin(3);
t148 = cos(pkin(8));
t218 = pkin(2) * t148;
t89 = t169 - t218;
t44 = t119 * t198 + t149 * t89;
t153 = cos(qJ(2));
t150 = sin(qJ(2));
t197 = t148 * t150;
t166 = pkin(1) * (t146 * t153 + t197);
t78 = qJD(1) * t166;
t215 = qJD(1) * pkin(1);
t186 = t150 * t215;
t107 = t146 * t186;
t185 = t153 * t215;
t80 = t148 * t185 - t107;
t226 = -qJD(5) * t44 - t149 * t193 - t152 * t78 + t199 * t80;
t43 = -t119 * t199 + t152 * t89;
t225 = qJD(5) * t43 - t149 * t78 + t152 * t193 - t198 * t80;
t220 = qJD(4) - t80;
t143 = qJD(1) + qJD(2);
t106 = -t143 * t147 + qJD(5);
t95 = pkin(2) * t143 + t185;
t53 = t146 * t95 + t148 * t186;
t47 = qJ(4) * t143 + t53;
t37 = -qJD(3) * t147 + t145 * t47;
t224 = t106 * (-Ifges(6,5) * t149 - Ifges(6,6) * t152) / 0.2e1 + t37 * (mrSges(6,1) * t152 - mrSges(6,2) * t149);
t216 = Ifges(6,4) * t152;
t217 = Ifges(6,4) * t149;
t223 = -t149 * (-Ifges(6,2) * t152 - t217) / 0.2e1 + t152 * (-Ifges(6,1) * t149 - t216) / 0.2e1;
t52 = t148 * t95 - t107;
t176 = qJD(4) - t52;
t25 = t143 * t169 + t176;
t38 = qJD(3) * t145 + t147 * t47;
t8 = -t149 * t38 + t152 * t25;
t222 = qJD(5) * t8;
t141 = t145 ^ 2;
t196 = t147 ^ 2 + t141;
t221 = t196 * mrSges(5,3);
t219 = pkin(1) * t153;
t144 = qJ(1) + qJ(2);
t135 = sin(t144);
t126 = pkin(2) * t135;
t136 = cos(t144);
t127 = pkin(2) * t136;
t140 = qJDD(1) + qJDD(2);
t195 = qJD(2) * t150;
t184 = pkin(1) * t195;
t92 = -qJD(1) * t184 + qJDD(1) * t219;
t75 = pkin(2) * t140 + t92;
t194 = qJD(2) * t153;
t93 = (qJD(1) * t194 + qJDD(1) * t150) * pkin(1);
t36 = t146 * t75 + t148 * t93;
t214 = t145 * t37;
t22 = qJ(4) * t140 + qJD(4) * t143 + t36;
t16 = -qJDD(3) * t147 + t145 * t22;
t213 = t147 * t16;
t17 = qJDD(3) * t145 + t147 * t22;
t212 = t147 * t17;
t162 = (t152 * Ifges(6,1) - t217) * t145;
t211 = t149 * (Ifges(6,5) * t106 + t143 * t162);
t161 = (-t149 * Ifges(6,2) + t216) * t145;
t210 = t152 * (Ifges(6,6) * t106 + t143 * t161);
t209 = t16 * t145;
t134 = pkin(8) + t144;
t123 = sin(t134);
t208 = t123 * t145;
t207 = t123 * t147;
t124 = cos(t134);
t206 = t124 * t145;
t205 = t124 * t147;
t204 = t140 * t145;
t203 = t140 * t147;
t202 = t143 * t145;
t201 = t145 * t149;
t200 = t145 * t152;
t128 = pkin(2) + t219;
t84 = pkin(1) * t197 + t146 * t128;
t192 = qJD(5) * t143;
t191 = qJD(5) * t145;
t190 = qJD(5) * t152;
t104 = qJDD(5) - t203;
t70 = (t140 * t152 - t149 * t192) * t145;
t71 = (-t140 * t149 - t143 * t190) * t145;
t189 = Ifges(6,5) * t70 + Ifges(6,6) * t71 + Ifges(6,3) * t104;
t188 = mrSges(6,3) * t201;
t187 = mrSges(6,3) * t200;
t183 = t124 * pkin(3) + t123 * qJ(4) + t127;
t35 = -t146 * t93 + t148 * t75;
t83 = -t146 * t150 * pkin(1) + t128 * t148;
t9 = t149 * t25 + t152 * t38;
t180 = t9 * mrSges(6,3) * t190;
t82 = -mrSges(5,1) * t203 + mrSges(5,2) * t204;
t179 = g(2) * t124 + g(3) * t123;
t178 = -t149 * t8 + t152 * t9;
t177 = t123 * pkin(3) - qJ(4) * t124 + t126;
t175 = qJDD(4) - t35;
t174 = -mrSges(5,1) * t147 + mrSges(5,2) * t145;
t173 = mrSges(6,1) * t149 + mrSges(6,2) * t152;
t172 = t147 * t38 + t214;
t62 = -mrSges(6,2) * t106 - t143 * t188;
t63 = mrSges(6,1) * t106 - t143 * t187;
t171 = t149 * t63 - t152 * t62;
t170 = pkin(4) * t205 + pkin(7) * t206 + t183;
t81 = t148 * pkin(1) * t194 - t146 * t184;
t51 = t169 - t83;
t76 = qJ(4) + t84;
t21 = t149 * t51 + t198 * t76;
t20 = t152 * t51 - t199 * t76;
t72 = qJD(4) + t81;
t168 = (t16 * t76 + t37 * t72) * t145;
t160 = pkin(4) * t207 + pkin(7) * t208 + t177;
t64 = -t123 * t199 - t124 * t152;
t65 = t123 * t198 - t124 * t149;
t159 = -t135 * mrSges(3,1) - t123 * mrSges(4,1) - mrSges(5,1) * t207 - mrSges(6,1) * t65 - t136 * mrSges(3,2) - t124 * mrSges(4,2) - mrSges(6,2) * t64 - mrSges(6,3) * t208;
t66 = -t123 * t152 + t124 * t199;
t67 = t123 * t149 + t124 * t198;
t156 = -t136 * mrSges(3,1) - t124 * mrSges(4,1) - mrSges(5,1) * t205 - t67 * mrSges(6,1) + t135 * mrSges(3,2) + t66 * mrSges(6,2) - mrSges(6,3) * t206 + (mrSges(4,2) - mrSges(5,3)) * t123;
t24 = -pkin(3) * t140 + t175;
t15 = t140 * t169 + t175;
t3 = t149 * t15 + t152 * t17 + t222;
t4 = -qJD(5) * t9 - t149 * t17 + t15 * t152;
t88 = t173 * t145;
t155 = t24 * t174 + t70 * (-Ifges(6,5) * t147 + t162) / 0.2e1 + t71 * (-Ifges(6,6) * t147 + t161) / 0.2e1 + t104 * (-Ifges(6,3) * t147 + (Ifges(6,5) * t152 - Ifges(6,6) * t149) * t145) / 0.2e1 + t16 * t88 + t92 * mrSges(3,1) - t93 * mrSges(3,2) + t35 * mrSges(4,1) - (Ifges(6,4) * t70 + Ifges(6,2) * t71 + Ifges(6,6) * t104) * t201 / 0.2e1 + (Ifges(5,4) * t145 + Ifges(5,2) * t147) * t203 + (Ifges(5,1) * t145 + Ifges(5,4) * t147) * t204 + (Ifges(6,1) * t70 + Ifges(6,4) * t71 + Ifges(6,5) * t104) * t200 / 0.2e1 - t147 * t189 / 0.2e1 + t4 * (-mrSges(6,1) * t147 - t187) + t3 * (mrSges(6,2) * t147 - t188) + t188 * t222 + mrSges(5,3) * t212 + t224 * t191 + t223 * t141 * t192 - (t210 + t211) * t191 / 0.2e1 + (Ifges(3,3) + Ifges(4,3)) * t140;
t154 = cos(qJ(1));
t151 = sin(qJ(1));
t138 = t154 * pkin(1);
t137 = t151 * pkin(1);
t125 = -pkin(3) - t218;
t85 = t174 * t143;
t79 = qJD(2) * t166;
t77 = -pkin(3) - t83;
t69 = t173 * t202;
t46 = -pkin(3) * t143 + t176;
t40 = -mrSges(6,2) * t104 + mrSges(6,3) * t71;
t39 = mrSges(6,1) * t104 - mrSges(6,3) * t70;
t23 = -mrSges(6,1) * t71 + mrSges(6,2) * t70;
t6 = -qJD(5) * t21 + t152 * t79 - t199 * t72;
t5 = qJD(5) * t20 + t149 * t79 + t198 * t72;
t1 = [(t23 * t76 + t69 * t72 - t180) * t145 + ((-t140 * t150 - t143 * t194) * mrSges(3,2) + (t140 * t153 - t143 * t195) * mrSges(3,1) + (-g(2) * t154 - g(3) * t151 + t150 * t93 + t153 * t92) * m(3)) * pkin(1) + m(6) * (t20 * t4 + t21 * t3 + t5 * t9 + t6 * t8 + t168) + (t151 * mrSges(2,2) - t154 * mrSges(2,1) - m(5) * (t138 + t183) - m(6) * (t138 + t170) - m(4) * (t127 + t138) + mrSges(5,2) * t206 + t156) * g(2) + t155 + m(5) * (t24 * t77 + t46 * t79 + (t17 * t76 + t38 * t72) * t147 + t168) + (t140 * t83 - t143 * t79) * mrSges(4,1) + (-t140 * t84 - t143 * t81 - t36) * mrSges(4,2) + t77 * t82 + t79 * t85 + t5 * t62 + t6 * t63 + t20 * t39 + t21 * t40 + (-m(5) * (t137 + t177) + mrSges(5,2) * t208 - t151 * mrSges(2,1) - t154 * mrSges(2,2) - m(6) * (t137 + t160) - m(4) * (t126 + t137) + t159) * g(3) + (g(3) * t124 + t209 + (t140 * t76 + t143 * t72) * t196) * mrSges(5,3) + m(4) * (t35 * t83 + t36 * t84 - t52 * t79 + t53 * t81) + Ifges(2,3) * qJDD(1); (t124 * mrSges(5,3) + t159) * g(3) + t156 * g(2) + t155 + t125 * t82 - t78 * t85 - t36 * mrSges(4,2) + t43 * t39 + t44 * t40 + (t78 * mrSges(4,1) + t80 * mrSges(4,2) + (mrSges(3,1) * t150 + mrSges(3,2) * t153) * t215 + t220 * t221) * t143 + (t179 * mrSges(5,2) + t16 * mrSges(5,3) + t119 * t23 + t220 * t69 - t180) * t145 + ((mrSges(4,1) * t148 - mrSges(4,2) * t146) * pkin(2) + t119 * t221) * t140 + t225 * t62 + t226 * t63 + (-t160 * g(3) - t170 * g(2) + t3 * t44 + t4 * t43 + (qJD(4) * t37 + t119 * t16) * t145 - t214 * t80 + t225 * t9 + t226 * t8) * m(6) + (-t177 * g(3) - t183 * g(2) + t125 * t24 + (t209 + t212) * t119 - t46 * t78 + t220 * t172) * m(5) + (-t126 * g(3) - t127 * g(2) + (t146 * t36 + t148 * t35) * pkin(2) + t52 * t78 - t53 * t80) * m(4); m(4) * qJDD(3) - t147 * t23 + (-t149 * t39 + t152 * t40 + (-t149 * t62 - t152 * t63) * qJD(5)) * t145 + m(5) * (t145 * t17 - t213) + m(6) * (-t213 + (-t149 * t4 + t152 * t3 + (-t149 * t9 - t152 * t8) * qJD(5)) * t145) + (-m(4) - m(5) - m(6)) * g(1); t149 * t40 + t152 * t39 - t171 * qJD(5) + (t178 * qJD(5) + t149 * t3 + t152 * t4 + t179) * m(6) + (t179 + t24) * m(5) + t82 + (-t145 * t69 + t171 * t147 - m(5) * t172 - m(6) * (t198 * t9 - t199 * t8 + t214) - t143 * t221) * t143; -t3 * mrSges(6,2) + t4 * mrSges(6,1) - t8 * t62 + t9 * t63 + g(1) * t88 - g(2) * (mrSges(6,1) * t64 - mrSges(6,2) * t65) - g(3) * (mrSges(6,1) * t66 + mrSges(6,2) * t67) + (t211 / 0.2e1 + t210 / 0.2e1 - t223 * t202 + t178 * mrSges(6,3) - t224) * t202 + t189;];
tau = t1;
