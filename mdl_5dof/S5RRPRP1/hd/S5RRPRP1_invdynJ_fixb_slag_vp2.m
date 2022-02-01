% Calculate vector of inverse dynamics joint torques for
% S5RRPRP1
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% m [6x1]
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
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:19:34
% EndTime: 2022-01-20 10:19:41
% DurationCPUTime: 2.61s
% Computational Cost: add. (2026->339), mult. (3250->419), div. (0->0), fcn. (1649->12), ass. (0->153)
t224 = Ifges(6,4) / 0.2e1 + Ifges(5,4) / 0.2e1;
t145 = cos(qJ(4));
t202 = mrSges(5,1) + mrSges(6,1);
t223 = -t202 * t145 - mrSges(4,1);
t201 = mrSges(5,2) + mrSges(6,2);
t222 = -mrSges(6,3) - mrSges(5,3);
t221 = Ifges(5,1) + Ifges(6,1);
t219 = Ifges(5,5) + Ifges(6,5);
t218 = Ifges(5,2) + Ifges(6,2);
t217 = Ifges(5,6) + Ifges(6,6);
t136 = qJDD(1) + qJDD(2);
t139 = sin(pkin(8));
t140 = cos(pkin(8));
t143 = sin(qJ(2));
t194 = pkin(1) * qJD(1);
t173 = t143 * t194;
t146 = cos(qJ(2));
t209 = pkin(1) * t146;
t85 = -qJD(2) * t173 + qJDD(1) * t209;
t63 = pkin(2) * t136 + t85;
t180 = qJD(2) * t146;
t86 = (qJD(1) * t180 + qJDD(1) * t143) * pkin(1);
t26 = t139 * t63 + t140 * t86;
t18 = pkin(7) * t136 + t26;
t216 = qJD(3) * qJD(4) + t18;
t215 = (Ifges(5,4) + Ifges(6,4)) * t145;
t214 = m(3) * pkin(1);
t138 = qJ(1) + qJ(2);
t133 = cos(t138);
t121 = pkin(2) * t133;
t208 = pkin(2) * t139;
t207 = pkin(2) * t140;
t206 = pkin(4) * t145;
t129 = pkin(8) + t138;
t116 = sin(t129);
t205 = g(1) * t116;
t117 = cos(t129);
t204 = g(2) * t117;
t142 = sin(qJ(4));
t178 = qJD(4) * t142;
t137 = qJD(1) + qJD(2);
t172 = t146 * t194;
t94 = pkin(2) * t137 + t172;
t45 = t139 * t94 + t140 * t173;
t33 = pkin(7) * t137 + t45;
t5 = t142 * qJDD(3) + t145 * t216 - t33 * t178;
t203 = t145 * t5;
t185 = t137 * t142;
t90 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t185;
t91 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t185;
t200 = t90 + t91;
t184 = t137 * t145;
t92 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t184;
t93 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t184;
t199 = t92 + t93;
t198 = Ifges(5,4) * t142;
t196 = Ifges(6,4) * t142;
t123 = pkin(2) + t209;
t182 = t140 * t143;
t72 = pkin(1) * t182 + t139 * t123;
t65 = pkin(7) + t72;
t193 = t145 * t65;
t183 = t139 * t143;
t69 = (t140 * t146 - t183) * qJD(2) * pkin(1);
t192 = t145 * t69;
t190 = -qJ(5) - t65;
t179 = qJD(3) * t142;
t29 = t145 * t33 + t179;
t189 = qJD(4) * t29;
t188 = qJD(4) * t65;
t118 = pkin(7) + t208;
t181 = -qJ(5) - t118;
t177 = qJD(4) * t145;
t130 = t145 * qJD(5);
t176 = m(4) + m(5) + m(6);
t171 = pkin(4) * t178;
t170 = t117 * pkin(3) + t116 * pkin(7) + t121;
t122 = pkin(3) + t206;
t169 = t137 * t178;
t77 = t136 * t145 - t169;
t78 = t136 * t142 + t137 * t177;
t30 = -t77 * mrSges(6,1) + t78 * mrSges(6,2);
t25 = -t139 * t86 + t140 * t63;
t166 = qJ(5) * t137 + t33;
t102 = t139 * t173;
t44 = t140 * t94 - t102;
t165 = qJD(4) * t190;
t71 = -pkin(1) * t183 + t123 * t140;
t164 = qJD(4) * t181;
t163 = t201 * t204;
t64 = -pkin(3) - t71;
t141 = -qJ(5) - pkin(7);
t160 = -t116 * t141 + t117 * t122 + t121;
t131 = t145 * qJD(3);
t15 = -t166 * t142 + t131;
t13 = qJD(4) * pkin(4) + t15;
t28 = -t142 * t33 + t131;
t159 = t28 * mrSges(5,3) + t13 * mrSges(6,3);
t16 = t166 * t145 + t179;
t158 = t29 * mrSges(5,3) + t16 * mrSges(6,3);
t100 = -mrSges(5,1) * t145 + mrSges(5,2) * t142;
t157 = -mrSges(6,1) * t145 + mrSges(6,2) * t142;
t156 = Ifges(5,2) * t145 + t198;
t155 = Ifges(6,2) * t145 + t196;
t154 = -t13 * t142 + t145 * t16;
t153 = -t142 * t28 + t145 * t29;
t17 = -pkin(3) * t136 - t25;
t152 = pkin(1) * (t139 * t146 + t182);
t67 = qJD(2) * t152;
t132 = sin(t138);
t150 = -t133 * mrSges(3,1) + t132 * mrSges(3,2) + (mrSges(4,2) + t222) * t116 + t223 * t117;
t149 = t133 * mrSges(3,2) + (t176 * pkin(2) + mrSges(3,1)) * t132 - t201 * t116 * t142 + t222 * t117;
t27 = -t122 * t137 + qJD(5) - t44;
t3 = qJ(5) * t77 + t137 * t130 + t5;
t32 = -pkin(3) * t137 - t44;
t59 = Ifges(6,6) * qJD(4) + t155 * t137;
t60 = Ifges(5,6) * qJD(4) + t156 * t137;
t103 = Ifges(6,4) * t184;
t61 = Ifges(6,1) * t185 + Ifges(6,5) * qJD(4) + t103;
t104 = Ifges(5,4) * t184;
t62 = Ifges(5,1) * t185 + Ifges(5,5) * qJD(4) + t104;
t8 = -pkin(4) * t77 + qJDD(5) + t17;
t148 = t3 * t145 * mrSges(6,3) + t85 * mrSges(3,1) + t25 * mrSges(4,1) - t86 * mrSges(3,2) + mrSges(5,3) * t203 + t17 * t100 + t8 * t157 + (-t217 * t142 + t219 * t145) * qJD(4) ^ 2 / 0.2e1 - (t60 + t59) * t178 / 0.2e1 + (t145 * t221 - t196 - t198) * t169 / 0.2e1 + (Ifges(3,3) + Ifges(4,3)) * t136 + (t32 * (mrSges(5,1) * t142 + mrSges(5,2) * t145) + t27 * (mrSges(6,1) * t142 + mrSges(6,2) * t145)) * qJD(4) + (t156 / 0.2e1 + t155 / 0.2e1 + t142 * t224 + t218 * t145 / 0.2e1) * t77 + (t142 * t221 + t215 / 0.2e1 + t145 * t224) * t78 + (t62 + t61 + (-t142 * t218 + t215) * t137) * t177 / 0.2e1 + (t142 * t219 + t145 * t217) * qJDD(4);
t147 = cos(qJ(1));
t144 = sin(qJ(1));
t135 = t147 * pkin(1);
t134 = t145 * qJ(5);
t128 = t145 * qJDD(3);
t119 = -pkin(3) - t207;
t111 = t117 * pkin(7);
t99 = -t122 - t207;
t84 = t118 * t145 + t134;
t83 = t181 * t142;
t76 = t100 * t137;
t75 = t157 * t137;
t68 = t140 * t172 - t102;
t66 = qJD(1) * t152;
t54 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t78;
t53 = qJDD(4) * mrSges(6,1) - mrSges(6,3) * t78;
t52 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t77;
t51 = -qJDD(4) * mrSges(6,2) + mrSges(6,3) * t77;
t50 = -qJD(5) * t142 + t145 * t164;
t49 = t142 * t164 + t130;
t48 = t64 - t206;
t46 = t67 + t171;
t39 = t134 + t193;
t38 = t190 * t142;
t31 = -mrSges(5,1) * t77 + mrSges(5,2) * t78;
t11 = (-qJD(5) - t69) * t142 + t145 * t165;
t10 = t142 * t165 + t130 + t192;
t6 = -t142 * t18 + t128 - t189;
t2 = -t33 * t177 + qJDD(4) * pkin(4) - qJ(5) * t78 + t128 + (-qJD(5) * t137 - t216) * t142;
t1 = [(t136 * t71 - t137 * t67) * mrSges(4,1) + m(4) * (t25 * t71 + t26 * t72 - t44 * t67 + t45 * t69) + m(6) * (t10 * t16 + t11 * t13 + t2 * t38 + t27 * t46 + t3 * t39 + t48 * t8) + (-m(5) * t111 + mrSges(2,2) * t147 + (m(6) * t141 + mrSges(4,2)) * t117 + (mrSges(2,1) + (m(3) + t176) * pkin(1)) * t144 + (m(5) * pkin(3) + m(6) * t122 - t223) * t116 + t149) * g(1) + (t144 * mrSges(2,2) - m(4) * (t121 + t135) - m(6) * (t135 + t160) - m(5) * (t135 + t170) + (-mrSges(2,1) - t214) * t147 + t150) * g(2) + (t143 * t86 + t146 * t85) * t214 + m(5) * (-t28 * t65 * t177 + t17 * t64 + t29 * t192 + t5 * t193 + t32 * t67) + t148 + (t65 * t52 + t69 * t93 + (-t65 * t91 - t159) * qJD(4)) * t145 + ((-t136 * t143 - t137 * t180) * mrSges(3,2) + (-qJD(2) * t137 * t143 + t136 * t146) * mrSges(3,1)) * pkin(1) + (m(5) * (-t29 * t188 - t28 * t69 - t6 * t65) - t69 * t91 - t65 * t54 - t93 * t188 + t163 + (-qJD(4) * t16 - t2) * mrSges(6,3) + (-t6 - t189) * mrSges(5,3)) * t142 + (-t136 * t72 - t137 * t69 - t26) * mrSges(4,2) + t11 * t90 + t10 * t92 + t46 * t75 + t67 * t76 + t48 * t30 + t39 * t51 + t38 * t53 + t64 * t31 + Ifges(2,3) * qJDD(1); (-t6 * mrSges(5,3) - t2 * mrSges(6,3) - t118 * t54 + t200 * t68 + t163 + (pkin(4) * t75 - t118 * t93 - t158) * qJD(4)) * t142 + (t118 * t52 - t199 * t68 + t202 * t205 + (-t118 * t91 - t159) * qJD(4)) * t145 + (t116 * mrSges(4,1) + t117 * mrSges(4,2) + t149) * g(1) + (t66 * mrSges(4,1) + t68 * mrSges(4,2) + (mrSges(3,1) * t143 + mrSges(3,2) * t146) * t194) * t137 + (-t75 - t76) * t66 + (-t136 * t208 - t26) * mrSges(4,2) + t150 * g(2) + t148 + t119 * t31 + t99 * t30 + t83 * t53 + t84 * t51 + t50 * t90 + t49 * t92 + t136 * mrSges(4,1) * t207 + ((t116 * t122 + t117 * t141) * g(1) + t13 * t50 + t16 * t49 + t2 * t83 + t3 * t84 + t8 * t99 - t160 * g(2) - t154 * t68 + (t171 - t66) * t27) * m(6) + ((t139 * t26 + t140 * t25) * pkin(2) + t44 * t66 - t45 * t68 - t121 * g(2)) * m(4) + (t119 * t17 + (pkin(3) * t116 - t111) * g(1) - t170 * g(2) + (-t142 * t6 + t203 + (-t142 * t29 - t145 * t28) * qJD(4)) * t118 - t153 * t68 - t32 * t66) * m(5); m(4) * qJDD(3) + (t53 + t54) * t145 + (t51 + t52) * t142 + (-t200 * t142 + t199 * t145) * qJD(4) + m(5) * (t153 * qJD(4) + t142 * t5 + t145 * t6) + m(6) * (t154 * qJD(4) + t142 * t3 + t145 * t2) - t176 * g(3); t6 * mrSges(5,1) + t2 * mrSges(6,1) - t5 * mrSges(5,2) - t3 * mrSges(6,2) - t15 * t92 - t28 * t93 + t29 * t91 + t219 * t78 + t217 * t77 + (Ifges(5,3) + Ifges(6,3)) * qJDD(4) + (t100 + t157) * g(3) + (t53 + (-g(3) * t145 + t2) * m(6)) * pkin(4) + (-m(6) * (-t13 + t15) + t90) * t16 + ((-t103 / 0.2e1 - t104 / 0.2e1 - t61 / 0.2e1 - t62 / 0.2e1 - t32 * mrSges(5,2) - t27 * mrSges(6,2) + (-Ifges(6,5) / 0.2e1 - Ifges(5,5) / 0.2e1) * qJD(4) + t159) * t145 + (t59 / 0.2e1 + t60 / 0.2e1 - t32 * mrSges(5,1) - t27 * mrSges(6,1) + t224 * t185 + (Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1) * qJD(4) + (-m(6) * t27 - t75) * pkin(4) + (-Ifges(6,1) / 0.2e1 - Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t184 + t158) * t142) * t137 + (g(1) * t117 + g(2) * t116) * (t201 * t145 + (m(6) * pkin(4) + t202) * t142); (t142 * t90 - t145 * t92) * t137 + (-t154 * t137 + t204 - t205 + t8) * m(6) + t30;];
tau = t1;
