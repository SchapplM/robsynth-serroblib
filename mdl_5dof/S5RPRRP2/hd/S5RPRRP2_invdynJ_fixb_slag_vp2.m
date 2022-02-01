% Calculate vector of inverse dynamics joint torques for
% S5RPRRP2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2022-01-23 09:28
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:27:52
% EndTime: 2022-01-23 09:27:56
% DurationCPUTime: 2.43s
% Computational Cost: add. (1896->315), mult. (3271->384), div. (0->0), fcn. (1667->12), ass. (0->137)
t218 = Ifges(6,4) / 0.2e1 + Ifges(5,4) / 0.2e1;
t136 = cos(qJ(4));
t191 = mrSges(5,1) + mrSges(6,1);
t217 = -t136 * t191 - mrSges(4,1);
t190 = mrSges(5,2) + mrSges(6,2);
t216 = -mrSges(5,3) - mrSges(6,3);
t215 = Ifges(5,1) + Ifges(6,1);
t213 = Ifges(6,5) + Ifges(5,5);
t212 = Ifges(5,2) + Ifges(6,2);
t211 = Ifges(6,6) + Ifges(5,6);
t127 = qJDD(1) + qJDD(3);
t134 = sin(qJ(3));
t137 = cos(qJ(3));
t130 = sin(pkin(8));
t198 = pkin(1) * t130;
t131 = cos(pkin(8));
t114 = pkin(1) * t131 + pkin(2);
t91 = t114 * qJD(1);
t208 = qJD(3) * t91 + qJDD(1) * t198;
t165 = qJD(1) * t198;
t209 = -qJD(3) * t165 + t114 * qJDD(1);
t17 = t209 * t134 + t208 * t137;
t13 = pkin(7) * t127 + t17;
t210 = qJD(2) * qJD(4) + t13;
t207 = (Ifges(5,4) + Ifges(6,4)) * t136;
t206 = m(3) * pkin(1);
t133 = sin(qJ(4));
t171 = qJD(2) * t133;
t128 = qJD(1) + qJD(3);
t47 = t134 * t91 + t137 * t165;
t37 = pkin(7) * t128 + t47;
t31 = t136 * t37 + t171;
t174 = t31 * qJD(4);
t121 = t136 * qJDD(2);
t6 = -t13 * t133 + t121 - t174;
t205 = -t6 - t174;
t65 = t114 * t137 - t134 * t198;
t124 = t136 * qJD(2);
t157 = qJ(5) * t128 + t37;
t24 = -t133 * t157 + t124;
t19 = qJD(4) * pkin(4) + t24;
t30 = -t133 * t37 + t124;
t204 = -t30 * mrSges(5,3) - t19 * mrSges(6,3);
t203 = m(5) * pkin(3);
t197 = pkin(4) * t136;
t129 = qJ(1) + pkin(8);
t122 = qJ(3) + t129;
t112 = sin(t122);
t196 = g(1) * t112;
t113 = cos(t122);
t195 = g(2) * t113;
t170 = qJD(4) * t133;
t5 = t133 * qJDD(2) + t210 * t136 - t170 * t37;
t194 = t136 * t5;
t132 = -qJ(5) - pkin(7);
t176 = t128 * t133;
t81 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t176;
t82 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t176;
t189 = t81 + t82;
t175 = t128 * t136;
t83 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t175;
t84 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t175;
t188 = t83 + t84;
t187 = Ifges(5,4) * t133;
t185 = Ifges(6,4) * t133;
t66 = t134 * t114 + t137 * t198;
t63 = pkin(7) + t66;
t182 = -qJ(5) - t63;
t180 = qJD(4) * t63;
t173 = t113 * pkin(3) + t112 * pkin(7);
t119 = cos(t129);
t138 = cos(qJ(1));
t172 = t138 * pkin(1) + pkin(2) * t119;
t169 = qJD(4) * t136;
t123 = t136 * qJD(5);
t168 = m(4) + m(5) + m(6);
t164 = pkin(4) * t170;
t163 = m(3) + t168;
t116 = pkin(3) + t197;
t161 = t128 * t170;
t71 = t127 * t136 - t161;
t72 = t127 * t133 + t128 * t169;
t28 = -t71 * mrSges(6,1) + t72 * mrSges(6,2);
t46 = -t134 * t165 + t137 * t91;
t158 = qJD(4) * t132;
t156 = -t112 * t132 + t113 * t116;
t155 = qJD(4) * t182;
t154 = -t190 * t112 * t133 + t216 * t113;
t153 = t190 * t195;
t62 = -pkin(3) - t65;
t25 = t136 * t157 + t171;
t148 = t31 * mrSges(5,3) + t25 * mrSges(6,3);
t93 = -t136 * mrSges(5,1) + mrSges(5,2) * t133;
t147 = -t136 * mrSges(6,1) + mrSges(6,2) * t133;
t146 = Ifges(5,2) * t136 + t187;
t145 = Ifges(6,2) * t136 + t185;
t144 = -t133 * t19 + t136 * t25;
t143 = -t133 * t30 + t136 * t31;
t18 = -t208 * t134 + t209 * t137;
t142 = (mrSges(4,2) + t216) * t112 + t217 * t113;
t14 = -pkin(3) * t127 - t18;
t57 = t66 * qJD(3);
t27 = -t116 * t128 + qJD(5) - t46;
t3 = qJ(5) * t71 + t123 * t128 + t5;
t36 = -pkin(3) * t128 - t46;
t58 = Ifges(6,6) * qJD(4) + t128 * t145;
t59 = Ifges(5,6) * qJD(4) + t128 * t146;
t99 = Ifges(6,4) * t175;
t60 = Ifges(6,1) * t176 + Ifges(6,5) * qJD(4) + t99;
t100 = Ifges(5,4) * t175;
t61 = Ifges(5,1) * t176 + Ifges(5,5) * qJD(4) + t100;
t8 = -pkin(4) * t71 + qJDD(5) + t14;
t140 = t3 * t136 * mrSges(6,3) + t18 * mrSges(4,1) + mrSges(5,3) * t194 + Ifges(4,3) * t127 + t14 * t93 + t8 * t147 + (-t211 * t133 + t213 * t136) * qJD(4) ^ 2 / 0.2e1 - (t59 + t58) * t170 / 0.2e1 + (t215 * t136 - t185 - t187) * t161 / 0.2e1 + (t27 * (mrSges(6,1) * t133 + mrSges(6,2) * t136) + t36 * (mrSges(5,1) * t133 + mrSges(5,2) * t136)) * qJD(4) + (t146 / 0.2e1 + t145 / 0.2e1 + t133 * t218 + t212 * t136 / 0.2e1) * t71 + (t215 * t133 + t207 / 0.2e1 + t136 * t218) * t72 + (t61 + t60 + (-t212 * t133 + t207) * t128) * t169 / 0.2e1 + (t213 * t133 + t211 * t136) * qJDD(4);
t135 = sin(qJ(1));
t125 = t136 * qJ(5);
t118 = sin(t129);
t107 = t113 * pkin(7);
t94 = pkin(7) * t136 + t125;
t92 = t132 * t133;
t70 = t93 * t128;
t69 = t147 * t128;
t68 = -qJD(5) * t133 + t136 * t158;
t67 = t133 * t158 + t123;
t56 = t65 * qJD(3);
t51 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t72;
t50 = qJDD(4) * mrSges(6,1) - mrSges(6,3) * t72;
t49 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t71;
t48 = -qJDD(4) * mrSges(6,2) + mrSges(6,3) * t71;
t45 = t62 - t197;
t40 = t57 + t164;
t39 = t136 * t63 + t125;
t38 = t182 * t133;
t29 = -mrSges(5,1) * t71 + mrSges(5,2) * t72;
t11 = (-qJD(5) - t56) * t133 + t136 * t155;
t10 = t133 * t155 + t136 * t56 + t123;
t2 = -t37 * t169 + qJDD(4) * pkin(4) - qJ(5) * t72 + t121 + (-qJD(5) * t128 - t210) * t133;
t1 = [(-m(5) * t107 + mrSges(2,2) * t138 + mrSges(3,2) * t119 + (m(6) * t132 + mrSges(4,2)) * t113 + (pkin(2) * t168 + mrSges(3,1)) * t118 + (pkin(1) * t163 + mrSges(2,1)) * t135 + (m(6) * t116 + t203 - t217) * t112 + t154) * g(1) + (-m(6) * (t156 + t172) - m(5) * (t172 + t173) + mrSges(2,2) * t135 - mrSges(3,1) * t119 + mrSges(3,2) * t118 - m(4) * t172 + (-mrSges(2,1) - t206) * t138 + t142) * g(2) + m(5) * (t14 * t62 + t36 * t57) + (-t66 * t127 - t56 * t128 - t17) * mrSges(4,2) + t11 * t81 + t10 * t83 + t40 * t69 + t57 * t70 + t45 * t28 + t39 * t48 + t38 * t50 + t62 * t29 + (-t84 * t180 + t153 + (-t25 * qJD(4) - t2) * mrSges(6,3) + t205 * mrSges(5,3) + (-m(5) * t30 - t82) * t56 + (m(5) * t205 - t51) * t63) * t133 + (m(5) * (-t180 * t30 + t31 * t56 + t5 * t63) + t56 * t84 + t63 * t49 - t82 * t180 + t204 * qJD(4)) * t136 + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * t131 * mrSges(3,1) - 0.2e1 * t130 * mrSges(3,2) + (t130 ^ 2 + t131 ^ 2) * t206) * pkin(1)) * qJDD(1) + m(4) * (t17 * t66 + t18 * t65 - t46 * t57 + t47 * t56) + m(6) * (t10 * t25 + t11 * t19 + t2 * t38 + t27 * t40 + t3 * t39 + t45 * t8) + (t65 * t127 - t57 * t128) * mrSges(4,1) + t140; (t50 + t51) * t136 + (t48 + t49) * t133 + (m(3) + m(4)) * qJDD(2) + (-t189 * t133 + t188 * t136) * qJD(4) + m(5) * (qJD(4) * t143 + t133 * t5 + t136 * t6) + m(6) * (qJD(4) * t144 + t133 * t3 + t136 * t2) - t163 * g(3); (t112 * mrSges(4,1) + t113 * mrSges(4,2) + t154) * g(1) + t142 * g(2) + (t128 * mrSges(4,1) - t69 - t70) * t47 - t14 * t203 - t116 * t28 + t94 * t48 + t68 * t81 + t67 * t83 + t92 * t50 - pkin(3) * t29 + (-t6 * mrSges(5,3) - t2 * mrSges(6,3) - pkin(7) * t51 + t189 * t46 + t153 + (pkin(4) * t69 - pkin(7) * t84 - t148) * qJD(4)) * t133 + (pkin(7) * t49 - t188 * t46 + t191 * t196 + (-pkin(7) * t82 + t204) * qJD(4)) * t136 + (t128 * t46 - t17) * mrSges(4,2) + t140 + ((t112 * t116 + t113 * t132) * g(1) - t156 * g(2) - t116 * t8 + t19 * t68 + t2 * t92 + t25 * t67 + t3 * t94 - t144 * t46 + (t164 - t47) * t27) * m(6) + ((pkin(3) * t112 - t107) * g(1) - t173 * g(2) + (-t6 * t133 + t194 + (-t133 * t31 - t136 * t30) * qJD(4)) * pkin(7) - t143 * t46 - t36 * t47) * m(5); t6 * mrSges(5,1) + t2 * mrSges(6,1) - t5 * mrSges(5,2) - t3 * mrSges(6,2) - t24 * t83 - t30 * t84 + t31 * t82 + t213 * t72 + t211 * t71 + (Ifges(6,3) + Ifges(5,3)) * qJDD(4) + (t147 + t93) * g(3) + (t50 + (-g(3) * t136 + t2) * m(6)) * pkin(4) + (-m(6) * (-t19 + t24) + t81) * t25 + ((-t99 / 0.2e1 - t100 / 0.2e1 - t60 / 0.2e1 - t61 / 0.2e1 - t27 * mrSges(6,2) - t36 * mrSges(5,2) + (-Ifges(6,5) / 0.2e1 - Ifges(5,5) / 0.2e1) * qJD(4) - t204) * t136 + (t58 / 0.2e1 + t59 / 0.2e1 - t27 * mrSges(6,1) - t36 * mrSges(5,1) + t218 * t176 + (Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1) * qJD(4) + (-m(6) * t27 - t69) * pkin(4) + (-Ifges(6,1) / 0.2e1 - Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t175 + t148) * t133) * t128 + (g(1) * t113 + g(2) * t112) * (t190 * t136 + (m(6) * pkin(4) + t191) * t133); (t133 * t81 - t136 * t83) * t128 + (-t128 * t144 + t195 - t196 + t8) * m(6) + t28;];
tau = t1;
