% Calculate vector of inverse dynamics joint torques for
% S4RRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR2_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR2_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR2_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR2_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR2_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:08
% EndTime: 2019-12-31 17:23:11
% DurationCPUTime: 1.68s
% Computational Cost: add. (2257->270), mult. (3531->383), div. (0->0), fcn. (2018->12), ass. (0->144)
t147 = sin(qJ(3));
t213 = t147 / 0.2e1;
t151 = cos(qJ(3));
t109 = -mrSges(4,1) * t151 + mrSges(4,2) * t147;
t144 = qJ(3) + qJ(4);
t132 = sin(t144);
t134 = cos(t144);
t173 = t134 * mrSges(5,1) - mrSges(5,2) * t132;
t223 = t109 - t173;
t145 = qJ(1) + qJ(2);
t133 = sin(t145);
t135 = cos(t145);
t218 = g(1) * t135 + g(2) * t133;
t222 = mrSges(3,2) - mrSges(5,3) - mrSges(4,3);
t141 = qJD(1) + qJD(2);
t188 = t141 * t147;
t101 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t188;
t187 = t141 * t151;
t102 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t187;
t148 = sin(qJ(2));
t197 = pkin(1) * qJD(1);
t103 = pkin(6) * t141 + t148 * t197;
t152 = cos(qJ(2));
t181 = t152 * t197;
t104 = -pkin(2) * t141 - t181;
t221 = (t141 * mrSges(3,2) + t147 * t101 - t151 * t102) * t152 - m(4) * (t104 * t148 + (t147 ^ 2 + t151 ^ 2) * t103 * t152);
t171 = mrSges(4,1) * t147 + mrSges(4,2) * t151;
t220 = t104 * t171 + qJD(3) * (Ifges(4,5) * t151 - Ifges(4,6) * t147) / 0.2e1;
t219 = -mrSges(3,1) + t223;
t217 = (-mrSges(3,1) + t109) * t141;
t146 = sin(qJ(4));
t150 = cos(qJ(4));
t97 = t146 * t151 + t147 * t150;
t82 = t97 * t141;
t214 = t82 / 0.2e1;
t154 = -pkin(7) - pkin(6);
t96 = -t146 * t147 + t150 * t151;
t81 = t96 * t141;
t211 = mrSges(5,3) * t81;
t210 = Ifges(5,4) * t82;
t209 = pkin(1) * t152;
t149 = sin(qJ(1));
t207 = g(1) * t149;
t205 = t82 * mrSges(5,3);
t126 = pkin(1) * t148 + pkin(6);
t204 = -pkin(7) - t126;
t200 = Ifges(4,4) * t147;
t199 = Ifges(4,4) * t151;
t198 = Ifges(4,2) * t151;
t196 = pkin(1) * qJD(2);
t195 = qJD(3) * pkin(3);
t175 = pkin(7) * t141 + t103;
t70 = t175 * t151;
t194 = t146 * t70;
t182 = qJD(3) * t151;
t139 = qJDD(1) + qJDD(2);
t178 = qJD(1) * t196;
t190 = pkin(1) * qJDD(1);
t94 = t148 * t190 + t152 * t178;
t85 = pkin(6) * t139 + t94;
t48 = -t103 * t182 - t147 * t85;
t193 = t147 * t48;
t192 = t150 * t70;
t183 = qJD(3) * t147;
t47 = -t103 * t183 + t151 * t85;
t191 = t151 * t47;
t184 = t135 * pkin(2) + t133 * pkin(6);
t180 = t152 * t196;
t179 = pkin(3) * t183;
t127 = pkin(3) * t151 + pkin(2);
t177 = qJD(3) * t154;
t174 = qJD(3) * t204;
t172 = t135 * t127 - t133 * t154;
t93 = -t148 * t178 + t152 * t190;
t69 = t175 * t147;
t170 = mrSges(5,1) * t132 + mrSges(5,2) * t134;
t169 = t198 + t200;
t61 = -t69 + t195;
t35 = t150 * t61 - t194;
t36 = t146 * t61 + t192;
t91 = t204 * t147;
t136 = t151 * pkin(7);
t92 = t126 * t151 + t136;
t51 = -t146 * t92 + t150 * t91;
t52 = t146 * t91 + t150 * t92;
t167 = t191 - t193;
t166 = t101 * t151 + t102 * t147;
t110 = t154 * t147;
t111 = pkin(6) * t151 + t136;
t64 = t110 * t150 - t111 * t146;
t65 = t110 * t146 + t111 * t150;
t164 = t147 * (Ifges(4,1) * t151 - t200);
t163 = t97 * qJD(4);
t162 = t96 * qJD(4);
t84 = -pkin(2) * t139 - t93;
t87 = t139 * t151 - t141 * t183;
t88 = t139 * t147 + t141 * t182;
t159 = -t166 * qJD(3) - t147 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t88) + t151 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t87);
t138 = qJDD(3) + qJDD(4);
t140 = qJD(3) + qJD(4);
t29 = t141 * t162 + t146 * t87 + t150 * t88;
t31 = qJDD(3) * pkin(3) - pkin(7) * t88 + t48;
t34 = pkin(7) * t87 + t47;
t3 = qJD(4) * t35 + t146 * t31 + t150 * t34;
t30 = -t141 * t163 - t146 * t88 + t150 * t87;
t4 = -qJD(4) * t36 - t146 * t34 + t150 * t31;
t40 = Ifges(5,2) * t81 + Ifges(5,6) * t140 + t210;
t76 = Ifges(5,4) * t81;
t41 = Ifges(5,1) * t82 + Ifges(5,5) * t140 + t76;
t83 = -t127 * t141 - t181;
t158 = t4 * mrSges(5,1) - t3 * mrSges(5,2) + t35 * t211 + t40 * t214 - t83 * (mrSges(5,1) * t82 + mrSges(5,2) * t81) + Ifges(5,3) * t138 - t82 * (Ifges(5,1) * t81 - t210) / 0.2e1 + Ifges(5,6) * t30 + Ifges(5,5) * t29 - t140 * (Ifges(5,5) * t81 - Ifges(5,6) * t82) / 0.2e1 - (-Ifges(5,2) * t82 + t41 + t76) * t81 / 0.2e1;
t157 = t222 * t133 + t219 * t135;
t156 = (m(4) * pkin(2) + m(5) * t127 - t219) * t133 + (-m(4) * pkin(6) + m(5) * t154 + t222) * t135;
t49 = -pkin(3) * t87 + t84;
t54 = t96 * qJD(3) + t162;
t55 = -t97 * qJD(3) - t163;
t79 = Ifges(4,6) * qJD(3) + t169 * t141;
t114 = Ifges(4,4) * t187;
t80 = Ifges(4,1) * t188 + Ifges(4,5) * qJD(3) + t114;
t155 = -t79 * t183 / 0.2e1 + t87 * t169 / 0.2e1 + t151 * (Ifges(4,4) * t88 + Ifges(4,2) * t87) / 0.2e1 + Ifges(3,3) * t139 + t140 * (Ifges(5,5) * t54 + Ifges(5,6) * t55) / 0.2e1 + t84 * t109 + t93 * mrSges(3,1) - t94 * mrSges(3,2) + t81 * (Ifges(5,4) * t54 + Ifges(5,2) * t55) / 0.2e1 + t83 * (-mrSges(5,1) * t55 + mrSges(5,2) * t54) + t54 * t41 / 0.2e1 + t55 * t40 / 0.2e1 + (Ifges(4,1) * t88 + Ifges(4,4) * t87) * t213 + (Ifges(5,1) * t54 + Ifges(5,4) * t55) * t214 + mrSges(4,3) * t191 + t88 * (Ifges(4,1) * t147 + t199) / 0.2e1 + (t141 * (-Ifges(4,2) * t147 + t199) + t80) * t182 / 0.2e1 + (-t35 * t54 + t36 * t55) * mrSges(5,3) + (0.2e1 * Ifges(4,5) * t213 + Ifges(4,6) * t151) * qJDD(3) + (t164 * t141 / 0.2e1 + t220) * qJD(3) + (t49 * mrSges(5,2) - t4 * mrSges(5,3) + Ifges(5,1) * t29 + Ifges(5,4) * t30 + Ifges(5,5) * t138) * t97 + (-t49 * mrSges(5,1) + t3 * mrSges(5,3) + Ifges(5,4) * t29 + Ifges(5,2) * t30 + Ifges(5,6) * t138) * t96;
t153 = cos(qJ(1));
t137 = t153 * pkin(1);
t108 = -t127 - t209;
t100 = t148 * t196 + t179;
t99 = t151 * t177;
t98 = t147 * t177;
t72 = t96 * t181;
t71 = t97 * t181;
t67 = -t147 * t180 + t151 * t174;
t66 = t147 * t174 + t151 * t180;
t60 = mrSges(5,1) * t140 - t205;
t59 = -mrSges(5,2) * t140 + t211;
t50 = -mrSges(4,1) * t87 + mrSges(4,2) * t88;
t45 = -mrSges(5,1) * t81 + mrSges(5,2) * t82;
t38 = -t150 * t69 - t194;
t37 = t146 * t69 - t192;
t33 = -t65 * qJD(4) - t146 * t98 + t150 * t99;
t32 = t64 * qJD(4) + t146 * t99 + t150 * t98;
t23 = -mrSges(5,2) * t138 + mrSges(5,3) * t30;
t22 = mrSges(5,1) * t138 - mrSges(5,3) * t29;
t9 = -t52 * qJD(4) - t146 * t66 + t150 * t67;
t8 = t51 * qJD(4) + t146 * t67 + t150 * t66;
t7 = -mrSges(5,1) * t30 + mrSges(5,2) * t29;
t1 = [m(5) * (t100 * t83 + t108 * t49 + t3 * t52 + t35 * t9 + t36 * t8 + t4 * t51) + t155 + ((mrSges(3,1) * t152 - mrSges(3,2) * t148) * t139 + (m(4) + m(5)) * t207 + (-g(2) * t153 + t148 * t94 + t152 * t93 + t207) * m(3) + (t148 * t217 - t221) * qJD(2)) * pkin(1) + (m(4) * t167 + t159) * t126 - mrSges(4,3) * t193 + t108 * t7 + t100 * t45 + t8 * t59 + t9 * t60 + t51 * t22 + t52 * t23 + (mrSges(2,1) * t149 + mrSges(2,2) * t153 + t156) * g(1) + (-m(5) * (t137 + t172) - mrSges(2,1) * t153 + mrSges(2,2) * t149 - m(4) * (t137 + t184) + t157) * g(2) + Ifges(2,3) * qJDD(1) + (m(4) * t84 + t50) * (-pkin(2) - t209); -m(5) * (-t35 * t71 + t36 * t72) + t155 - t127 * t7 + t64 * t22 + t65 * t23 - pkin(2) * t50 + t159 * pkin(6) + t156 * g(1) + (-m(4) * t184 - m(5) * t172 + t157) * g(2) + (-t48 * mrSges(4,3) + t45 * t195) * t147 + m(5) * (-t127 * t49 + t83 * t179 + t3 * t65 + t32 * t36 + t33 * t35 + t4 * t64) + (-t72 + t32) * t59 + (t71 + t33) * t60 + ((-m(5) * t83 - t217 - t45) * t148 + t221) * t197 + m(4) * (-pkin(2) * t84 + t167 * pkin(6)); (t79 * t213 + (t198 * t213 - t164 / 0.2e1) * t141 - (t114 + t80) * t151 / 0.2e1 - t220) * t141 + t223 * g(3) + t158 + t166 * t103 + (-t45 * t188 + t146 * t23 + t150 * t22 + (-g(3) * t151 + t146 * t3 + t150 * t4 + (-t141 * t83 + t218) * t147) * m(5) + (-t146 * t60 + t150 * t59 + (-t146 * t35 + t150 * t36) * m(5)) * qJD(4)) * pkin(3) - m(5) * (t35 * t37 + t36 * t38) + t36 * t205 + Ifges(4,6) * t87 + Ifges(4,5) * t88 - t38 * t59 - t37 * t60 - t47 * mrSges(4,2) + t48 * mrSges(4,1) + Ifges(4,3) * qJDD(3) + t218 * (t170 + t171); t158 + (t60 + t205) * t36 - g(3) * t173 - t35 * t59 + t218 * t170;];
tau = t1;
