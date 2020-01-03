% Calculate vector of inverse dynamics joint torques for
% S5RPRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:35:57
% EndTime: 2020-01-03 11:36:06
% DurationCPUTime: 2.18s
% Computational Cost: add. (2271->279), mult. (3897->391), div. (0->0), fcn. (2232->14), ass. (0->140)
t139 = sin(qJ(5));
t142 = cos(qJ(5));
t137 = cos(pkin(9));
t185 = qJD(4) * t137;
t191 = t137 * t142;
t138 = cos(pkin(8));
t120 = pkin(1) * t138 + pkin(2);
t100 = t120 * qJD(1);
t140 = sin(qJ(3));
t143 = cos(qJ(3));
t136 = sin(pkin(8));
t211 = pkin(1) * t136;
t178 = qJD(1) * t211;
t65 = t100 * t143 - t140 * t178;
t202 = t100 * t140;
t66 = t143 * t178 + t202;
t192 = t137 * t139;
t135 = sin(pkin(9));
t92 = -pkin(4) * t137 - pkin(7) * t135 - pkin(3);
t67 = -qJ(4) * t192 + t142 * t92;
t222 = t67 * qJD(5) - t139 * t66 + t142 * t185 - t65 * t191;
t68 = qJ(4) * t191 + t139 * t92;
t221 = -t68 * qJD(5) - t139 * t185 - t142 * t66 + t65 * t192;
t177 = qJD(3) * t211;
t220 = -qJD(1) * t177 + t120 * qJDD(1);
t164 = qJD(4) - t65;
t133 = qJD(1) + qJD(3);
t102 = -t133 * t137 + qJD(5);
t48 = qJ(4) * t133 + t66;
t35 = -t137 * qJD(2) + t135 * t48;
t219 = t35 * (mrSges(6,1) * t142 - mrSges(6,2) * t139) + t102 * (-Ifges(6,5) * t139 - Ifges(6,6) * t142) / 0.2e1;
t209 = Ifges(6,4) * t142;
t210 = Ifges(6,4) * t139;
t218 = t142 * (-Ifges(6,1) * t139 - t209) / 0.2e1 - t139 * (-Ifges(6,2) * t142 - t210) / 0.2e1;
t217 = m(3) * pkin(1);
t134 = qJ(1) + pkin(8);
t126 = qJ(3) + t134;
t118 = sin(t126);
t200 = t118 * t137;
t201 = t118 * t135;
t216 = pkin(4) * t200 + pkin(7) * t201;
t119 = cos(t126);
t198 = t119 * t137;
t199 = t119 * t135;
t215 = pkin(4) * t198 + pkin(7) * t199;
t30 = t92 * t133 + t164;
t36 = qJD(2) * t135 + t137 * t48;
t8 = -t139 * t36 + t142 * t30;
t214 = t8 * qJD(5);
t131 = t135 ^ 2;
t187 = t137 ^ 2 + t131;
t213 = t187 * t133;
t80 = t140 * t120 + t143 * t211;
t212 = m(3) + m(4);
t130 = qJDD(1) + qJDD(3);
t175 = qJDD(1) * t211;
t186 = qJD(3) * t143;
t28 = t100 * t186 + t220 * t140 + t143 * t175;
t19 = qJ(4) * t130 + qJD(4) * t133 + t28;
t13 = -t137 * qJDD(2) + t135 * t19;
t208 = t13 * t135;
t207 = t13 * t137;
t206 = t135 * t35;
t14 = qJDD(2) * t135 + t137 * t19;
t205 = t137 * t14;
t152 = (t142 * Ifges(6,1) - t210) * t135;
t204 = t139 * (Ifges(6,5) * t102 + t133 * t152);
t151 = (-t139 * Ifges(6,2) + t209) * t135;
t203 = t142 * (Ifges(6,6) * t102 + t133 * t151);
t197 = t130 * t135;
t196 = t130 * t137;
t195 = t133 * t135;
t194 = t135 * t139;
t193 = t135 * t142;
t190 = t119 * pkin(3) + t118 * qJ(4);
t124 = sin(t134);
t141 = sin(qJ(1));
t189 = t141 * pkin(1) + pkin(2) * t124;
t125 = cos(t134);
t144 = cos(qJ(1));
t188 = t144 * pkin(1) + pkin(2) * t125;
t184 = qJD(5) * t133;
t183 = qJD(5) * t135;
t182 = qJD(5) * t142;
t69 = (t130 * t142 - t139 * t184) * t135;
t70 = (-t130 * t139 - t133 * t182) * t135;
t99 = qJDD(5) - t196;
t181 = Ifges(6,5) * t69 + Ifges(6,6) * t70 + Ifges(6,3) * t99;
t180 = mrSges(6,3) * t194;
t179 = mrSges(6,3) * t193;
t176 = -mrSges(2,1) - t217;
t172 = t187 * t130;
t171 = t118 * pkin(3) - qJ(4) * t119;
t79 = t120 * t143 - t140 * t211;
t9 = t139 * t30 + t142 * t36;
t170 = t9 * mrSges(6,3) * t182;
t168 = t188 + t190;
t78 = -mrSges(5,1) * t196 + mrSges(5,2) * t197;
t167 = g(2) * t119 + g(3) * t118;
t166 = -t139 * t8 + t142 * t9;
t165 = g(3) * t119 + t208;
t163 = -mrSges(5,1) * t137 + mrSges(5,2) * t135;
t162 = mrSges(6,1) * t139 + mrSges(6,2) * t142;
t161 = t137 * t36 + t206;
t58 = -mrSges(6,2) * t102 - t133 * t180;
t59 = mrSges(6,1) * t102 - t133 * t179;
t160 = t139 * t59 - t142 * t58;
t29 = -qJD(3) * t202 - t140 * t175 + t220 * t143;
t72 = t120 * t186 - t140 * t177;
t60 = -t118 * t192 - t119 * t142;
t61 = t118 * t191 - t119 * t139;
t159 = -t118 * mrSges(4,1) - mrSges(5,1) * t200 - t61 * mrSges(6,1) - t119 * mrSges(4,2) - t60 * mrSges(6,2) - mrSges(6,3) * t201;
t158 = t171 + t189;
t49 = -t79 + t92;
t76 = qJ(4) + t80;
t21 = t139 * t49 + t76 * t191;
t20 = t142 * t49 - t76 * t192;
t71 = qJD(4) + t72;
t157 = (t13 * t76 + t35 * t71) * t135;
t150 = qJDD(4) - t29;
t62 = -t118 * t142 + t119 * t192;
t63 = t118 * t139 + t119 * t191;
t149 = -t119 * mrSges(4,1) - mrSges(5,1) * t198 - t63 * mrSges(6,1) + t62 * mrSges(6,2) - mrSges(6,3) * t199 + (mrSges(4,2) - mrSges(5,3)) * t118;
t25 = -pkin(3) * t130 + t150;
t12 = t92 * t130 + t150;
t3 = t12 * t139 + t14 * t142 + t214;
t4 = -t9 * qJD(5) + t12 * t142 - t139 * t14;
t82 = t162 * t135;
t146 = t3 * (mrSges(6,2) * t137 - t180) + mrSges(5,3) * t205 + Ifges(4,3) * t130 - (Ifges(6,4) * t69 + Ifges(6,2) * t70 + Ifges(6,6) * t99) * t194 / 0.2e1 + (Ifges(6,1) * t69 + Ifges(6,4) * t70 + Ifges(6,5) * t99) * t193 / 0.2e1 - t137 * t181 / 0.2e1 + t25 * t163 + t4 * (-mrSges(6,1) * t137 - t179) + t29 * mrSges(4,1) + t69 * (-Ifges(6,5) * t137 + t152) / 0.2e1 + t70 * (-Ifges(6,6) * t137 + t151) / 0.2e1 + t99 * (-Ifges(6,3) * t137 + (Ifges(6,5) * t142 - Ifges(6,6) * t139) * t135) / 0.2e1 + (Ifges(5,1) * t135 + Ifges(5,4) * t137) * t197 + (Ifges(5,4) * t135 + Ifges(5,2) * t137) * t196 + t180 * t214 + t13 * t82 + t219 * t183 + t218 * t131 * t184 - (t204 + t203) * t183 / 0.2e1;
t81 = t163 * t133;
t77 = -pkin(3) - t79;
t73 = t80 * qJD(3);
t64 = t162 * t195;
t47 = -pkin(3) * t133 + t164;
t38 = -mrSges(6,2) * t99 + mrSges(6,3) * t70;
t37 = mrSges(6,1) * t99 - mrSges(6,3) * t69;
t27 = -mrSges(6,1) * t70 + mrSges(6,2) * t69;
t6 = -t21 * qJD(5) + t142 * t73 - t71 * t192;
t5 = t20 * qJD(5) + t139 * t73 + t71 * t191;
t1 = [(-m(5) * t158 + mrSges(5,2) * t201 - mrSges(2,2) * t144 - mrSges(3,1) * t124 - mrSges(3,2) * t125 - m(6) * (t158 + t216) - m(4) * t189 + t176 * t141 + t159) * g(3) + (t76 * t27 + t71 * t64 - t170) * t135 + (t79 * t130 - t73 * t133) * mrSges(4,1) + m(4) * (t28 * t80 + t29 * t79 - t65 * t73 + t66 * t72) + t77 * t78 + t73 * t81 + t5 * t58 + t6 * t59 + t20 * t37 + t21 * t38 + t146 + (t76 * t172 + t71 * t213 + t165) * mrSges(5,3) + (Ifges(2,3) + Ifges(3,3) + (0.2e1 * t138 * mrSges(3,1) - 0.2e1 * t136 * mrSges(3,2) + (t136 ^ 2 + t138 ^ 2) * t217) * pkin(1)) * qJDD(1) + (mrSges(2,2) * t141 - mrSges(3,1) * t125 + mrSges(3,2) * t124 - m(5) * t168 + mrSges(5,2) * t199 - m(6) * (t168 + t215) - m(4) * t188 + t176 * t144 + t149) * g(2) + m(6) * (t20 * t4 + t21 * t3 + t5 * t9 + t6 * t8 + t157) + m(5) * (t25 * t77 + t47 * t73 + (t14 * t76 + t36 * t71) * t137 + t157) + (-t80 * t130 - t72 * t133 - t28) * mrSges(4,2); -t137 * t27 + t212 * qJDD(2) + (-t139 * t37 + t142 * t38 + (-t139 * t58 - t142 * t59) * qJD(5)) * t135 + m(5) * (t135 * t14 - t207) + m(6) * (-t207 + (-t139 * t4 + t142 * t3 + (-t139 * t9 - t142 * t8) * qJD(5)) * t135) + (-m(5) - m(6) - t212) * g(1); t221 * t59 + t222 * t58 + (qJ(4) * t172 + t164 * t213 + t165) * mrSges(5,3) + (t167 * mrSges(5,2) + qJ(4) * t27 + t164 * t64 - t170) * t135 + (t133 * t65 - t28) * mrSges(4,2) + t159 * g(3) - pkin(3) * t78 + t67 * t37 + t68 * t38 + t146 + t149 * g(2) + (t133 * mrSges(4,1) - t81) * t66 + (-t65 * t206 + (-t171 - t216) * g(3) + (-t190 - t215) * g(2) + t3 * t68 + t4 * t67 + (qJ(4) * t13 + qJD(4) * t35) * t135 + t222 * t9 + t221 * t8) * m(6) + (-t171 * g(3) - t190 * g(2) - t25 * pkin(3) + (t205 + t208) * qJ(4) - t47 * t66 + t164 * t161) * m(5); t139 * t38 + t142 * t37 - t160 * qJD(5) + (t166 * qJD(5) + t139 * t3 + t142 * t4 + t167) * m(6) + (t167 + t25) * m(5) + t78 + (-t135 * t64 + t160 * t137 - m(5) * t161 - m(6) * (t9 * t191 - t8 * t192 + t206) - mrSges(5,3) * t213) * t133; -t3 * mrSges(6,2) + t4 * mrSges(6,1) - t8 * t58 + t9 * t59 + g(1) * t82 - g(2) * (mrSges(6,1) * t60 - mrSges(6,2) * t61) - g(3) * (mrSges(6,1) * t62 + mrSges(6,2) * t63) + (t204 / 0.2e1 + t203 / 0.2e1 - t218 * t195 + t166 * mrSges(6,3) - t219) * t195 + t181;];
tau = t1;
