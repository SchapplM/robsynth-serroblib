% Calculate vector of inverse dynamics joint torques for
% S5RPRPR2
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
% Datum: 2022-01-23 09:19
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:18:49
% EndTime: 2022-01-23 09:18:52
% DurationCPUTime: 1.66s
% Computational Cost: add. (2279->248), mult. (3954->329), div. (0->0), fcn. (2390->16), ass. (0->124)
t145 = sin(qJ(3));
t148 = cos(qJ(3));
t140 = sin(pkin(8));
t185 = pkin(1) * t140;
t169 = qJD(1) * t185;
t142 = cos(pkin(8));
t119 = pkin(1) * t142 + pkin(2);
t99 = t119 * qJD(1);
t66 = -t145 * t169 + t148 * t99;
t160 = qJD(4) - t66;
t144 = sin(qJ(5));
t147 = cos(qJ(5));
t139 = sin(pkin(9));
t143 = -pkin(7) - qJ(4);
t95 = t143 * t139;
t141 = cos(pkin(9));
t130 = t141 * pkin(7);
t96 = qJ(4) * t141 + t130;
t57 = -t144 * t96 + t147 * t95;
t86 = -t139 * t144 + t141 * t147;
t197 = qJD(5) * t57 + t160 * t86;
t58 = t144 * t95 + t147 * t96;
t87 = t139 * t147 + t141 * t144;
t196 = -qJD(5) * t58 - t160 * t87;
t159 = -t141 * mrSges(5,1) + t139 * mrSges(5,2);
t195 = mrSges(4,1) - t159;
t168 = qJD(3) * t185;
t194 = -qJD(1) * t168 + t119 * qJDD(1);
t136 = pkin(9) + qJ(5);
t125 = sin(t136);
t127 = cos(t136);
t193 = t127 * mrSges(6,1) - t125 * mrSges(6,2);
t192 = mrSges(4,2) - mrSges(6,3) - mrSges(5,3);
t191 = -t193 - t195;
t190 = m(3) * pkin(1);
t137 = qJD(1) + qJD(3);
t172 = t139 ^ 2 + t141 ^ 2;
t189 = t172 * t137;
t80 = t145 * t119 + t148 * t185;
t71 = t87 * t137;
t187 = t71 / 0.2e1;
t186 = Ifges(6,4) * t71;
t184 = pkin(4) * t141;
t122 = t141 * qJDD(2);
t133 = qJDD(1) + qJDD(3);
t166 = qJDD(1) * t185;
t171 = qJD(3) * t148;
t39 = t145 * t194 + t148 * t166 + t99 * t171;
t30 = qJ(4) * t133 + qJD(4) * t137 + t39;
t19 = -t139 * t30 + t122;
t180 = t139 * t19;
t20 = qJDD(2) * t139 + t141 * t30;
t178 = t141 * t20;
t177 = t145 * t99;
t67 = t148 * t169 + t177;
t56 = qJ(4) * t137 + t67;
t51 = qJD(2) * t139 + t141 * t56;
t176 = t133 * t139;
t175 = t133 * t141;
t138 = qJ(1) + pkin(8);
t129 = qJ(3) + t138;
t116 = sin(t129);
t117 = cos(t129);
t174 = pkin(3) * t117 + qJ(4) * t116;
t128 = cos(t138);
t149 = cos(qJ(1));
t173 = pkin(1) * t149 + pkin(2) * t128;
t170 = m(4) + m(5) + m(6);
t167 = m(3) + t170;
t118 = pkin(3) + t184;
t82 = t86 * qJD(5);
t47 = t133 * t87 + t137 * t82;
t83 = t87 * qJD(5);
t48 = t133 * t86 - t137 * t83;
t12 = -t48 * mrSges(6,1) + mrSges(6,2) * t47;
t165 = -t116 * t143 + t117 * t118;
t164 = t172 * t133;
t79 = t119 * t148 - t145 * t185;
t70 = t86 * t137;
t162 = mrSges(6,1) * t70 - mrSges(6,2) * t71 + t137 * t195;
t77 = -pkin(3) - t79;
t78 = -mrSges(5,1) * t175 + mrSges(5,2) * t176;
t161 = -g(1) * t116 + g(2) * t117;
t124 = t141 * qJD(2);
t50 = -t139 * t56 + t124;
t157 = -t139 * t50 + t141 * t51;
t43 = t124 + (-pkin(7) * t137 - t56) * t139;
t44 = t130 * t137 + t51;
t9 = -t144 * t44 + t147 * t43;
t10 = t144 * t43 + t147 * t44;
t76 = qJ(4) + t80;
t59 = (-pkin(7) - t76) * t139;
t60 = t141 * t76 + t130;
t25 = -t144 * t60 + t147 * t59;
t26 = t144 * t59 + t147 * t60;
t40 = -qJD(3) * t177 - t145 * t166 + t148 * t194;
t72 = t119 * t171 - t145 * t168;
t156 = qJDD(4) - t40;
t153 = t116 * t192 + t117 * t191;
t152 = (m(5) * pkin(3) + m(6) * t118 - t191) * t116 + (-m(5) * qJ(4) + m(6) * t143 + t192) * t117;
t13 = t122 + (-pkin(7) * t133 - t30) * t139;
t14 = pkin(7) * t175 + t20;
t2 = qJD(5) * t9 + t13 * t144 + t14 * t147;
t24 = -t118 * t133 + t156;
t3 = -qJD(5) * t10 + t13 * t147 - t14 * t144;
t33 = -pkin(3) * t133 + t156;
t37 = Ifges(6,2) * t70 + Ifges(6,6) * qJD(5) + t186;
t68 = Ifges(6,4) * t70;
t38 = Ifges(6,1) * t71 + Ifges(6,5) * qJD(5) + t68;
t52 = -t118 * t137 + t160;
t151 = Ifges(4,3) * t133 + t52 * (mrSges(6,1) * t83 + mrSges(6,2) * t82) + mrSges(5,3) * t178 + t82 * t38 / 0.2e1 - t83 * t37 / 0.2e1 + t33 * t159 + t70 * (Ifges(6,4) * t82 - Ifges(6,2) * t83) / 0.2e1 + (Ifges(6,1) * t82 - Ifges(6,4) * t83) * t187 + t40 * mrSges(4,1) + qJD(5) * (Ifges(6,5) * t82 - Ifges(6,6) * t83) / 0.2e1 + (Ifges(5,1) * t139 + Ifges(5,4) * t141) * t176 + (Ifges(5,4) * t139 + Ifges(5,2) * t141) * t175 + (-t10 * t83 - t82 * t9) * mrSges(6,3) + (mrSges(6,2) * t24 - mrSges(6,3) * t3 + Ifges(6,1) * t47 + Ifges(6,4) * t48 + Ifges(6,5) * qJDD(5)) * t87 + (-mrSges(6,1) * t24 + mrSges(6,3) * t2 + Ifges(6,4) * t47 + Ifges(6,2) * t48 + Ifges(6,6) * qJDD(5)) * t86;
t146 = sin(qJ(1));
t126 = sin(t138);
t73 = t80 * qJD(3);
t69 = qJD(4) + t72;
t65 = t77 - t184;
t62 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t71;
t61 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t70;
t55 = -pkin(3) * t137 + t160;
t35 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t48;
t34 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t47;
t5 = -qJD(5) * t26 - t69 * t87;
t4 = qJD(5) * t25 + t69 * t86;
t1 = [m(4) * (t39 * t80 + t40 * t79 - t66 * t73 + t67 * t72) + m(6) * (t10 * t4 + t2 * t26 + t24 * t65 + t25 * t3 + t5 * t9 + t52 * t73) - t162 * t73 + t79 * t133 * mrSges(4,1) + t77 * t78 + t5 * t62 + t65 * t12 + t4 * t61 + t25 * t34 + t26 * t35 + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * t142 * mrSges(3,1) - 0.2e1 * t140 * mrSges(3,2) + (t140 ^ 2 + t142 ^ 2) * t190) * pkin(1)) * qJDD(1) + t151 + m(5) * (t33 * t77 + t55 * t73 + (t20 * t76 + t51 * t69) * t141 + (-t19 * t76 - t50 * t69) * t139) + (-t133 * t80 - t137 * t72 - t39) * mrSges(4,2) + (t164 * t76 + t189 * t69 - t180) * mrSges(5,3) + (mrSges(2,2) * t149 + mrSges(3,2) * t128 + (pkin(2) * t170 + mrSges(3,1)) * t126 + (pkin(1) * t167 + mrSges(2,1)) * t146 + t152) * g(1) + (-m(6) * (t165 + t173) + mrSges(2,2) * t146 - mrSges(3,1) * t128 + mrSges(3,2) * t126 - m(5) * (t173 + t174) - m(4) * t173 + (-mrSges(2,1) - t190) * t149 + t153) * g(2); t86 * t34 + t87 * t35 + t82 * t61 - t83 * t62 + (m(3) + m(4)) * qJDD(2) + m(5) * (t139 * t20 + t141 * t19) + m(6) * (t10 * t82 + t2 * t87 + t3 * t86 - t83 * t9) - t167 * g(3); t196 * t62 + t197 * t61 + (qJ(4) * t164 + t160 * t189 - t180) * mrSges(5,3) + t152 * g(1) + t153 * g(2) + (t137 * t66 - t39) * mrSges(4,2) - t118 * t12 - pkin(3) * t78 + t162 * t67 + t57 * t34 + t58 * t35 + t151 + (-t165 * g(2) + t10 * t197 - t118 * t24 + t196 * t9 + t2 * t58 + t3 * t57 - t52 * t67) * m(6) + (-t174 * g(2) - pkin(3) * t33 + (t178 - t180) * qJ(4) - t55 * t67 + t160 * t157) * m(5); -t172 * t137 ^ 2 * mrSges(5,3) - t70 * t61 + t71 * t62 + t12 + t78 + (-t10 * t70 + t71 * t9 + t161 + t24) * m(6) + (-t137 * t157 + t161 + t33) * m(5); Ifges(6,5) * t47 + Ifges(6,6) * t48 + Ifges(6,3) * qJDD(5) - t2 * mrSges(6,2) + t3 * mrSges(6,1) - t52 * (mrSges(6,1) * t71 + mrSges(6,2) * t70) - t71 * (Ifges(6,1) * t70 - t186) / 0.2e1 + t37 * t187 - qJD(5) * (Ifges(6,5) * t70 - Ifges(6,6) * t71) / 0.2e1 - t9 * t61 + t10 * t62 - g(3) * t193 + (t10 * t71 + t70 * t9) * mrSges(6,3) - (-Ifges(6,2) * t71 + t38 + t68) * t70 / 0.2e1 + (g(1) * t117 + g(2) * t116) * (mrSges(6,1) * t125 + mrSges(6,2) * t127);];
tau = t1;
