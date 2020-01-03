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
% Datum: 2020-01-03 11:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:33:38
% EndTime: 2020-01-03 11:33:46
% DurationCPUTime: 1.98s
% Computational Cost: add. (2279->254), mult. (3954->334), div. (0->0), fcn. (2390->16), ass. (0->129)
t146 = cos(pkin(8));
t122 = pkin(1) * t146 + pkin(2);
t100 = t122 * qJD(1);
t149 = sin(qJ(3));
t152 = cos(qJ(3));
t144 = sin(pkin(8));
t194 = pkin(1) * t144;
t176 = qJD(1) * t194;
t66 = t100 * t152 - t149 * t176;
t165 = qJD(4) - t66;
t143 = sin(pkin(9));
t145 = cos(pkin(9));
t164 = -mrSges(5,1) * t145 + t143 * mrSges(5,2);
t207 = mrSges(4,1) - t164;
t140 = pkin(9) + qJ(5);
t128 = sin(t140);
t130 = cos(t140);
t206 = mrSges(6,1) * t130 - t128 * mrSges(6,2);
t205 = -mrSges(4,2) + mrSges(6,3);
t148 = sin(qJ(5));
t151 = cos(qJ(5));
t147 = -pkin(7) - qJ(4);
t96 = t147 * t143;
t133 = t145 * pkin(7);
t97 = qJ(4) * t145 + t133;
t57 = -t148 * t97 + t151 * t96;
t86 = -t143 * t148 + t145 * t151;
t204 = t57 * qJD(5) + t165 * t86;
t58 = t148 * t96 + t151 * t97;
t87 = t143 * t151 + t145 * t148;
t203 = -t58 * qJD(5) - t165 * t87;
t175 = qJD(3) * t194;
t202 = -qJD(1) * t175 + t122 * qJDD(1);
t201 = -t206 - t207;
t200 = m(3) * pkin(1);
t141 = qJD(1) + qJD(3);
t178 = t143 ^ 2 + t145 ^ 2;
t199 = t178 * t141;
t80 = t149 * t122 + t152 * t194;
t71 = t87 * t141;
t197 = t71 / 0.2e1;
t196 = m(3) + m(4);
t195 = Ifges(6,4) * t71;
t193 = pkin(4) * t145;
t142 = qJ(1) + pkin(8);
t132 = qJ(3) + t142;
t120 = cos(t132);
t192 = g(3) * t120;
t125 = t145 * qJDD(2);
t137 = qJDD(1) + qJDD(3);
t173 = qJDD(1) * t194;
t177 = qJD(3) * t152;
t39 = t100 * t177 + t149 * t202 + t152 * t173;
t30 = qJ(4) * t137 + qJD(4) * t141 + t39;
t19 = -t143 * t30 + t125;
t187 = t143 * t19;
t20 = qJDD(2) * t143 + t145 * t30;
t186 = t145 * t20;
t184 = t100 * t149;
t67 = t152 * t176 + t184;
t56 = qJ(4) * t141 + t67;
t51 = qJD(2) * t143 + t145 * t56;
t119 = sin(t132);
t121 = pkin(3) + t193;
t185 = t119 * t121 + t120 * t147;
t183 = t137 * t143;
t182 = t137 * t145;
t181 = pkin(3) * t120 + qJ(4) * t119;
t129 = sin(t142);
t150 = sin(qJ(1));
t180 = pkin(1) * t150 + pkin(2) * t129;
t131 = cos(t142);
t153 = cos(qJ(1));
t179 = pkin(1) * t153 + pkin(2) * t131;
t174 = -mrSges(2,1) - t200;
t82 = t86 * qJD(5);
t47 = t137 * t87 + t141 * t82;
t83 = t87 * qJD(5);
t48 = t137 * t86 - t141 * t83;
t12 = -t48 * mrSges(6,1) + mrSges(6,2) * t47;
t172 = -t119 * t147 + t120 * t121;
t171 = t178 * t137;
t170 = pkin(3) * t119 - qJ(4) * t120;
t79 = t122 * t152 - t149 * t194;
t70 = t86 * t141;
t168 = mrSges(6,1) * t70 - mrSges(6,2) * t71 + t141 * t207;
t77 = -pkin(3) - t79;
t78 = -mrSges(5,1) * t182 + mrSges(5,2) * t183;
t167 = g(2) * t120 + g(3) * t119;
t166 = -t187 + t192;
t127 = t145 * qJD(2);
t50 = -t143 * t56 + t127;
t161 = -t143 * t50 + t145 * t51;
t43 = t127 + (-pkin(7) * t141 - t56) * t143;
t44 = t133 * t141 + t51;
t9 = -t148 * t44 + t151 * t43;
t10 = t148 * t43 + t151 * t44;
t76 = qJ(4) + t80;
t59 = (-pkin(7) - t76) * t143;
t60 = t145 * t76 + t133;
t25 = -t148 * t60 + t151 * t59;
t26 = t148 * t59 + t151 * t60;
t40 = -qJD(3) * t184 - t149 * t173 + t152 * t202;
t72 = t122 * t177 - t149 * t175;
t160 = qJDD(4) - t40;
t157 = t119 * t201 + t120 * t205;
t156 = (-mrSges(5,3) - t205) * t119 + t201 * t120;
t13 = t125 + (-pkin(7) * t137 - t30) * t143;
t14 = pkin(7) * t182 + t20;
t2 = qJD(5) * t9 + t13 * t148 + t14 * t151;
t24 = -t121 * t137 + t160;
t3 = -qJD(5) * t10 + t13 * t151 - t14 * t148;
t33 = -pkin(3) * t137 + t160;
t37 = Ifges(6,2) * t70 + Ifges(6,6) * qJD(5) + t195;
t68 = Ifges(6,4) * t70;
t38 = Ifges(6,1) * t71 + Ifges(6,5) * qJD(5) + t68;
t52 = -t121 * t141 + t165;
t155 = Ifges(4,3) * t137 + t52 * (mrSges(6,1) * t83 + mrSges(6,2) * t82) + mrSges(5,3) * t186 + t82 * t38 / 0.2e1 - t83 * t37 / 0.2e1 + t33 * t164 + t70 * (Ifges(6,4) * t82 - Ifges(6,2) * t83) / 0.2e1 + (Ifges(6,1) * t82 - Ifges(6,4) * t83) * t197 + t40 * mrSges(4,1) + qJD(5) * (Ifges(6,5) * t82 - Ifges(6,6) * t83) / 0.2e1 + (Ifges(5,1) * t143 + Ifges(5,4) * t145) * t183 + (Ifges(5,4) * t143 + Ifges(5,2) * t145) * t182 + (-t10 * t83 - t82 * t9) * mrSges(6,3) + (mrSges(6,2) * t24 - mrSges(6,3) * t3 + Ifges(6,1) * t47 + Ifges(6,4) * t48 + Ifges(6,5) * qJDD(5)) * t87 + (-mrSges(6,1) * t24 + mrSges(6,3) * t2 + Ifges(6,4) * t47 + Ifges(6,2) * t48 + Ifges(6,6) * qJDD(5)) * t86;
t73 = t80 * qJD(3);
t69 = qJD(4) + t72;
t65 = t77 - t193;
t62 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t71;
t61 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t70;
t55 = -pkin(3) * t141 + t165;
t35 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t48;
t34 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t47;
t5 = -qJD(5) * t26 - t69 * t87;
t4 = qJD(5) * t25 + t69 * t86;
t1 = [t155 + (mrSges(2,2) * t150 - mrSges(3,1) * t131 + mrSges(3,2) * t129 - m(6) * (t172 + t179) - m(4) * t179 - m(5) * (t179 + t181) + t174 * t153 + t156) * g(2) + (-m(5) * (t170 + t180) - mrSges(2,2) * t153 - mrSges(3,1) * t129 - mrSges(3,2) * t131 - m(6) * (t180 + t185) - m(4) * t180 + t174 * t150 + t157) * g(3) - t168 * t73 + m(6) * (t10 * t4 + t2 * t26 + t24 * t65 + t25 * t3 + t5 * t9 + t52 * t73) + m(4) * (t39 * t80 + t40 * t79 - t66 * t73 + t67 * t72) + t77 * t78 + t4 * t61 + t5 * t62 + t65 * t12 + t25 * t34 + t26 * t35 + t79 * t137 * mrSges(4,1) + (t171 * t76 + t199 * t69 + t166) * mrSges(5,3) + (Ifges(2,3) + Ifges(3,3) + (0.2e1 * t146 * mrSges(3,1) - 0.2e1 * t144 * mrSges(3,2) + (t144 ^ 2 + t146 ^ 2) * t200) * pkin(1)) * qJDD(1) + m(5) * (t33 * t77 + t55 * t73 + (t20 * t76 + t51 * t69) * t145 + (-t19 * t76 - t50 * t69) * t143) + (-t137 * t80 - t141 * t72 - t39) * mrSges(4,2); t86 * t34 + t87 * t35 + t82 * t61 - t83 * t62 + t196 * qJDD(2) + m(5) * (t143 * t20 + t145 * t19) + m(6) * (t10 * t82 + t2 * t87 + t3 * t86 - t83 * t9) + (-m(5) - m(6) - t196) * g(1); t155 + t203 * t62 + t204 * t61 - t121 * t12 - pkin(3) * t78 + (qJ(4) * t171 + t165 * t199 + t166) * mrSges(5,3) + t57 * t34 + t58 * t35 + t168 * t67 + (t141 * t66 - t39) * mrSges(4,2) + t157 * g(3) + t156 * g(2) + (-t172 * g(2) - t185 * g(3) + t10 * t204 - t121 * t24 + t2 * t58 + t203 * t9 + t3 * t57 - t52 * t67) * m(6) + (-pkin(3) * t33 + (t186 - t187) * qJ(4) - t55 * t67 - t170 * g(3) - t181 * g(2) + t165 * t161) * m(5); -t178 * t141 ^ 2 * mrSges(5,3) - t70 * t61 + t71 * t62 + t12 + t78 + (-t10 * t70 + t9 * t71 + t167 + t24) * m(6) + (-t141 * t161 + t167 + t33) * m(5); Ifges(6,5) * t47 + Ifges(6,6) * t48 + Ifges(6,3) * qJDD(5) - t2 * mrSges(6,2) + t3 * mrSges(6,1) - t52 * (mrSges(6,1) * t71 + mrSges(6,2) * t70) - t71 * (Ifges(6,1) * t70 - t195) / 0.2e1 + t37 * t197 - qJD(5) * (Ifges(6,5) * t70 - Ifges(6,6) * t71) / 0.2e1 - t9 * t61 + t10 * t62 - g(1) * t206 + (t10 * t71 + t70 * t9) * mrSges(6,3) - (-Ifges(6,2) * t71 + t38 + t68) * t70 / 0.2e1 + (g(2) * t119 - t192) * (mrSges(6,1) * t128 + mrSges(6,2) * t130);];
tau = t1;
