% Calculate vector of inverse dynamics joint torques for
% S5RPPRP2
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP2_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP2_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP2_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP2_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP2_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:05
% EndTime: 2019-12-31 17:49:12
% DurationCPUTime: 4.01s
% Computational Cost: add. (1537->278), mult. (3306->342), div. (0->0), fcn. (2122->12), ass. (0->117)
t172 = mrSges(5,1) + mrSges(6,1);
t169 = Ifges(5,1) + Ifges(6,1);
t167 = Ifges(6,4) + Ifges(5,5);
t88 = sin(pkin(7));
t69 = pkin(1) * t88 + qJ(3);
t59 = qJD(1) * qJD(3) + qJDD(1) * t69;
t168 = -Ifges(5,4) + Ifges(6,5);
t87 = sin(pkin(8));
t89 = cos(pkin(8));
t156 = mrSges(4,3) * (t87 ^ 2 + t89 ^ 2);
t140 = cos(qJ(4));
t92 = sin(qJ(4));
t62 = t140 * t87 + t92 * t89;
t146 = t62 / 0.2e1;
t171 = -m(4) - m(3);
t170 = -m(5) - m(6);
t166 = Ifges(6,6) - Ifges(5,6);
t116 = t140 * t89;
t100 = -t92 * t87 + t116;
t55 = t100 * qJD(4);
t30 = qJD(1) * t55 + qJDD(1) * t62;
t16 = -qJDD(4) * mrSges(6,1) + t30 * mrSges(6,2);
t165 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t30 + t16;
t122 = qJDD(1) * t87;
t56 = t62 * qJD(4);
t31 = qJD(1) * t56 - qJDD(1) * t116 + t122 * t92;
t18 = -mrSges(6,2) * t31 + qJDD(4) * mrSges(6,3);
t164 = -qJDD(4) * mrSges(5,2) - mrSges(5,3) * t31 + t18;
t53 = t100 * qJD(1);
t138 = Ifges(6,5) * t53;
t49 = Ifges(5,4) * t53;
t54 = t62 * qJD(1);
t163 = t167 * qJD(4) + t169 * t54 - t138 + t49;
t134 = t53 * mrSges(5,3);
t135 = t53 * mrSges(6,2);
t44 = qJD(4) * mrSges(6,3) + t135;
t127 = -qJD(4) * mrSges(5,2) + t134 + t44;
t132 = t54 * mrSges(5,3);
t133 = t54 * mrSges(6,2);
t126 = t172 * qJD(4) - t132 - t133;
t75 = t89 * qJDD(2);
t39 = -t59 * t87 + t75;
t40 = t87 * qJDD(2) + t89 * t59;
t161 = -t39 * t87 + t40 * t89;
t86 = qJ(1) + pkin(7);
t79 = sin(t86);
t81 = cos(t86);
t160 = g(1) * t81 + g(2) * t79;
t85 = pkin(8) + qJ(4);
t78 = sin(t85);
t80 = cos(t85);
t107 = mrSges(5,1) * t80 - mrSges(5,2) * t78;
t108 = -mrSges(4,1) * t89 + mrSges(4,2) * t87;
t157 = -m(4) * pkin(2) - mrSges(3,1) - t107 + t108;
t35 = t75 + (-pkin(6) * qJDD(1) - t59) * t87;
t121 = qJDD(1) * t89;
t36 = pkin(6) * t121 + t40;
t125 = pkin(6) * qJD(1);
t66 = t69 * qJD(1);
t77 = t89 * qJD(2);
t37 = t77 + (-t66 - t125) * t87;
t46 = t87 * qJD(2) + t89 * t66;
t38 = t125 * t89 + t46;
t9 = t140 * t38 + t92 * t37;
t4 = -qJD(4) * t9 + t140 * t35 - t92 * t36;
t7 = qJD(4) * qJ(5) + t9;
t155 = -m(6) * t7 - t127;
t117 = t140 * t37;
t131 = t92 * t38;
t8 = t117 - t131;
t6 = -qJD(4) * pkin(4) + qJD(5) - t8;
t154 = -m(6) * t6 + t126;
t153 = -m(4) * qJ(3) + mrSges(3,2) - mrSges(6,2) - mrSges(4,3) - mrSges(5,3);
t151 = t53 / 0.2e1;
t150 = -t53 / 0.2e1;
t148 = t54 / 0.2e1;
t90 = cos(pkin(7));
t145 = pkin(1) * t90;
t93 = sin(qJ(1));
t144 = pkin(1) * t93;
t94 = cos(qJ(1));
t82 = t94 * pkin(1);
t141 = pkin(6) + t69;
t139 = Ifges(5,4) * t54;
t120 = m(4) - t170;
t118 = qJD(4) * t117 + t140 * t36 + t92 * t35;
t73 = -pkin(2) - t145;
t72 = pkin(3) * t89 + pkin(2);
t113 = t31 * mrSges(5,1) + t30 * mrSges(5,2);
t112 = t31 * mrSges(6,1) - t30 * mrSges(6,3);
t110 = -mrSges(4,1) * t121 + mrSges(4,2) * t122;
t105 = t80 * mrSges(6,1) + t78 * mrSges(6,3);
t104 = pkin(4) * t80 + qJ(5) * t78;
t103 = -(-t66 * t87 + t77) * t87 + t46 * t89;
t65 = -t72 - t145;
t57 = t141 * t87;
t58 = t141 * t89;
t101 = -t140 * t57 - t92 * t58;
t24 = t140 * t58 - t92 * t57;
t52 = qJD(1) * t65 + qJD(3);
t50 = qJDD(1) * t65 + qJDD(3);
t98 = m(6) * t104 + t105;
t91 = -pkin(6) - qJ(3);
t64 = qJDD(1) * t73 + qJDD(3);
t48 = Ifges(6,5) * t54;
t26 = -mrSges(6,1) * t53 - mrSges(6,3) * t54;
t25 = pkin(4) * t54 - qJ(5) * t53;
t20 = t53 * Ifges(5,2) + Ifges(5,6) * qJD(4) + t139;
t19 = Ifges(6,6) * qJD(4) - t53 * Ifges(6,3) + t48;
t14 = -pkin(4) * t100 - qJ(5) * t62 + t65;
t13 = pkin(4) * t56 - qJ(5) * t55 - qJD(5) * t62;
t12 = -pkin(4) * t53 - qJ(5) * t54 + t52;
t5 = pkin(4) * t31 - qJ(5) * t30 - qJD(5) * t54 + t50;
t3 = -qJD(4) * t131 + t118;
t2 = -qJDD(4) * pkin(4) + qJDD(5) - t4;
t1 = qJDD(4) * qJ(5) + (qJD(5) - t131) * qJD(4) + t118;
t10 = [(t1 * t100 + t2 * t62 + t55 * t6 - t56 * t7) * mrSges(6,2) + (t100 * t3 - t4 * t62 - t55 * t8 - t56 * t9) * mrSges(5,3) + t5 * (-mrSges(6,1) * t100 - mrSges(6,3) * t62) + (m(5) * t65 - mrSges(5,1) * t100 + mrSges(5,2) * t62) * t50 + (-t100 * t166 + t167 * t62) * qJDD(4) / 0.2e1 + (-t100 * t168 + t169 * t62) * t30 / 0.2e1 - t100 * (Ifges(6,5) * t30 + Ifges(6,6) * qJDD(4)) / 0.2e1 + t100 * (Ifges(5,4) * t30 + Ifges(5,6) * qJDD(4)) / 0.2e1 + m(6) * (t12 * t13 + t14 * t5) + (-mrSges(2,1) * t94 + mrSges(2,2) * t93 + t171 * t82 + t170 * (t81 * t72 - t79 * t91 + t82) + (-t98 + t157) * t81 + t153 * t79) * g(2) + (mrSges(2,1) * t93 + mrSges(2,2) * t94 - t171 * t144 + t170 * (-t81 * t91 - t144) + t153 * t81 + (-m(6) * (-t104 - t72) + t105 + m(5) * t72 - t157) * t79) * g(1) + (Ifges(4,4) * t87 + Ifges(4,2) * t89) * t121 + (Ifges(4,1) * t87 + Ifges(4,4) * t89) * t122 + t59 * t156 + (m(3) * (t88 ^ 2 + t90 ^ 2) * pkin(1) ^ 2 + Ifges(3,3) + Ifges(2,3) + 0.2e1 * (t90 * mrSges(3,1) - t88 * mrSges(3,2)) * pkin(1)) * qJDD(1) + (-(Ifges(6,3) + Ifges(5,2)) * t100 + 0.2e1 * t168 * t146) * t31 - (-m(5) * t4 + m(6) * t2 + t165) * t101 + (t166 * t56 + t167 * t55) * qJD(4) / 0.2e1 + (t167 * qJDD(4) + t169 * t30) * t146 + (t168 * t56 + t169 * t55) * t148 + t163 * t55 / 0.2e1 + (m(5) * t3 + m(6) * t1 + t164) * t24 + t161 * mrSges(4,3) + m(4) * (t103 * qJD(3) + t161 * t69 + t64 * t73) + (-m(5) * t8 - t154) * (qJD(3) * t62 + qJD(4) * t24) + (m(5) * t9 - t155) * (qJD(3) * t100 + qJD(4) * t101) + t14 * t112 + t64 * t108 + t73 * t110 + t13 * t26 + (Ifges(6,3) * t150 - Ifges(5,2) * t151 + t12 * mrSges(6,1) + t19 / 0.2e1 + t52 * mrSges(5,1) - t20 / 0.2e1) * t56 + (t52 * mrSges(5,2) - t12 * mrSges(6,3) + Ifges(5,4) * t151 + Ifges(6,5) * t150) * t55 + t65 * t113; m(3) * qJDD(2) + t164 * t62 - t165 * t100 - t126 * t56 + t127 * t55 + m(4) * (t39 * t89 + t40 * t87) + m(5) * (t100 * t4 + t3 * t62 + t55 * t9 - t56 * t8) + m(6) * (t1 * t62 - t100 * t2 + t55 * t7 + t56 * t6) + (-m(3) - t120) * g(3); t126 * t54 - t127 * t53 - qJD(1) ^ 2 * t156 + t110 + t112 + t113 + (-g(1) * t79 + g(2) * t81) * t120 + (-t53 * t7 - t54 * t6 + t5) * m(6) + (-t53 * t9 + t54 * t8 + t50) * m(5) + (-qJD(1) * t103 + t64) * m(4); -(t166 * t54 + t167 * t53) * qJD(4) / 0.2e1 - (t169 * t53 - t139 + t19 + t48) * t54 / 0.2e1 - t12 * (mrSges(6,1) * t54 - mrSges(6,3) * t53) - t52 * (mrSges(5,1) * t54 + mrSges(5,2) * t53) + (t134 + t155) * t8 + (Ifges(6,3) * t54 + t138) * t151 + t20 * t148 + t7 * t133 - t6 * t135 + ((-m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3)) * t80 + (m(6) * pkin(4) + t172) * t78) * t160 + t166 * t31 + t167 * t30 + (t132 + t154) * t9 + (-pkin(4) * t2 + qJ(5) * t1 + qJD(5) * t7 - t12 * t25) * m(6) + (Ifges(6,2) + Ifges(5,3)) * qJDD(4) + (-Ifges(5,2) * t54 + t163 + t49) * t150 + (-t107 - t98) * g(3) + qJD(5) * t44 - t25 * t26 - pkin(4) * t16 + qJ(5) * t18 + t1 * mrSges(6,3) - t2 * mrSges(6,1) - t3 * mrSges(5,2) + t4 * mrSges(5,1); -qJD(4) * t44 + t54 * t26 + (g(3) * t80 - t7 * qJD(4) + t12 * t54 - t160 * t78 + t2) * m(6) + t16;];
tau = t10;
