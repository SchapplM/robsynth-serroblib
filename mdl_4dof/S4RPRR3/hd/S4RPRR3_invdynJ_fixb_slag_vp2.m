% Calculate vector of inverse dynamics joint torques for
% S4RPRR3
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR3_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR3_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR3_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR3_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR3_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:05
% EndTime: 2019-12-31 16:49:10
% DurationCPUTime: 2.37s
% Computational Cost: add. (1244->234), mult. (2674->337), div. (0->0), fcn. (1585->12), ass. (0->112)
t85 = qJ(3) + qJ(4);
t79 = sin(t85);
t80 = cos(t85);
t108 = t80 * mrSges(5,1) - mrSges(5,2) * t79;
t89 = sin(qJ(3));
t92 = cos(qJ(3));
t66 = -mrSges(4,1) * t92 + mrSges(4,2) * t89;
t148 = -t108 + t66;
t84 = qJ(1) + pkin(7);
t74 = sin(t84);
t75 = cos(t84);
t143 = g(1) * t75 + g(2) * t74;
t114 = qJD(1) * t92;
t147 = -Ifges(4,5) * qJD(3) / 0.2e1 - Ifges(4,4) * t114 / 0.2e1;
t142 = m(3) + m(5) + m(4);
t146 = pkin(1) * t142 + mrSges(2,1);
t128 = Ifges(4,4) * t89;
t145 = t128 / 0.2e1;
t112 = qJD(3) * t89;
t86 = sin(pkin(7));
t68 = pkin(1) * t86 + pkin(5);
t60 = t68 * qJDD(1);
t62 = t68 * qJD(1);
t78 = t92 * qJD(2);
t26 = qJD(3) * t78 + t89 * qJDD(2) - t112 * t62 + t92 * t60;
t113 = qJD(2) * t89;
t41 = t62 * t92 + t113;
t27 = -t41 * qJD(3) + t92 * qJDD(2) - t60 * t89;
t144 = t26 * t92 - t27 * t89;
t102 = Ifges(4,2) * t92 + t128;
t117 = Ifges(4,6) * qJD(3);
t141 = -t41 * mrSges(4,3) - qJD(1) * t102 / 0.2e1 - t117 / 0.2e1;
t115 = qJD(1) * t89;
t40 = -t62 * t89 + t78;
t140 = -t40 * mrSges(4,3) + Ifges(4,1) * t115 / 0.2e1 - t147;
t71 = pkin(3) * t92 + pkin(2);
t139 = -m(4) * pkin(2) - m(5) * t71 - mrSges(3,1) + t148;
t138 = -m(4) * pkin(5) + m(5) * (-pkin(6) - pkin(5)) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t88 = sin(qJ(4));
t91 = cos(qJ(4));
t56 = t88 * t92 + t89 * t91;
t50 = t56 * qJD(1);
t136 = t50 / 0.2e1;
t87 = cos(pkin(7));
t135 = pkin(1) * t87;
t131 = pkin(6) + t68;
t55 = -t88 * t89 + t91 * t92;
t49 = t55 * qJD(1);
t129 = mrSges(5,3) * t49;
t127 = Ifges(4,4) * t92;
t106 = pkin(6) * qJD(1) + t62;
t36 = t106 * t92 + t113;
t124 = t36 * t88;
t123 = t36 * t91;
t120 = t50 * mrSges(5,3);
t119 = t50 * Ifges(5,4);
t111 = qJD(3) * t92;
t110 = qJD(1) * qJD(3);
t109 = pkin(3) * t112;
t69 = -pkin(2) - t135;
t107 = qJD(3) * t131;
t104 = mrSges(4,1) * t89 + mrSges(4,2) * t92;
t103 = mrSges(5,1) * t79 + mrSges(5,2) * t80;
t35 = -t106 * t89 + t78;
t34 = qJD(3) * pkin(3) + t35;
t7 = t34 * t91 - t124;
t8 = t34 * t88 + t123;
t52 = t131 * t89;
t53 = t131 * t92;
t28 = -t52 * t91 - t53 * t88;
t29 = -t52 * t88 + t53 * t91;
t61 = t69 * qJDD(1);
t59 = -t71 - t135;
t101 = t56 * qJD(4);
t100 = t55 * qJD(4);
t58 = qJDD(1) * t89 + t110 * t92;
t13 = qJDD(3) * pkin(3) - pkin(6) * t58 + t27;
t57 = qJDD(1) * t92 - t110 * t89;
t16 = pkin(6) * t57 + t26;
t2 = qJD(4) * t7 + t13 * t88 + t16 * t91;
t21 = qJD(1) * t100 + t57 * t88 + t58 * t91;
t22 = -qJD(1) * t101 + t57 * t91 - t58 * t88;
t83 = qJD(3) + qJD(4);
t24 = t49 * Ifges(5,2) + t83 * Ifges(5,6) + t119;
t46 = Ifges(5,4) * t49;
t25 = t50 * Ifges(5,1) + t83 * Ifges(5,5) + t46;
t3 = -qJD(4) * t8 + t13 * t91 - t16 * t88;
t51 = t59 * qJD(1);
t82 = qJDD(3) + qJDD(4);
t96 = t3 * mrSges(5,1) - t2 * mrSges(5,2) + t7 * t129 + t24 * t136 - t51 * (mrSges(5,1) * t50 + mrSges(5,2) * t49) - t50 * (Ifges(5,1) * t49 - t119) / 0.2e1 + Ifges(5,6) * t22 + Ifges(5,5) * t21 - t83 * (Ifges(5,5) * t49 - Ifges(5,6) * t50) / 0.2e1 + Ifges(5,3) * t82 - (-Ifges(5,2) * t50 + t25 + t46) * t49 / 0.2e1;
t93 = cos(qJ(1));
t90 = sin(qJ(1));
t65 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t114;
t64 = t69 * qJD(1);
t63 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t115;
t45 = t92 * t107;
t44 = t89 * t107;
t43 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t58;
t42 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t57;
t39 = mrSges(5,1) * t83 - t120;
t38 = -mrSges(5,2) * t83 + t129;
t37 = -pkin(3) * t57 + t61;
t33 = -qJD(3) * t56 - t101;
t32 = qJD(3) * t55 + t100;
t31 = -mrSges(5,1) * t49 + mrSges(5,2) * t50;
t15 = t35 * t91 - t124;
t14 = -t35 * t88 - t123;
t12 = -mrSges(5,2) * t82 + mrSges(5,3) * t22;
t11 = mrSges(5,1) * t82 - mrSges(5,3) * t21;
t5 = -qJD(4) * t29 + t44 * t88 - t45 * t91;
t4 = qJD(4) * t28 - t44 * t91 - t45 * t88;
t1 = [(mrSges(2,2) * t93 + t138 * t75 - t139 * t74 + t146 * t90) * g(1) + (mrSges(2,2) * t90 + t138 * t74 + t139 * t75 - t146 * t93) * g(2) + (-t32 * t7 + t33 * t8) * mrSges(5,3) + t57 * t102 / 0.2e1 + (mrSges(5,2) * t37 - mrSges(5,3) * t3 + Ifges(5,1) * t21 + Ifges(5,4) * t22 + Ifges(5,5) * t82) * t56 + t58 * Ifges(4,1) * t89 + (-mrSges(5,1) * t37 + mrSges(5,3) * t2 + Ifges(5,4) * t21 + Ifges(5,2) * t22 + Ifges(5,6) * t82) * t55 + (Ifges(5,1) * t32 + Ifges(5,4) * t33) * t136 + (t92 * (-Ifges(4,2) * t89 + t127) + t89 * (Ifges(4,1) * t92 - t128)) * t110 / 0.2e1 + qJD(3) ^ 2 * (Ifges(4,5) * t92 - Ifges(4,6) * t89) / 0.2e1 + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t87 - 0.2e1 * mrSges(3,2) * t86 + m(3) * (t86 ^ 2 + t87 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + qJDD(3) * (Ifges(4,5) * t89 + Ifges(4,6) * t92) + t83 * (Ifges(5,5) * t32 + Ifges(5,6) * t33) / 0.2e1 + (t92 * t42 - t89 * t43 + m(4) * ((-t40 * t92 - t41 * t89) * qJD(3) + t144) - t63 * t111 - t65 * t112) * t68 + t144 * mrSges(4,3) + t140 * t111 + t141 * t112 + (m(4) * t69 + t66) * t61 + t69 * (-mrSges(4,1) * t57 + mrSges(4,2) * t58) + t59 * (-mrSges(5,1) * t22 + mrSges(5,2) * t21) + t49 * (Ifges(5,4) * t32 + Ifges(5,2) * t33) / 0.2e1 + t51 * (-mrSges(5,1) * t33 + mrSges(5,2) * t32) + t4 * t38 + t5 * t39 + t28 * t11 + t29 * t12 + t32 * t25 / 0.2e1 + t33 * t24 / 0.2e1 + t31 * t109 + m(5) * (t51 * t109 + t2 * t29 + t28 * t3 + t37 * t59 + t4 * t8 + t5 * t7) + t64 * t104 * qJD(3) + t57 * t145 + t92 * (Ifges(4,4) * t58 + Ifges(4,2) * t57) / 0.2e1 + t58 * t127 / 0.2e1; m(3) * qJDD(2) + t55 * t11 + t56 * t12 + t32 * t38 + t33 * t39 + t89 * t42 + t92 * t43 + (-t89 * t63 + t92 * t65) * qJD(3) + m(4) * (t26 * t89 + t27 * t92 + (-t40 * t89 + t41 * t92) * qJD(3)) + m(5) * (t2 * t56 + t3 * t55 + t32 * t8 + t33 * t7) - t142 * g(3); t148 * g(3) - m(5) * (t14 * t7 + t15 * t8) + t96 - t40 * t65 + Ifges(4,6) * t57 + Ifges(4,5) * t58 + t41 * t63 - t15 * t38 - t14 * t39 - t26 * mrSges(4,2) + t27 * mrSges(4,1) + ((-t64 * mrSges(4,2) - t140 + t147) * t92 + (-t64 * mrSges(4,1) + t117 / 0.2e1 + (t145 + (Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t92) * qJD(1) + (-m(5) * t51 - t31) * pkin(3) - t141) * t89) * qJD(1) + (t91 * t11 + t88 * t12 + (-g(3) * t92 + t143 * t89 + t2 * t88 + t3 * t91) * m(5) + (t38 * t91 - t39 * t88 + (-t7 * t88 + t8 * t91) * m(5)) * qJD(4)) * pkin(3) + Ifges(4,3) * qJDD(3) + t8 * t120 + t143 * (t103 + t104); (t39 + t120) * t8 - g(3) * t108 + t96 - t7 * t38 + t143 * t103;];
tau = t1;
