% Calculate vector of inverse dynamics joint torques for
% S4RPPP1
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:46
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S4RPPP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPP1_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPP1_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPP1_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:25
% EndTime: 2018-11-14 13:45:27
% DurationCPUTime: 1.43s
% Computational Cost: add. (564->209), mult. (1415->259), div. (0->0), fcn. (924->10), ass. (0->100)
t127 = mrSges(4,1) + mrSges(3,3);
t126 = mrSges(5,2) + mrSges(4,3);
t124 = m(4) + m(5);
t125 = t124 * qJ(3) - mrSges(3,2) + t126;
t75 = cos(pkin(4));
t123 = pkin(1) * t75;
t72 = sin(pkin(6));
t73 = sin(pkin(4));
t122 = t72 * t73;
t74 = cos(pkin(6));
t121 = t73 * t74;
t120 = -mrSges(3,1) - mrSges(5,3);
t46 = (qJ(2) * qJDD(1) + qJD(1) * qJD(2)) * t73;
t105 = qJDD(1) * t75;
t95 = pkin(1) * t105;
t14 = t74 * t46 + t72 * t95;
t118 = ((-mrSges(3,1) + mrSges(4,2)) * t75 + t127 * t122) * qJD(1);
t117 = (t126 * t75 + (mrSges(4,1) + mrSges(5,1)) * t121) * qJD(1);
t111 = qJD(1) * t73;
t93 = qJ(2) * t111;
t110 = qJD(1) * t75;
t96 = pkin(1) * t110;
t30 = t72 * t96 + t74 * t93;
t113 = qJ(2) * t73;
t35 = t74 * t113 + t72 * t123;
t106 = qJDD(1) * t73;
t91 = t74 * t106;
t116 = mrSges(5,1) * t91 + mrSges(5,2) * t105;
t92 = t72 * t106;
t115 = mrSges(4,1) * t92 + mrSges(4,2) * t105;
t76 = sin(qJ(1));
t77 = cos(qJ(1));
t114 = t77 * pkin(1) + t76 * t113;
t112 = qJ(3) * t75;
t109 = qJD(2) * t73;
t108 = qJD(4) * t75;
t48 = t72 * t93;
t107 = qJD(3) + t48;
t31 = t72 * t46;
t104 = qJDD(3) + t31;
t103 = m(3) + t124;
t102 = qJD(1) * qJD(3);
t71 = pkin(4) - pkin(6);
t64 = cos(t71) / 0.2e1;
t70 = pkin(4) + pkin(6);
t68 = cos(t70);
t101 = t64 + t68 / 0.2e1;
t63 = -sin(t71) / 0.2e1;
t67 = sin(t70);
t100 = t67 / 0.2e1 + t63;
t99 = 0.2e1 * t75;
t94 = -pkin(1) * t74 - pkin(2);
t57 = -pkin(1) * t106 + qJDD(2);
t90 = -m(3) * t57 - mrSges(3,2) * t92;
t89 = -t76 * pkin(1) + t77 * t113;
t88 = -qJ(3) * t72 - pkin(1);
t86 = t94 * t75;
t84 = m(5) * qJ(4) - mrSges(4,2) - t120;
t83 = -qJD(3) * t72 - qJD(4) * t74;
t82 = -pkin(2) * t74 + t88;
t81 = pkin(3) * t121 + t112;
t20 = -t77 * t101 + t72 * t76;
t22 = t76 * t101 + t77 * t72;
t80 = -g(1) * t22 - g(2) * t20 - g(3) * (-t67 / 0.2e1 + t63);
t26 = t82 * t73;
t79 = (-pkin(2) - qJ(4)) * t74 + t88;
t12 = t79 * t73;
t78 = pkin(3) * t122 + (-qJ(4) + t94) * t75;
t58 = t72 * t113;
t52 = mrSges(4,2) * t91;
t50 = mrSges(5,1) * t92;
t45 = qJD(3) * t75 + t74 * t109;
t44 = t72 * t109 - t108;
t38 = (mrSges(5,1) * t122 - t75 * mrSges(5,3)) * qJD(1);
t37 = (-t75 * mrSges(3,2) + mrSges(3,3) * t121) * qJD(1);
t34 = t74 * t123 - t58;
t33 = t83 * t73;
t29 = (-t72 * mrSges(5,2) - mrSges(5,3) * t74) * t111;
t28 = (mrSges(4,2) * t74 - t72 * mrSges(4,3)) * t111;
t27 = t74 * t96 - t48;
t25 = t58 + t86;
t24 = -t112 - t35;
t23 = -t76 * t100 + t74 * t77;
t21 = t77 * t100 + t76 * t74;
t17 = -qJ(3) * t110 - t30;
t16 = qJD(1) * t26 + qJD(2);
t15 = qJD(1) * t86 + t107;
t13 = t74 * t95 - t31;
t11 = t81 + t35;
t10 = t58 + t78;
t9 = qJDD(1) * t86 + t104;
t8 = qJD(1) * t12 + qJD(2);
t7 = t81 * qJD(1) + qJD(4) + t30;
t6 = qJDD(2) + (t82 * qJDD(1) - t72 * t102) * t73;
t5 = (-qJ(3) * qJDD(1) - t102) * t75 - t14;
t4 = t78 * qJD(1) + t107;
t3 = t81 * qJDD(1) + t75 * t102 + qJDD(4) + t14;
t2 = qJDD(2) + (t83 * qJD(1) + t79 * qJDD(1)) * t73;
t1 = -qJD(1) * t108 + t78 * qJDD(1) + t104;
t18 = [t10 * t50 + t26 * t52 + t11 * t116 + t25 * t115 + t33 * t29 + t44 * t38 + t117 * t45 + m(3) * (t13 * t34 + t14 * t35) + m(4) * (-t17 * t45 + t24 * t5 + t25 * t9 + t26 * t6) + m(5) * (t1 * t10 + t11 * t3 + t12 * t2 + t33 * t8 + t4 * t44 + t45 * t7) + (-m(3) * t114 - t77 * mrSges(2,1) + t76 * mrSges(2,2) - t84 * t23 - t125 * t22 - t124 * (t23 * pkin(2) + t114)) * g(2) + (-m(3) * t89 + t76 * mrSges(2,1) + t77 * mrSges(2,2) + t84 * t21 + t125 * t20 + t124 * (t21 * pkin(2) - t89)) * g(1) + (t13 * mrSges(3,1) - t14 * mrSges(3,2) + t9 * mrSges(4,2) + t3 * mrSges(5,2) - t5 * mrSges(4,3) - t1 * mrSges(5,3)) * t75 + (t90 * pkin(1) + (-t57 * mrSges(3,1) - t5 * mrSges(4,1) + t3 * mrSges(5,1) + t6 * mrSges(4,2) + t14 * mrSges(3,3) - t2 * mrSges(5,3) + (m(3) * t30 + t37) * qJD(2)) * t74 + (t9 * mrSges(4,1) + t1 * mrSges(5,1) + t57 * mrSges(3,2) - t2 * mrSges(5,2) - t13 * mrSges(3,3) - t6 * mrSges(4,3) + (-m(4) * t16 - t28) * qJD(3) + (-m(3) * t27 + m(4) * t15 + t118) * qJD(2)) * t72 + (g(1) * t77 + g(2) * t76) * (-m(5) * pkin(3) - mrSges(5,1) - t127)) * t73 + (Ifges(2,3) + (t34 * mrSges(3,1) - t35 * mrSges(3,2) - t24 * mrSges(4,3) - t10 * mrSges(5,3) + (Ifges(3,3) + Ifges(5,1) + Ifges(4,1)) * t75) * t75 + ((-t12 * mrSges(5,2) - t34 * mrSges(3,3) - t26 * mrSges(4,3) + (Ifges(5,3) + Ifges(4,2) + Ifges(3,1)) * t122 + (-Ifges(4,4) + Ifges(3,5) + Ifges(5,5)) * t99) * t72 + (-t24 * mrSges(4,1) + t35 * mrSges(3,3) - t12 * mrSges(5,3) + (pkin(1) * mrSges(3,1) + (Ifges(5,2) + Ifges(4,3) + Ifges(3,2)) * t74) * t73 + 0.2e1 * (Ifges(3,4) + Ifges(4,6) - Ifges(5,6)) * t122 + (Ifges(3,6) - Ifges(5,4) - Ifges(4,5)) * t99) * t74) * t73) * qJDD(1); t52 - t103 * t75 * g(3) + m(4) * t6 + m(5) * t2 + ((t120 * t74 - t126 * t72) * qJDD(1) + ((-t37 - t117) * t74 + (-t38 - t118) * t72 - m(5) * (t4 * t72 + t7 * t74) - m(4) * (t15 * t72 - t17 * t74) - m(3) * (-t27 * t72 + t30 * t74)) * qJD(1) + (-g(1) * t76 + g(2) * t77) * t103) * t73 - t90; -mrSges(5,3) * t105 + t50 + (t1 + t80) * m(5) + (t80 + t9) * m(4) + (-t117 * t75 + (t28 + t29) * t122 - m(4) * (-t16 * t122 - t17 * t75) - m(5) * (-t8 * t122 + t7 * t75)) * qJD(1) + t115; (t29 * t121 + t75 * t38) * qJD(1) + (t3 - g(1) * t23 - g(2) * t21 - g(3) * (t64 - t68 / 0.2e1) - (-t8 * t121 - t4 * t75) * qJD(1)) * m(5) + t116;];
tau  = t18;
