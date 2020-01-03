% Calculate vector of inverse dynamics joint torques for
% S4RPRP7
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRP7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP7_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP7_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_invdynJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP7_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP7_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP7_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:03
% EndTime: 2019-12-31 16:47:06
% DurationCPUTime: 1.99s
% Computational Cost: add. (530->185), mult. (1021->233), div. (0->0), fcn. (395->4), ass. (0->86)
t134 = Ifges(5,1) + Ifges(4,1);
t43 = sin(qJ(1));
t133 = g(1) * t43;
t42 = sin(qJ(3));
t44 = cos(qJ(3));
t85 = qJD(1) * qJD(3);
t25 = qJDD(1) * t42 + t44 * t85;
t10 = -mrSges(5,2) * t25 + qJDD(3) * mrSges(5,3);
t88 = t44 * qJD(1);
t77 = mrSges(5,2) * t88;
t123 = -mrSges(4,3) * t88 - t77 + (mrSges(4,1) + mrSges(5,1)) * qJD(3);
t24 = qJDD(1) * t44 - t42 * t85;
t89 = t42 * qJD(1);
t27 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t89;
t78 = mrSges(5,2) * t89;
t30 = qJD(3) * mrSges(5,3) - t78;
t8 = -qJDD(3) * mrSges(5,1) + t24 * mrSges(5,2);
t132 = (-qJDD(3) * mrSges(4,2) - mrSges(4,3) * t25 + t10) * t42 + (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t24 - t8) * t44 + (-t123 * t42 + (t27 + t30) * t44) * qJD(3);
t126 = Ifges(5,4) + Ifges(4,5);
t125 = Ifges(5,6) - Ifges(4,6);
t102 = Ifges(5,5) * t42;
t104 = Ifges(4,4) * t42;
t131 = t134 * t44 + t102 - t104;
t64 = t42 * mrSges(5,1) - t44 * mrSges(5,3);
t66 = mrSges(4,1) * t42 + t44 * mrSges(4,2);
t129 = t64 + t66;
t45 = cos(qJ(1));
t119 = g(2) * t45 - t133;
t128 = -m(4) - m(5);
t46 = -pkin(1) - pkin(5);
t33 = qJD(1) * t46 + qJD(2);
t100 = t33 * t42;
t17 = qJD(3) * qJ(4) + t100;
t32 = qJDD(1) * t46 + qJDD(2);
t93 = qJD(3) * t42;
t4 = t32 * t44 - t33 * t93;
t3 = -qJDD(3) * pkin(3) + qJDD(4) - t4;
t127 = -t17 * qJD(3) + t3;
t124 = t131 * qJD(1) + qJD(3) * t126;
t103 = Ifges(4,4) * t44;
t122 = t42 * (Ifges(5,3) * t44 - t102) + t44 * (-Ifges(4,1) * t42 - t103);
t121 = t125 * t44 - t126 * t42;
t84 = qJDD(1) * qJ(2);
t86 = qJD(1) * qJD(2);
t34 = t84 + t86;
t92 = qJD(3) * t44;
t5 = t42 * t32 + t33 * t92;
t68 = t4 * t44 + t42 * t5;
t99 = t33 * t44;
t118 = -t99 + qJD(4);
t113 = -mrSges(2,1) + mrSges(3,2) - mrSges(5,2) - mrSges(4,3);
t57 = pkin(3) * t42 - qJ(4) * t44;
t112 = -m(5) * t57 + mrSges(2,2) - mrSges(3,3) - t129;
t47 = qJD(1) ^ 2;
t110 = t44 / 0.2e1;
t2 = qJDD(3) * qJ(4) + qJD(3) * qJD(4) + t5;
t107 = t2 * t42;
t101 = Ifges(5,5) * t44;
t95 = t45 * pkin(1) + t43 * qJ(2);
t87 = qJDD(1) * mrSges(3,2);
t73 = (t34 + t86) * qJ(2);
t72 = -t85 / 0.2e1;
t70 = -qJD(4) * t44 + qJD(2);
t67 = mrSges(4,1) * t44 - mrSges(4,2) * t42;
t65 = t44 * mrSges(5,1) + t42 * mrSges(5,3);
t61 = -Ifges(4,2) * t42 + t103;
t58 = pkin(3) * t44 + qJ(4) * t42;
t11 = -qJD(3) * pkin(3) + t118;
t56 = t11 * t42 + t17 * t44;
t26 = qJ(2) + t57;
t12 = t26 * qJD(1);
t55 = t12 * t65;
t54 = t42 * (-Ifges(4,2) * t44 - t104);
t51 = qJ(2) * t67;
t50 = -t47 * qJ(2) + t119;
t48 = qJD(3) * t56 - t3 * t44 + t107;
t36 = -qJDD(1) * pkin(1) + qJDD(2);
t35 = Ifges(5,5) * t88;
t22 = t66 * qJD(1);
t21 = t64 * qJD(1);
t20 = t58 * qJD(1);
t14 = Ifges(4,6) * qJD(3) + qJD(1) * t61;
t13 = Ifges(5,6) * qJD(3) + Ifges(5,3) * t89 + t35;
t6 = qJD(3) * t58 + t70;
t1 = pkin(3) * t25 - qJ(4) * t24 + qJD(1) * t70 + t84;
t7 = [t51 * t85 + (t126 * qJDD(3) + t134 * t24) * t110 - t26 * t24 * mrSges(5,3) + qJ(2) * t24 * mrSges(4,2) - t68 * mrSges(4,3) + (t44 * (-Ifges(5,1) * t42 + t101) + t122) * t85 / 0.2e1 - t124 * t93 / 0.2e1 + (t125 * t42 + t126 * t44) * qJDD(3) / 0.2e1 + (-t14 / 0.2e1 + t13 / 0.2e1) * t92 + (t101 / 0.2e1 - t61 / 0.2e1 + qJ(2) * mrSges(4,1) + t26 * mrSges(5,1) + (Ifges(5,3) + Ifges(4,2) / 0.2e1) * t42 + (Ifges(5,5) - Ifges(4,4)) * t110) * t25 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + m(4) * t73 + m(5) * (t1 * t26 + t12 * t6) + (-m(3) * t95 + t128 * (t45 * pkin(5) + t95) + t113 * t45 + t112 * t43) * g(2) - pkin(1) * t87 + m(3) * (-pkin(1) * t36 + t73) + t1 * t64 + t131 * t24 / 0.2e1 + (m(4) * t68 + m(5) * t48 + t128 * t133 + t132) * t46 + t36 * mrSges(3,2) + qJD(2) * t22 + t6 * t21 + (t55 + t121 * qJD(3) / 0.2e1) * qJD(3) - t42 * (Ifges(4,4) * t24 + Ifges(4,6) * qJDD(3)) / 0.2e1 + (0.2e1 * mrSges(3,3) + t66) * t34 + t54 * t72 + (-t11 * t93 + t127 * t44 - t107) * mrSges(5,2) + ((m(3) * pkin(1) - t113) * t43 + ((-m(3) + t128) * qJ(2) + t112) * t45) * g(1) + t42 * (Ifges(5,5) * t24 + Ifges(5,6) * qJDD(3)) / 0.2e1; t87 - t47 * mrSges(3,3) + (-t21 - t22) * qJD(1) + (-qJD(1) * t12 + t119 + t48) * m(5) + (t50 + t68) * m(4) + (t36 + t50) * m(3) + t132; t129 * g(3) - (-Ifges(5,1) * t89 + t13 + t35) * t88 / 0.2e1 + t118 * t30 + (-t3 * pkin(3) + g(3) * t57 + t2 * qJ(4) + t17 * qJD(4) + t119 * t58 - t12 * t20 - t56 * t33) * m(5) + t121 * t72 + t123 * t100 + t124 * t89 / 0.2e1 + t125 * t25 + t126 * t24 + (Ifges(5,2) + Ifges(4,3)) * qJDD(3) - t27 * t99 + t14 * t88 / 0.2e1 - qJD(1) * t55 + qJ(4) * t10 - t20 * t21 + t4 * mrSges(4,1) - t5 * mrSges(4,2) - pkin(3) * t8 + t17 * t77 + t11 * t78 + t2 * mrSges(5,3) - t3 * mrSges(5,1) + (t54 / 0.2e1 - t51 - t122 / 0.2e1) * t47 + t119 * (t67 + t65); t21 * t88 - qJD(3) * t30 + (-g(3) * t42 - t119 * t44 + t12 * t88 + t127) * m(5) + t8;];
tau = t7;
