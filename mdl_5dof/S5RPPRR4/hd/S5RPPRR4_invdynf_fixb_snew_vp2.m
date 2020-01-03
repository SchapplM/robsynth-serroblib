% Calculate vector of cutting forces with Newton-Euler
% S5RPPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% f_new [3x6]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPPRR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR4_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR4_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR4_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:30:35
% EndTime: 2020-01-03 11:30:42
% DurationCPUTime: 2.04s
% Computational Cost: add. (17152->147), mult. (46826->210), div. (0->0), fcn. (31990->10), ass. (0->89)
t79 = sin(pkin(8));
t81 = cos(pkin(8));
t101 = -pkin(2) * t81 - qJ(3) * t79;
t113 = qJD(1) * t79;
t84 = sin(qJ(1));
t87 = cos(qJ(1));
t115 = -t87 * g(2) - t84 * g(3);
t88 = qJD(1) ^ 2;
t96 = -t88 * qJ(2) + qJDD(2) - t115;
t123 = -0.2e1 * qJD(3) * t113 + (-pkin(1) + t101) * qJDD(1) + t96;
t76 = t79 ^ 2;
t122 = (t81 ^ 2 + t76) * mrSges(3,3);
t78 = sin(pkin(9));
t120 = t78 * t79;
t97 = t81 * mrSges(4,2) - mrSges(4,3) * t120;
t56 = t97 * qJD(1);
t80 = cos(pkin(9));
t119 = t79 * t80;
t98 = -t81 * mrSges(4,1) - mrSges(4,3) * t119;
t57 = t98 * qJD(1);
t100 = t56 * t78 + t57 * t80;
t102 = mrSges(4,1) * t78 + mrSges(4,2) * t80;
t107 = -t84 * g(2) + t87 * g(3);
t58 = -t88 * pkin(1) + qJDD(1) * qJ(2) + t107;
t116 = -t81 * g(1) - t79 * t58;
t103 = -t81 * mrSges(3,1) + t79 * mrSges(3,2);
t61 = t103 * qJD(1);
t109 = 0.2e1 * qJD(1) * qJD(2);
t60 = t101 * qJD(1);
t35 = t79 * t109 + t60 * t113 + qJDD(3) - t116;
t83 = sin(qJ(4));
t86 = cos(qJ(4));
t94 = (-t78 * t83 + t80 * t86) * t79;
t50 = qJD(1) * t94;
t95 = (-t78 * t86 - t80 * t83) * t79;
t39 = -t50 * qJD(4) + qJDD(1) * t95;
t49 = qJD(1) * t95;
t40 = t49 * qJD(4) + qJDD(1) * t94;
t112 = t81 * qJD(1);
t68 = qJD(4) - t112;
t41 = -t68 * mrSges(5,2) + t49 * mrSges(5,3);
t42 = t68 * mrSges(5,1) - t50 * mrSges(5,3);
t108 = qJDD(1) * t120;
t111 = t78 ^ 2 * t76 * t88;
t99 = -pkin(3) * t81 - pkin(6) * t119;
t59 = t99 * qJD(1);
t91 = t80 * t59 * t113 + pkin(3) * t108 - pkin(6) * t111 + t35;
t82 = sin(qJ(5));
t85 = cos(qJ(5));
t34 = t82 * t49 + t85 * t50;
t18 = -t34 * qJD(5) + t85 * t39 - t82 * t40;
t33 = t85 * t49 - t82 * t50;
t19 = t33 * qJD(5) + t82 * t39 + t85 * t40;
t66 = qJD(5) + t68;
t29 = -t66 * mrSges(6,2) + t33 * mrSges(6,3);
t30 = t66 * mrSges(6,1) - t34 * mrSges(6,3);
t45 = t68 * pkin(4) - t50 * pkin(7);
t48 = t49 ^ 2;
t92 = t18 * mrSges(6,1) + t33 * t29 - m(6) * (-t39 * pkin(4) - t48 * pkin(7) + t50 * t45 + t91) - t19 * mrSges(6,2) - t34 * t30;
t90 = -m(5) * t91 + t39 * mrSges(5,1) - t40 * mrSges(5,2) + t49 * t41 - t50 * t42 + t92;
t89 = m(4) * t35 - t90;
t10 = -t89 + ((-mrSges(3,3) - t102) * qJDD(1) + (-0.2e1 * m(3) * qJD(2) - t100 - t61) * qJD(1)) * t79 + m(3) * t116;
t104 = -t79 * g(1) + (t109 + t58) * t81;
t117 = t123 * t80;
t36 = t60 * t112 + t104;
t53 = t102 * t113;
t23 = t99 * qJDD(1) + (-t36 + (-pkin(3) * t76 * t80 + pkin(6) * t79 * t81) * t88) * t78 + t117;
t110 = t123 * t78 + t80 * t36;
t24 = -pkin(3) * t111 - pkin(6) * t108 + t59 * t112 + t110;
t106 = t86 * t23 - t83 * t24;
t67 = -t81 * qJDD(1) + qJDD(4);
t13 = (t49 * t68 - t40) * pkin(7) + (t49 * t50 + t67) * pkin(4) + t106;
t118 = t83 * t23 + t86 * t24;
t14 = -t48 * pkin(4) + t39 * pkin(7) - t68 * t45 + t118;
t26 = -t33 * mrSges(6,1) + t34 * mrSges(6,2);
t64 = qJDD(5) + t67;
t11 = m(6) * (t85 * t13 - t82 * t14) - t19 * mrSges(6,3) + t64 * mrSges(6,1) - t34 * t26 + t66 * t29;
t12 = m(6) * (t82 * t13 + t85 * t14) + t18 * mrSges(6,3) - t64 * mrSges(6,2) + t33 * t26 - t66 * t30;
t37 = -t49 * mrSges(5,1) + t50 * mrSges(5,2);
t7 = m(5) * t106 + t67 * mrSges(5,1) - t40 * mrSges(5,3) + t85 * t11 + t82 * t12 - t50 * t37 + t68 * t41;
t8 = m(5) * t118 - t67 * mrSges(5,2) + t39 * mrSges(5,3) - t82 * t11 + t85 * t12 + t49 * t37 - t68 * t42;
t5 = m(4) * (-t78 * t36 + t117) + t83 * t8 + t86 * t7 + t98 * qJDD(1) + (-t53 * t119 - t81 * t56) * qJD(1);
t6 = m(4) * t110 + t86 * t8 - t83 * t7 + t97 * qJDD(1) + (-t53 * t120 + t81 * t57) * qJD(1);
t4 = m(3) * t104 + t80 * t6 - t78 * t5 + (qJDD(1) * mrSges(3,3) + qJD(1) * t61) * t81;
t121 = t81 * t10 + t79 * t4;
t93 = m(3) * (-qJDD(1) * pkin(1) + t96) + t80 * t5 + t78 * t6;
t2 = m(2) * t115 + (-mrSges(2,2) + t122) * t88 + (mrSges(2,1) - t103) * qJDD(1) - t93;
t1 = m(2) * t107 - t88 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t79 * t10 + t81 * t4;
t3 = [(-m(1) - m(2)) * g(1) + t121, t1, t4, t6, t8, t12; -m(1) * g(2) + t84 * t1 + t87 * t2, t2, t10, t5, t7, t11; -m(1) * g(3) - t87 * t1 + t84 * t2, -m(2) * g(1) + t121, t103 * qJDD(1) - t88 * t122 + t93, (t100 * qJD(1) + t102 * qJDD(1)) * t79 + t89, -t90, -t92;];
f_new = t3;
