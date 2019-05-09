% Calculate vector of cutting forces with Newton-Euler
% S6RPPRPR2
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x7]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 14:04
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPPRPR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:01:49
% EndTime: 2019-05-05 14:01:52
% DurationCPUTime: 1.43s
% Computational Cost: add. (14765->165), mult. (33081->203), div. (0->0), fcn. (21816->10), ass. (0->88)
t86 = sin(qJ(1));
t88 = cos(qJ(1));
t107 = t86 * g(1) - t88 * g(2);
t67 = qJDD(1) * pkin(1) + t107;
t105 = -t88 * g(1) - t86 * g(2);
t90 = qJD(1) ^ 2;
t68 = -t90 * pkin(1) + t105;
t81 = sin(pkin(9));
t83 = cos(pkin(9));
t106 = t83 * t67 - t81 * t68;
t103 = qJDD(3) - t106;
t124 = cos(qJ(4));
t80 = sin(pkin(10));
t82 = cos(pkin(10));
t85 = sin(qJ(4));
t128 = -t82 * t124 + t80 * t85;
t61 = t128 * qJD(1);
t113 = t61 * qJD(4);
t114 = pkin(7) * qJDD(1);
t111 = qJD(1) * qJD(3);
t79 = -g(3) + qJDD(2);
t116 = -0.2e1 * t80 * t111 + t82 * t79;
t125 = pkin(3) * t90;
t117 = t81 * t67 + t83 * t68;
t39 = -t90 * pkin(2) + qJDD(1) * qJ(3) + t117;
t27 = (t82 * t125 - t114 - t39) * t80 + t116;
t109 = t80 * t79 + (0.2e1 * t111 + t39) * t82;
t78 = t82 ^ 2;
t30 = t82 * t114 - t78 * t125 + t109;
t104 = t124 * t27 - t85 * t30;
t99 = t124 * t80 + t82 * t85;
t62 = t99 * qJD(1);
t42 = t61 * pkin(4) - t62 * qJ(5);
t89 = qJD(4) ^ 2;
t19 = -qJDD(4) * pkin(4) - t89 * qJ(5) + t62 * t42 + qJDD(5) - t104;
t49 = t99 * qJDD(1) - t113;
t14 = (t61 * t62 - qJDD(4)) * pkin(8) + (t49 + t113) * pkin(5) + t19;
t112 = t62 * qJD(4);
t48 = t128 * qJDD(1) + t112;
t56 = t62 * pkin(5) - qJD(4) * pkin(8);
t60 = t61 ^ 2;
t127 = -2 * qJD(5);
t115 = t80 ^ 2 + t78;
t92 = (-pkin(3) * t82 - pkin(2)) * qJDD(1) + (-t115 * pkin(7) - qJ(3)) * t90 + t103;
t91 = pkin(4) * t112 + t62 * t127 + (-t49 + t113) * qJ(5) + t92;
t17 = -t60 * pkin(5) - t62 * t56 + (pkin(4) + pkin(8)) * t48 + t91;
t84 = sin(qJ(6));
t87 = cos(qJ(6));
t50 = -t84 * qJD(4) + t87 * t61;
t29 = t50 * qJD(6) + t87 * qJDD(4) + t84 * t48;
t51 = t87 * qJD(4) + t84 * t61;
t33 = -t50 * mrSges(7,1) + t51 * mrSges(7,2);
t59 = qJD(6) + t62;
t35 = -t59 * mrSges(7,2) + t50 * mrSges(7,3);
t47 = qJDD(6) + t49;
t12 = m(7) * (t87 * t14 - t84 * t17) - t29 * mrSges(7,3) + t47 * mrSges(7,1) - t51 * t33 + t59 * t35;
t28 = -t51 * qJD(6) - t84 * qJDD(4) + t87 * t48;
t36 = t59 * mrSges(7,1) - t51 * mrSges(7,3);
t13 = m(7) * (t84 * t14 + t87 * t17) + t28 * mrSges(7,3) - t47 * mrSges(7,2) + t50 * t33 - t59 * t36;
t55 = t62 * mrSges(6,1) + qJD(4) * mrSges(6,2);
t100 = t84 * t12 - t87 * t13 - m(6) * (t48 * pkin(4) + t91) + t62 * t55 + t49 * mrSges(6,3);
t54 = t61 * mrSges(6,1) - qJD(4) * mrSges(6,3);
t118 = -qJD(4) * mrSges(5,2) - t61 * mrSges(5,3) - t54;
t122 = mrSges(5,1) - mrSges(6,2);
t53 = qJD(4) * mrSges(5,1) - t62 * mrSges(5,3);
t93 = m(5) * t92 + t49 * mrSges(5,2) + t118 * t61 + t122 * t48 + t62 * t53 - t100;
t130 = -m(4) * (-qJDD(1) * pkin(2) - t90 * qJ(3) + t103) - t93;
t129 = t115 * mrSges(4,3);
t121 = -mrSges(5,3) - mrSges(6,1);
t120 = t124 * t30 + t85 * t27;
t44 = -t61 * mrSges(6,2) - t62 * mrSges(6,3);
t119 = -t61 * mrSges(5,1) - t62 * mrSges(5,2) - t44;
t95 = -t89 * pkin(4) + qJDD(4) * qJ(5) - t61 * t42 + t120;
t97 = -t28 * mrSges(7,1) - t50 * t35 + m(7) * (-t48 * pkin(5) - t60 * pkin(8) + ((2 * qJD(5)) + t56) * qJD(4) + t95) + t29 * mrSges(7,2) + t51 * t36;
t94 = -m(6) * (qJD(4) * t127 - t95) + t97;
t10 = m(5) * t120 + t119 * t61 + t121 * t48 + (-mrSges(5,2) + mrSges(6,3)) * qJDD(4) + (-t53 + t55) * qJD(4) + t94;
t102 = -t82 * mrSges(4,1) + t80 * mrSges(4,2);
t101 = qJDD(1) * mrSges(4,3) + t90 * t102;
t98 = -m(6) * t19 - t87 * t12 - t84 * t13;
t9 = m(5) * t104 + t118 * qJD(4) + t122 * qJDD(4) + t119 * t62 + t121 * t49 + t98;
t6 = m(4) * t116 + t85 * t10 + t124 * t9 + (-m(4) * t39 - t101) * t80;
t7 = m(4) * t109 + t124 * t10 + t101 * t82 - t85 * t9;
t110 = m(3) * t79 + t82 * t6 + t80 * t7;
t8 = m(3) * t106 + (-mrSges(3,2) + t129) * t90 + (mrSges(3,1) - t102) * qJDD(1) + t130;
t3 = m(3) * t117 - t90 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t80 * t6 + t82 * t7;
t2 = m(2) * t105 - t90 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t83 * t3 - t81 * t8;
t1 = m(2) * t107 + qJDD(1) * mrSges(2,1) - t90 * mrSges(2,2) + t81 * t3 + t83 * t8;
t4 = [-m(1) * g(1) - t86 * t1 + t88 * t2, t2, t3, t7, t10, -t48 * mrSges(6,2) - t61 * t54 - t100, t13; -m(1) * g(2) + t88 * t1 + t86 * t2, t1, t8, t6, t9, t48 * mrSges(6,1) - qJDD(4) * mrSges(6,3) - qJD(4) * t55 + t61 * t44 - t94, t12; (-m(1) - m(2)) * g(3) + t110, -m(2) * g(3) + t110, t110, t102 * qJDD(1) - t90 * t129 - t130, t93, t49 * mrSges(6,1) + qJDD(4) * mrSges(6,2) + qJD(4) * t54 + t62 * t44 - t98, t97;];
f_new  = t4;
