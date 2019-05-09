% Calculate vector of cutting forces with Newton-Euler
% S6PRRRPR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-05-05 07:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRRPR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:58:00
% EndTime: 2019-05-05 06:58:10
% DurationCPUTime: 5.29s
% Computational Cost: add. (74181->179), mult. (158114->245), div. (0->0), fcn. (117527->14), ass. (0->100)
t101 = sin(qJ(2));
t105 = cos(qJ(2));
t93 = sin(pkin(11));
t96 = cos(pkin(11));
t80 = t93 * g(1) - t96 * g(2);
t97 = cos(pkin(6));
t130 = t80 * t97;
t81 = -t96 * g(1) - t93 * g(2);
t91 = -g(3) + qJDD(1);
t94 = sin(pkin(6));
t133 = -t101 * t81 + (t91 * t94 + t130) * t105;
t100 = sin(qJ(3));
t104 = cos(qJ(3));
t106 = qJD(2) ^ 2;
t112 = -qJDD(2) * pkin(2) - t133;
t125 = qJD(2) * t100;
t123 = qJD(2) * qJD(3);
t79 = t104 * qJDD(2) - t100 * t123;
t85 = qJD(3) * pkin(3) - pkin(9) * t125;
t90 = t104 ^ 2;
t109 = -t79 * pkin(3) + t85 * t125 + (-pkin(9) * t90 - pkin(8)) * t106 + t112;
t102 = cos(qJ(6));
t103 = cos(qJ(4));
t99 = sin(qJ(4));
t73 = (t100 * t103 + t104 * t99) * qJD(2);
t118 = t104 * t123;
t78 = t100 * qJDD(2) + t118;
t51 = -t73 * qJD(4) + t103 * t79 - t99 * t78;
t89 = qJD(3) + qJD(4);
t67 = t89 * pkin(4) - t73 * qJ(5);
t72 = (-t100 * t99 + t103 * t104) * qJD(2);
t71 = t72 ^ 2;
t107 = -t51 * pkin(4) - t71 * qJ(5) + t73 * t67 + qJDD(5) + t109;
t131 = 2 * qJD(5);
t127 = t101 * t94;
t121 = t101 * t130 + t105 * t81 + t91 * t127;
t56 = -t106 * pkin(2) + qJDD(2) * pkin(8) + t121;
t69 = -t94 * t80 + t97 * t91;
t117 = -t100 * t56 + t104 * t69;
t35 = (-t78 + t118) * pkin(9) + (t100 * t104 * t106 + qJDD(3)) * pkin(3) + t117;
t128 = t100 * t69 + t104 * t56;
t38 = -t90 * t106 * pkin(3) + t79 * pkin(9) - qJD(3) * t85 + t128;
t119 = t103 * t35 - t99 * t38;
t52 = t72 * qJD(4) + t103 * t78 + t99 * t79;
t88 = qJDD(3) + qJDD(4);
t23 = (t72 * t89 - t52) * qJ(5) + (t72 * t73 + t88) * pkin(4) + t119;
t129 = t103 * t38 + t99 * t35;
t25 = -t71 * pkin(4) + t51 * qJ(5) - t67 * t89 + t129;
t92 = sin(pkin(12));
t95 = cos(pkin(12));
t61 = t95 * t72 - t92 * t73;
t122 = t61 * t131 + t92 * t23 + t95 * t25;
t62 = t92 * t72 + t95 * t73;
t44 = -t61 * pkin(5) - t62 * pkin(10);
t87 = t89 ^ 2;
t20 = -pkin(5) * t87 + pkin(10) * t88 + t61 * t44 + t122;
t36 = t95 * t51 - t92 * t52;
t37 = t92 * t51 + t95 * t52;
t21 = (-t61 * t89 - t37) * pkin(10) + (t62 * t89 - t36) * pkin(5) + t107;
t98 = sin(qJ(6));
t48 = t102 * t89 - t98 * t62;
t29 = t48 * qJD(6) + t102 * t37 + t98 * t88;
t34 = qJDD(6) - t36;
t49 = t102 * t62 + t98 * t89;
t39 = -t48 * mrSges(7,1) + t49 * mrSges(7,2);
t58 = qJD(6) - t61;
t40 = -t58 * mrSges(7,2) + t48 * mrSges(7,3);
t17 = m(7) * (t102 * t21 - t20 * t98) - t29 * mrSges(7,3) + t34 * mrSges(7,1) - t49 * t39 + t58 * t40;
t28 = -t49 * qJD(6) + t102 * t88 - t98 * t37;
t41 = t58 * mrSges(7,1) - t49 * mrSges(7,3);
t18 = m(7) * (t102 * t20 + t21 * t98) + t28 * mrSges(7,3) - t34 * mrSges(7,2) + t48 * t39 - t58 * t41;
t53 = -t89 * mrSges(6,2) + t61 * mrSges(6,3);
t54 = t89 * mrSges(6,1) - t62 * mrSges(6,3);
t113 = -m(6) * t107 + t36 * mrSges(6,1) - t37 * mrSges(6,2) - t102 * t17 - t98 * t18 + t61 * t53 - t62 * t54;
t66 = -t89 * mrSges(5,2) + t72 * mrSges(5,3);
t68 = t89 * mrSges(5,1) - t73 * mrSges(5,3);
t110 = -m(5) * t109 + t51 * mrSges(5,1) - t52 * mrSges(5,2) + t72 * t66 - t73 * t68 + t113;
t82 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t125;
t124 = qJD(2) * t104;
t83 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t124;
t132 = (t100 * t82 - t104 * t83) * qJD(2) + m(4) * (-t106 * pkin(8) + t112) - t79 * mrSges(4,1) + t78 * mrSges(4,2) - t110;
t12 = m(3) * t133 + qJDD(2) * mrSges(3,1) - t106 * mrSges(3,2) - t132;
t126 = t105 * t12;
t43 = -t61 * mrSges(6,1) + t62 * mrSges(6,2);
t13 = m(6) * t122 - t88 * mrSges(6,2) + t36 * mrSges(6,3) + t102 * t18 - t98 * t17 + t61 * t43 - t89 * t54;
t116 = -t95 * t23 + t92 * t25;
t111 = m(7) * (-t88 * pkin(5) - t87 * pkin(10) + (t131 + t44) * t62 + t116) - t28 * mrSges(7,1) + t29 * mrSges(7,2) - t48 * t40 + t49 * t41;
t14 = m(6) * (-0.2e1 * qJD(5) * t62 - t116) - t37 * mrSges(6,3) + t88 * mrSges(6,1) - t62 * t43 + t89 * t53 - t111;
t63 = -t72 * mrSges(5,1) + t73 * mrSges(5,2);
t10 = m(5) * t129 - t88 * mrSges(5,2) + t51 * mrSges(5,3) + t95 * t13 - t92 * t14 + t72 * t63 - t89 * t68;
t77 = (-mrSges(4,1) * t104 + mrSges(4,2) * t100) * qJD(2);
t9 = m(5) * t119 + t88 * mrSges(5,1) - t52 * mrSges(5,3) + t92 * t13 + t95 * t14 - t73 * t63 + t89 * t66;
t7 = m(4) * t117 + qJDD(3) * mrSges(4,1) - t78 * mrSges(4,3) + qJD(3) * t83 + t99 * t10 + t103 * t9 - t77 * t125;
t8 = m(4) * t128 - qJDD(3) * mrSges(4,2) + t79 * mrSges(4,3) - qJD(3) * t82 + t103 * t10 + t77 * t124 - t99 * t9;
t4 = m(3) * t121 - t106 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t100 * t7 + t104 * t8;
t6 = m(3) * t69 + t100 * t8 + t104 * t7;
t120 = m(2) * t91 + t94 * t126 + t4 * t127 + t97 * t6;
t2 = m(2) * t81 - t101 * t12 + t105 * t4;
t1 = m(2) * t80 - t94 * t6 + (t101 * t4 + t126) * t97;
t3 = [-m(1) * g(1) - t1 * t93 + t2 * t96, t2, t4, t8, t10, t13, t18; -m(1) * g(2) + t1 * t96 + t2 * t93, t1, t12, t7, t9, t14, t17; -m(1) * g(3) + t120, t120, t6, t132, -t110, -t113, t111;];
f_new  = t3;
