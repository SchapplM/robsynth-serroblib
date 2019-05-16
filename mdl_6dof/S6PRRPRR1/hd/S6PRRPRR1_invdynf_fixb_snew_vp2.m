% Calculate vector of cutting forces with Newton-Euler
% S6PRRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-05-05 04:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRPRR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:14:19
% EndTime: 2019-05-05 04:14:28
% DurationCPUTime: 5.11s
% Computational Cost: add. (70516->178), mult. (154524->245), div. (0->0), fcn. (114999->14), ass. (0->99)
t102 = sin(qJ(2));
t106 = cos(qJ(2));
t94 = sin(pkin(11));
t97 = cos(pkin(11));
t81 = t94 * g(1) - t97 * g(2);
t98 = cos(pkin(6));
t131 = t81 * t98;
t82 = -t97 * g(1) - t94 * g(2);
t92 = -g(3) + qJDD(1);
t95 = sin(pkin(6));
t133 = -t102 * t82 + t106 * (t92 * t95 + t131);
t101 = sin(qJ(3));
t105 = cos(qJ(3));
t107 = qJD(2) ^ 2;
t113 = -qJDD(2) * pkin(2) - t133;
t126 = qJD(2) * t101;
t124 = qJD(2) * qJD(3);
t80 = t105 * qJDD(2) - t101 * t124;
t83 = qJD(3) * pkin(3) - qJ(4) * t126;
t91 = t105 ^ 2;
t110 = -t80 * pkin(3) + qJDD(4) + t83 * t126 + (-qJ(4) * t91 - pkin(8)) * t107 + t113;
t103 = cos(qJ(6));
t120 = t105 * t124;
t79 = t101 * qJDD(2) + t120;
t93 = sin(pkin(12));
t96 = cos(pkin(12));
t61 = -t93 * t79 + t96 * t80;
t73 = (t101 * t96 + t105 * t93) * qJD(2);
t67 = qJD(3) * pkin(4) - t73 * pkin(9);
t72 = (-t101 * t93 + t105 * t96) * qJD(2);
t71 = t72 ^ 2;
t108 = -t61 * pkin(4) - t71 * pkin(9) + t73 * t67 + t110;
t100 = sin(qJ(5));
t104 = cos(qJ(5));
t128 = t102 * t95;
t122 = t102 * t131 + t106 * t82 + t92 * t128;
t53 = -t107 * pkin(2) + qJDD(2) * pkin(8) + t122;
t68 = -t95 * t81 + t98 * t92;
t119 = -t101 * t53 + t105 * t68;
t37 = (-t79 + t120) * qJ(4) + (t101 * t105 * t107 + qJDD(3)) * pkin(3) + t119;
t129 = t101 * t68 + t105 * t53;
t39 = -t91 * t107 * pkin(3) + t80 * qJ(4) - qJD(3) * t83 + t129;
t118 = -0.2e1 * qJD(4) * t73 + t96 * t37 - t93 * t39;
t62 = t96 * t79 + t93 * t80;
t23 = (qJD(3) * t72 - t62) * pkin(9) + (t72 * t73 + qJDD(3)) * pkin(4) + t118;
t123 = 0.2e1 * qJD(4) * t72 + t93 * t37 + t96 * t39;
t25 = -t71 * pkin(4) + t61 * pkin(9) - qJD(3) * t67 + t123;
t130 = t100 * t23 + t104 * t25;
t55 = -t100 * t73 + t104 * t72;
t56 = t100 * t72 + t104 * t73;
t44 = -t55 * pkin(5) - t56 * pkin(10);
t90 = qJD(3) + qJD(5);
t88 = t90 ^ 2;
t89 = qJDD(3) + qJDD(5);
t20 = -t88 * pkin(5) + t89 * pkin(10) + t55 * t44 + t130;
t32 = -t56 * qJD(5) - t100 * t62 + t104 * t61;
t33 = t55 * qJD(5) + t100 * t61 + t104 * t62;
t21 = (-t55 * t90 - t33) * pkin(10) + (t56 * t90 - t32) * pkin(5) + t108;
t99 = sin(qJ(6));
t47 = t103 * t90 - t99 * t56;
t27 = t47 * qJD(6) + t103 * t33 + t99 * t89;
t31 = qJDD(6) - t32;
t48 = t103 * t56 + t99 * t90;
t38 = -t47 * mrSges(7,1) + t48 * mrSges(7,2);
t54 = qJD(6) - t55;
t40 = -t54 * mrSges(7,2) + t47 * mrSges(7,3);
t17 = m(7) * (t103 * t21 - t99 * t20) - t27 * mrSges(7,3) + t31 * mrSges(7,1) - t48 * t38 + t54 * t40;
t26 = -t48 * qJD(6) + t103 * t89 - t99 * t33;
t41 = t54 * mrSges(7,1) - t48 * mrSges(7,3);
t18 = m(7) * (t103 * t20 + t99 * t21) + t26 * mrSges(7,3) - t31 * mrSges(7,2) + t47 * t38 - t54 * t41;
t50 = -t90 * mrSges(6,2) + t55 * mrSges(6,3);
t51 = t90 * mrSges(6,1) - t56 * mrSges(6,3);
t114 = -m(6) * t108 + t32 * mrSges(6,1) - t33 * mrSges(6,2) - t103 * t17 - t99 * t18 + t55 * t50 - t56 * t51;
t65 = -qJD(3) * mrSges(5,2) + t72 * mrSges(5,3);
t66 = qJD(3) * mrSges(5,1) - t73 * mrSges(5,3);
t111 = -m(5) * t110 + t61 * mrSges(5,1) - t62 * mrSges(5,2) + t72 * t65 - t73 * t66 + t114;
t84 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t126;
t125 = qJD(2) * t105;
t85 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t125;
t132 = (t101 * t84 - t105 * t85) * qJD(2) + m(4) * (-t107 * pkin(8) + t113) - t80 * mrSges(4,1) + t79 * mrSges(4,2) - t111;
t12 = m(3) * t133 + qJDD(2) * mrSges(3,1) - t107 * mrSges(3,2) - t132;
t127 = t106 * t12;
t43 = -t55 * mrSges(6,1) + t56 * mrSges(6,2);
t13 = m(6) * t130 - t89 * mrSges(6,2) + t32 * mrSges(6,3) + t103 * t18 - t99 * t17 + t55 * t43 - t90 * t51;
t116 = -t100 * t25 + t104 * t23;
t112 = m(7) * (-t89 * pkin(5) - t88 * pkin(10) + t56 * t44 - t116) - t26 * mrSges(7,1) + t27 * mrSges(7,2) - t47 * t40 + t48 * t41;
t14 = m(6) * t116 + t89 * mrSges(6,1) - t33 * mrSges(6,3) - t56 * t43 + t90 * t50 - t112;
t59 = -t72 * mrSges(5,1) + t73 * mrSges(5,2);
t10 = m(5) * t123 - qJDD(3) * mrSges(5,2) + t61 * mrSges(5,3) - qJD(3) * t66 - t100 * t14 + t104 * t13 + t72 * t59;
t78 = (-mrSges(4,1) * t105 + mrSges(4,2) * t101) * qJD(2);
t9 = m(5) * t118 + qJDD(3) * mrSges(5,1) - t62 * mrSges(5,3) + qJD(3) * t65 + t100 * t13 + t104 * t14 - t73 * t59;
t7 = m(4) * t119 + qJDD(3) * mrSges(4,1) - t79 * mrSges(4,3) + qJD(3) * t85 + t93 * t10 - t126 * t78 + t96 * t9;
t8 = m(4) * t129 - qJDD(3) * mrSges(4,2) + t80 * mrSges(4,3) - qJD(3) * t84 + t96 * t10 + t125 * t78 - t93 * t9;
t4 = m(3) * t122 - t107 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t101 * t7 + t105 * t8;
t6 = m(3) * t68 + t101 * t8 + t105 * t7;
t121 = m(2) * t92 + t95 * t127 + t4 * t128 + t98 * t6;
t2 = m(2) * t82 - t102 * t12 + t106 * t4;
t1 = m(2) * t81 - t95 * t6 + (t102 * t4 + t127) * t98;
t3 = [-m(1) * g(1) - t94 * t1 + t97 * t2, t2, t4, t8, t10, t13, t18; -m(1) * g(2) + t97 * t1 + t94 * t2, t1, t12, t7, t9, t14, t17; -m(1) * g(3) + t121, t121, t6, t132, -t111, -t114, t112;];
f_new  = t3;
