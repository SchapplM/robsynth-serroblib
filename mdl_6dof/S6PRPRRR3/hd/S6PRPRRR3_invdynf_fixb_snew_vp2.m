% Calculate vector of cutting forces with Newton-Euler
% S6PRPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-05-05 00:44
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRPRRR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 00:39:57
% EndTime: 2019-05-05 00:40:07
% DurationCPUTime: 4.75s
% Computational Cost: add. (65284->167), mult. (147159->223), div. (0->0), fcn. (115525->14), ass. (0->99)
t102 = qJD(2) ^ 2;
t91 = cos(pkin(12));
t85 = t91 ^ 2;
t88 = sin(pkin(12));
t124 = t88 ^ 2 + t85;
t132 = mrSges(4,3) * t124;
t101 = cos(qJ(2));
t89 = sin(pkin(11));
t92 = cos(pkin(11));
t76 = t89 * g(1) - t92 * g(2);
t93 = cos(pkin(6));
t129 = t76 * t93;
t77 = -t92 * g(1) - t89 * g(2);
t87 = -g(3) + qJDD(1);
t90 = sin(pkin(6));
t97 = sin(qJ(2));
t131 = t101 * (t87 * t90 + t129) - t97 * t77;
t130 = pkin(3) * t102;
t128 = t90 * t97;
t100 = cos(qJ(4));
t122 = pkin(8) * qJDD(2);
t120 = qJD(2) * qJD(3);
t68 = -t90 * t76 + t93 * t87;
t125 = -0.2e1 * t88 * t120 + t91 * t68;
t118 = t101 * t77 + t87 * t128 + t97 * t129;
t53 = -t102 * pkin(2) + qJDD(2) * qJ(3) + t118;
t38 = (t130 * t91 - t122 - t53) * t88 + t125;
t119 = t88 * t68 + (0.2e1 * t120 + t53) * t91;
t41 = t122 * t91 - t130 * t85 + t119;
t96 = sin(qJ(4));
t116 = t100 * t38 - t96 * t41;
t112 = t100 * t91 - t88 * t96;
t70 = t112 * qJD(2);
t121 = t70 * qJD(4);
t111 = t100 * t88 + t91 * t96;
t62 = qJDD(2) * t111 + t121;
t71 = t111 * qJD(2);
t23 = (-t62 + t121) * pkin(9) + (t70 * t71 + qJDD(4)) * pkin(4) + t116;
t126 = t100 * t41 + t96 * t38;
t61 = -t71 * qJD(4) + qJDD(2) * t112;
t67 = qJD(4) * pkin(4) - t71 * pkin(9);
t69 = t70 ^ 2;
t25 = -t69 * pkin(4) + t61 * pkin(9) - qJD(4) * t67 + t126;
t95 = sin(qJ(5));
t99 = cos(qJ(5));
t127 = t95 * t23 + t99 * t25;
t109 = qJDD(3) - t131;
t104 = (-pkin(3) * t91 - pkin(2)) * qJDD(2) + (-pkin(8) * t124 - qJ(3)) * t102 + t109;
t103 = -t61 * pkin(4) - t69 * pkin(9) + t71 * t67 + t104;
t55 = t99 * t70 - t95 * t71;
t56 = t95 * t70 + t99 * t71;
t44 = -t55 * pkin(5) - t56 * pkin(10);
t86 = qJD(4) + qJD(5);
t82 = t86 ^ 2;
t83 = qJDD(4) + qJDD(5);
t20 = -t82 * pkin(5) + t83 * pkin(10) + t55 * t44 + t127;
t32 = -t56 * qJD(5) + t99 * t61 - t95 * t62;
t33 = t55 * qJD(5) + t95 * t61 + t99 * t62;
t21 = (-t55 * t86 - t33) * pkin(10) + (t56 * t86 - t32) * pkin(5) + t103;
t94 = sin(qJ(6));
t98 = cos(qJ(6));
t47 = -t94 * t56 + t98 * t86;
t27 = t47 * qJD(6) + t98 * t33 + t94 * t83;
t31 = qJDD(6) - t32;
t48 = t98 * t56 + t94 * t86;
t34 = -t47 * mrSges(7,1) + t48 * mrSges(7,2);
t54 = qJD(6) - t55;
t39 = -t54 * mrSges(7,2) + t47 * mrSges(7,3);
t17 = m(7) * (-t94 * t20 + t98 * t21) - t27 * mrSges(7,3) + t31 * mrSges(7,1) - t48 * t34 + t54 * t39;
t26 = -t48 * qJD(6) - t94 * t33 + t98 * t83;
t40 = t54 * mrSges(7,1) - t48 * mrSges(7,3);
t18 = m(7) * (t98 * t20 + t94 * t21) + t26 * mrSges(7,3) - t31 * mrSges(7,2) + t47 * t34 - t54 * t40;
t51 = -t86 * mrSges(6,2) + t55 * mrSges(6,3);
t52 = t86 * mrSges(6,1) - t56 * mrSges(6,3);
t108 = -m(6) * t103 + t32 * mrSges(6,1) - t33 * mrSges(6,2) - t98 * t17 - t94 * t18 + t55 * t51 - t56 * t52;
t65 = -qJD(4) * mrSges(5,2) + t70 * mrSges(5,3);
t66 = qJD(4) * mrSges(5,1) - t71 * mrSges(5,3);
t106 = -m(5) * t104 + t61 * mrSges(5,1) - t62 * mrSges(5,2) + t70 * t65 - t71 * t66 + t108;
t105 = m(4) * (-qJDD(2) * pkin(2) - t102 * qJ(3) + t109) - t106;
t115 = -t91 * mrSges(4,1) + t88 * mrSges(4,2);
t12 = -t105 + m(3) * t131 + (-mrSges(3,2) + t132) * t102 + (mrSges(3,1) - t115) * qJDD(2);
t123 = t101 * t12;
t43 = -t55 * mrSges(6,1) + t56 * mrSges(6,2);
t13 = m(6) * t127 - t83 * mrSges(6,2) + t32 * mrSges(6,3) - t94 * t17 + t98 * t18 + t55 * t43 - t86 * t52;
t114 = t99 * t23 - t95 * t25;
t107 = m(7) * (-t83 * pkin(5) - t82 * pkin(10) + t56 * t44 - t114) - t26 * mrSges(7,1) + t27 * mrSges(7,2) - t47 * t39 + t48 * t40;
t14 = m(6) * t114 + t83 * mrSges(6,1) - t33 * mrSges(6,3) - t56 * t43 + t86 * t51 - t107;
t59 = -t70 * mrSges(5,1) + t71 * mrSges(5,2);
t10 = m(5) * t126 - qJDD(4) * mrSges(5,2) + t61 * mrSges(5,3) - qJD(4) * t66 + t99 * t13 - t95 * t14 + t70 * t59;
t110 = qJDD(2) * mrSges(4,3) + t102 * t115;
t9 = m(5) * t116 + qJDD(4) * mrSges(5,1) - t62 * mrSges(5,3) + qJD(4) * t65 + t95 * t13 + t99 * t14 - t71 * t59;
t7 = m(4) * t125 + t96 * t10 + t100 * t9 + (-m(4) * t53 - t110) * t88;
t8 = m(4) * t119 + t100 * t10 + t110 * t91 - t96 * t9;
t4 = m(3) * t118 - t102 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t88 * t7 + t91 * t8;
t6 = m(3) * t68 + t91 * t7 + t88 * t8;
t117 = m(2) * t87 + t90 * t123 + t4 * t128 + t93 * t6;
t2 = m(2) * t77 + t101 * t4 - t97 * t12;
t1 = m(2) * t76 - t90 * t6 + (t4 * t97 + t123) * t93;
t3 = [-m(1) * g(1) - t89 * t1 + t92 * t2, t2, t4, t8, t10, t13, t18; -m(1) * g(2) + t92 * t1 + t89 * t2, t1, t12, t7, t9, t14, t17; -m(1) * g(3) + t117, t117, t6, qJDD(2) * t115 - t102 * t132 + t105, -t106, -t108, t107;];
f_new  = t3;
