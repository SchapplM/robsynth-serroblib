% Calculate vector of cutting forces with Newton-Euler
% S6RRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-05-06 13:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRPR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 13:01:17
% EndTime: 2019-05-06 13:01:27
% DurationCPUTime: 3.22s
% Computational Cost: add. (36977->204), mult. (86629->261), div. (0->0), fcn. (62470->10), ass. (0->100)
t102 = sin(qJ(2));
t105 = cos(qJ(2));
t107 = qJD(1) ^ 2;
t103 = sin(qJ(1));
t106 = cos(qJ(1));
t125 = t103 * g(1) - t106 * g(2);
t119 = -qJDD(1) * pkin(1) - t125;
t129 = qJD(1) * t102;
t127 = qJD(1) * qJD(2);
t89 = t105 * qJDD(1) - t102 * t127;
t90 = qJD(2) * pkin(2) - qJ(3) * t129;
t97 = t105 ^ 2;
t112 = -t89 * pkin(2) + qJDD(3) + t90 * t129 + (-qJ(3) * t97 - pkin(7)) * t107 + t119;
t88 = t102 * qJDD(1) + t105 * t127;
t98 = sin(pkin(10));
t99 = cos(pkin(10));
t72 = -t98 * t88 + t99 * t89;
t83 = (t102 * t99 + t105 * t98) * qJD(1);
t76 = qJD(2) * pkin(3) - t83 * pkin(8);
t82 = (-t102 * t98 + t105 * t99) * qJD(1);
t81 = t82 ^ 2;
t110 = -t72 * pkin(3) - t81 * pkin(8) + t83 * t76 + t112;
t100 = sin(qJ(6));
t104 = cos(qJ(6));
t101 = sin(qJ(4));
t142 = cos(qJ(4));
t65 = t101 * t83 - t142 * t82;
t96 = qJD(2) + qJD(4);
t140 = t65 * t96;
t145 = -2 * qJD(5);
t73 = t99 * t88 + t98 * t89;
t38 = -t65 * qJD(4) + t101 * t72 + t142 * t73;
t66 = t101 * t82 + t142 * t83;
t108 = (-t38 + t140) * qJ(5) + t110 + (t96 * pkin(4) + t145) * t66;
t121 = -t106 * g(1) - t103 * g(2);
t85 = -t107 * pkin(1) + qJDD(1) * pkin(7) + t121;
t130 = t102 * t85;
t141 = pkin(2) * t107;
t50 = qJDD(2) * pkin(2) - t88 * qJ(3) - t130 + (qJ(3) * t127 + t102 * t141 - g(3)) * t105;
t124 = -t102 * g(3) + t105 * t85;
t51 = t89 * qJ(3) - qJD(2) * t90 - t141 * t97 + t124;
t123 = -0.2e1 * qJD(3) * t83 + t99 * t50 - t98 * t51;
t23 = (qJD(2) * t82 - t73) * pkin(8) + (t82 * t83 + qJDD(2)) * pkin(3) + t123;
t126 = 0.2e1 * qJD(3) * t82 + t98 * t50 + t99 * t51;
t26 = -t81 * pkin(3) + t72 * pkin(8) - qJD(2) * t76 + t126;
t122 = -t101 * t26 + t142 * t23;
t44 = t65 * pkin(4) - t66 * qJ(5);
t94 = t96 ^ 2;
t95 = qJDD(2) + qJDD(4);
t19 = -t95 * pkin(4) - t94 * qJ(5) + t66 * t44 + qJDD(5) - t122;
t14 = (t65 * t66 - t95) * pkin(9) + (t38 + t140) * pkin(5) + t19;
t37 = t66 * qJD(4) + t101 * t73 - t142 * t72;
t60 = t66 * pkin(5) - t96 * pkin(9);
t64 = t65 ^ 2;
t17 = -t66 * t60 - t64 * pkin(5) + t108 + (pkin(4) + pkin(9)) * t37;
t54 = -t100 * t96 + t104 * t65;
t29 = t54 * qJD(6) + t100 * t37 + t104 * t95;
t36 = qJDD(6) + t38;
t55 = t100 * t65 + t104 * t96;
t39 = -t54 * mrSges(7,1) + t55 * mrSges(7,2);
t63 = qJD(6) + t66;
t40 = -t63 * mrSges(7,2) + t54 * mrSges(7,3);
t12 = m(7) * (-t100 * t17 + t104 * t14) - t29 * mrSges(7,3) + t36 * mrSges(7,1) - t55 * t39 + t63 * t40;
t28 = -t55 * qJD(6) - t100 * t95 + t104 * t37;
t41 = t63 * mrSges(7,1) - t55 * mrSges(7,3);
t13 = m(7) * (t100 * t14 + t104 * t17) + t28 * mrSges(7,3) - t36 * mrSges(7,2) + t54 * t39 - t63 * t41;
t59 = t66 * mrSges(6,1) + t96 * mrSges(6,2);
t118 = t100 * t12 - t104 * t13 - m(6) * (t37 * pkin(4) + t108) + t38 * mrSges(6,3) + t66 * t59;
t58 = t65 * mrSges(6,1) - t96 * mrSges(6,3);
t131 = t96 * mrSges(5,2) + t65 * mrSges(5,3) + t58;
t135 = mrSges(6,2) - mrSges(5,1);
t57 = t96 * mrSges(5,1) - t66 * mrSges(5,3);
t111 = m(5) * t110 + t38 * mrSges(5,2) - t131 * t65 - t135 * t37 + t66 * t57 - t118;
t74 = -qJD(2) * mrSges(4,2) + t82 * mrSges(4,3);
t75 = qJD(2) * mrSges(4,1) - t83 * mrSges(4,3);
t109 = m(4) * t112 - t72 * mrSges(4,1) + t73 * mrSges(4,2) - t82 * t74 + t83 * t75 + t111;
t91 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t129;
t128 = qJD(1) * t105;
t92 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t128;
t147 = t109 + (t102 * t91 - t105 * t92) * qJD(1) - t89 * mrSges(3,1) + t88 * mrSges(3,2) + m(3) * (-t107 * pkin(7) + t119);
t133 = t101 * t23 + t142 * t26;
t114 = -t94 * pkin(4) + t95 * qJ(5) - t65 * t44 + t133;
t116 = -t28 * mrSges(7,1) - t54 * t40 + m(7) * (-t37 * pkin(5) - t64 * pkin(9) + ((2 * qJD(5)) + t60) * t96 + t114) + t29 * mrSges(7,2) + t55 * t41;
t113 = -m(6) * (t96 * t145 - t114) + t116;
t46 = -t65 * mrSges(6,2) - t66 * mrSges(6,3);
t132 = -t65 * mrSges(5,1) - t66 * mrSges(5,2) - t46;
t134 = -mrSges(5,3) - mrSges(6,1);
t10 = m(5) * t133 + (-t57 + t59) * t96 + (-mrSges(5,2) + mrSges(6,3)) * t95 + t132 * t65 + t134 * t37 + t113;
t69 = -t82 * mrSges(4,1) + t83 * mrSges(4,2);
t117 = -m(6) * t19 - t100 * t13 - t104 * t12;
t9 = m(5) * t122 - t131 * t96 + t132 * t66 + t134 * t38 - t135 * t95 + t117;
t6 = m(4) * t123 + qJDD(2) * mrSges(4,1) - t73 * mrSges(4,3) + qJD(2) * t74 + t101 * t10 + t142 * t9 - t83 * t69;
t7 = m(4) * t126 - qJDD(2) * mrSges(4,2) + t72 * mrSges(4,3) - qJD(2) * t75 + t10 * t142 - t101 * t9 + t82 * t69;
t87 = (-mrSges(3,1) * t105 + mrSges(3,2) * t102) * qJD(1);
t4 = m(3) * (-t105 * g(3) - t130) - t88 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t87 * t129 + qJD(2) * t92 + t98 * t7 + t99 * t6;
t5 = m(3) * t124 - qJDD(2) * mrSges(3,2) + t89 * mrSges(3,3) - qJD(2) * t91 + t128 * t87 - t98 * t6 + t99 * t7;
t144 = t102 * t5 + t105 * t4;
t8 = m(2) * t125 + qJDD(1) * mrSges(2,1) - t107 * mrSges(2,2) - t147;
t1 = m(2) * t121 - t107 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t102 * t4 + t105 * t5;
t2 = [-m(1) * g(1) + t106 * t1 - t103 * t8, t1, t5, t7, t10, -t37 * mrSges(6,2) - t65 * t58 - t118, t13; -m(1) * g(2) + t103 * t1 + t106 * t8, t8, t4, t6, t9, t37 * mrSges(6,1) - t95 * mrSges(6,3) + t65 * t46 - t96 * t59 - t113, t12; (-m(1) - m(2)) * g(3) + t144, -m(2) * g(3) + t144, t147, t109, t111, t38 * mrSges(6,1) + t95 * mrSges(6,2) + t66 * t46 + t96 * t58 - t117, t116;];
f_new  = t2;
