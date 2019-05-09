% Calculate vector of cutting forces with Newton-Euler
% S6PRRRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-05-05 10:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRRRP5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 10:00:46
% EndTime: 2019-05-05 10:00:57
% DurationCPUTime: 5.00s
% Computational Cost: add. (69911->180), mult. (143594->244), div. (0->0), fcn. (113916->14), ass. (0->99)
t105 = qJD(2) ^ 2;
t100 = sin(qJ(2));
t104 = cos(qJ(2));
t93 = sin(pkin(6));
t124 = t104 * t93;
t91 = sin(pkin(12));
t94 = cos(pkin(12));
t82 = g(1) * t91 - g(2) * t94;
t96 = cos(pkin(6));
t130 = t82 * t96;
t83 = -g(1) * t94 - g(2) * t91;
t90 = -g(3) + qJDD(1);
t111 = -t100 * t83 + t104 * t130 + t90 * t124;
t92 = sin(pkin(7));
t133 = pkin(9) * t92;
t54 = qJDD(2) * pkin(2) + t105 * t133 + t111;
t69 = -t82 * t93 + t90 * t96;
t95 = cos(pkin(7));
t136 = t54 * t95 + t69 * t92;
t103 = cos(qJ(3));
t122 = qJD(2) * qJD(3);
t99 = sin(qJ(3));
t76 = (-qJDD(2) * t103 + t99 * t122) * t92;
t125 = t100 * t93;
t117 = t100 * t130 + t104 * t83 + t90 * t125;
t55 = -pkin(2) * t105 + qJDD(2) * t133 + t117;
t135 = t136 * t103 - t99 * t55;
t102 = cos(qJ(4));
t123 = qJD(2) * t92;
t115 = t103 * t123;
t118 = t103 * t55 + t136 * t99;
t74 = (-pkin(3) * t103 - pkin(10) * t99) * t123;
t88 = qJD(2) * t95 + qJD(3);
t86 = t88 ^ 2;
t87 = qJDD(2) * t95 + qJDD(3);
t29 = -pkin(3) * t86 + pkin(10) * t87 + t74 * t115 + t118;
t65 = t95 * t69;
t75 = (qJDD(2) * t99 + t103 * t122) * t92;
t31 = pkin(3) * t76 - pkin(10) * t75 + t65 + (-t54 + (pkin(3) * t99 - pkin(10) * t103) * t88 * qJD(2)) * t92;
t98 = sin(qJ(4));
t112 = t102 * t31 - t98 * t29;
t116 = t99 * t123;
t67 = t102 * t88 - t98 * t116;
t68 = t102 * t116 + t88 * t98;
t57 = -t67 * pkin(4) - pkin(11) * t68;
t70 = qJDD(4) + t76;
t81 = qJD(4) - t115;
t80 = t81 ^ 2;
t21 = -pkin(4) * t70 - pkin(11) * t80 + t68 * t57 - t112;
t101 = cos(qJ(5));
t50 = t67 * qJD(4) + t102 * t75 + t87 * t98;
t97 = sin(qJ(5));
t60 = t101 * t68 + t81 * t97;
t34 = -t60 * qJD(5) + t101 * t70 - t50 * t97;
t59 = t101 * t81 - t68 * t97;
t35 = t59 * qJD(5) + t101 * t50 + t70 * t97;
t66 = qJD(5) - t67;
t44 = pkin(5) * t66 - qJ(6) * t60;
t45 = mrSges(7,1) * t66 - mrSges(7,3) * t60;
t58 = t59 ^ 2;
t119 = m(7) * (-t34 * pkin(5) - t58 * qJ(6) + t60 * t44 + qJDD(6) + t21) + t35 * mrSges(7,2) + t60 * t45;
t42 = -mrSges(7,2) * t66 + mrSges(7,3) * t59;
t43 = -mrSges(6,2) * t66 + mrSges(6,3) * t59;
t46 = mrSges(6,1) * t66 - mrSges(6,3) * t60;
t134 = m(6) * t21 + t35 * mrSges(6,2) - (t43 + t42) * t59 - (mrSges(6,1) + mrSges(7,1)) * t34 + t60 * t46 + t119;
t127 = t102 * t29 + t98 * t31;
t22 = -pkin(4) * t80 + pkin(11) * t70 + t67 * t57 + t127;
t28 = -pkin(3) * t87 - pkin(10) * t86 + t74 * t116 - t135;
t49 = -qJD(4) * t68 + t102 * t87 - t75 * t98;
t25 = (-t67 * t81 - t50) * pkin(11) + (t68 * t81 - t49) * pkin(4) + t28;
t128 = t101 * t22 + t97 * t25;
t113 = t101 * t25 - t22 * t97;
t48 = qJDD(5) - t49;
t121 = m(7) * (-0.2e1 * qJD(6) * t60 + (t59 * t66 - t35) * qJ(6) + (t59 * t60 + t48) * pkin(5) + t113) + t66 * t42 + t48 * mrSges(7,1);
t39 = -mrSges(7,1) * t59 + mrSges(7,2) * t60;
t120 = m(7) * (-pkin(5) * t58 + qJ(6) * t34 + 0.2e1 * qJD(6) * t59 - t44 * t66 + t128) + t34 * mrSges(7,3) + t59 * t39;
t40 = -mrSges(6,1) * t59 + mrSges(6,2) * t60;
t13 = m(6) * t113 + t48 * mrSges(6,1) + t66 * t43 + (-t40 - t39) * t60 + (-mrSges(6,3) - mrSges(7,3)) * t35 + t121;
t14 = m(6) * t128 + t34 * mrSges(6,3) + t59 * t40 + (-t46 - t45) * t66 + (-mrSges(6,2) - mrSges(7,2)) * t48 + t120;
t56 = -t67 * mrSges(5,1) + mrSges(5,2) * t68;
t62 = mrSges(5,1) * t81 - mrSges(5,3) * t68;
t12 = m(5) * t127 - t70 * mrSges(5,2) + t49 * mrSges(5,3) + t101 * t14 - t97 * t13 + t67 * t56 - t81 * t62;
t61 = -mrSges(5,2) * t81 + t67 * mrSges(5,3);
t15 = m(5) * t112 + t70 * mrSges(5,1) - t50 * mrSges(5,3) - t68 * t56 + t81 * t61 - t134;
t71 = mrSges(4,1) * t88 - mrSges(4,3) * t116;
t72 = -mrSges(4,2) * t88 + mrSges(4,3) * t115;
t10 = m(4) * (-t54 * t92 + t65) + t75 * mrSges(4,2) + t76 * mrSges(4,1) + t98 * t12 + t102 * t15 + (-t103 * t72 + t71 * t99) * t123;
t106 = m(5) * t28 - t49 * mrSges(5,1) + t50 * mrSges(5,2) + t101 * t13 + t97 * t14 - t67 * t61 + t68 * t62;
t73 = (-mrSges(4,1) * t103 + mrSges(4,2) * t99) * t123;
t11 = m(4) * t135 + t87 * mrSges(4,1) - t75 * mrSges(4,3) - t73 * t116 + t88 * t72 - t106;
t9 = m(4) * t118 - t87 * mrSges(4,2) - t76 * mrSges(4,3) + t102 * t12 + t73 * t115 - t98 * t15 - t88 * t71;
t110 = t103 * t11 + t9 * t99;
t4 = m(3) * t111 + qJDD(2) * mrSges(3,1) - t105 * mrSges(3,2) - t92 * t10 + t110 * t95;
t6 = m(3) * t69 + t10 * t95 + t110 * t92;
t8 = m(3) * t117 - t105 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t103 * t9 - t99 * t11;
t114 = m(2) * t90 + t4 * t124 + t8 * t125 + t96 * t6;
t2 = m(2) * t83 - t100 * t4 + t104 * t8;
t1 = m(2) * t82 - t6 * t93 + (t100 * t8 + t104 * t4) * t96;
t3 = [-m(1) * g(1) - t1 * t91 + t2 * t94, t2, t8, t9, t12, t14, -t48 * mrSges(7,2) - t66 * t45 + t120; -m(1) * g(2) + t1 * t94 + t2 * t91, t1, t4, t11, t15, t13, -t35 * mrSges(7,3) - t60 * t39 + t121; -m(1) * g(3) + t114, t114, t6, t10, t106, t134, -t34 * mrSges(7,1) - t59 * t42 + t119;];
f_new  = t3;
