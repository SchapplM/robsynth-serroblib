% Calculate vector of cutting forces with Newton-Euler
% S6RPRPPR1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
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
% Datum: 2019-05-05 16:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRPPR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:23:49
% EndTime: 2019-05-05 16:23:55
% DurationCPUTime: 3.31s
% Computational Cost: add. (42463->177), mult. (94917->240), div. (0->0), fcn. (63726->12), ass. (0->93)
t131 = -2 * qJD(4);
t126 = cos(pkin(10));
t101 = cos(qJ(3));
t104 = qJD(1) ^ 2;
t122 = qJD(1) * qJD(3);
t115 = t101 * t122;
t102 = cos(qJ(1));
t99 = sin(qJ(1));
t118 = t99 * g(1) - t102 * g(2);
t77 = qJDD(1) * pkin(1) + t118;
t113 = -g(1) * t102 - g(2) * t99;
t79 = -pkin(1) * t104 + t113;
t94 = sin(pkin(9));
t96 = cos(pkin(9));
t127 = t94 * t77 + t96 * t79;
t52 = -pkin(2) * t104 + qJDD(1) * pkin(7) + t127;
t91 = -g(3) + qJDD(2);
t98 = sin(qJ(3));
t117 = t101 * t91 - t98 * t52;
t80 = qJDD(1) * t98 + t115;
t35 = (-t80 + t115) * qJ(4) + (t101 * t104 * t98 + qJDD(3)) * pkin(3) + t117;
t128 = t101 * t52 + t98 * t91;
t81 = qJDD(1) * t101 - t98 * t122;
t125 = qJD(1) * t98;
t82 = qJD(3) * pkin(3) - qJ(4) * t125;
t90 = t101 ^ 2;
t38 = -pkin(3) * t104 * t90 + qJ(4) * t81 - qJD(3) * t82 + t128;
t93 = sin(pkin(10));
t71 = (t101 * t93 + t126 * t98) * qJD(1);
t130 = t126 * t35 + t71 * t131 - t93 * t38;
t116 = t96 * t77 - t94 * t79;
t110 = -qJDD(1) * pkin(2) - t116;
t107 = -t81 * pkin(3) + qJDD(4) + t82 * t125 + (-qJ(4) * t90 - pkin(7)) * t104 + t110;
t100 = cos(qJ(6));
t103 = qJD(3) ^ 2;
t123 = qJD(1) * t101;
t70 = -t126 * t123 + t93 * t125;
t119 = t126 * t38 + t70 * t131 + t93 * t35;
t54 = pkin(4) * t70 - qJ(5) * t71;
t23 = -pkin(4) * t103 + qJDD(3) * qJ(5) - t54 * t70 + t119;
t58 = -t126 * t81 + t80 * t93;
t59 = t126 * t80 + t93 * t81;
t26 = (qJD(3) * t70 - t59) * qJ(5) + (qJD(3) * t71 + t58) * pkin(4) + t107;
t92 = sin(pkin(11));
t95 = cos(pkin(11));
t64 = qJD(3) * t92 + t71 * t95;
t114 = -0.2e1 * qJD(5) * t64 - t92 * t23 + t95 * t26;
t50 = qJDD(3) * t92 + t59 * t95;
t63 = qJD(3) * t95 - t71 * t92;
t17 = (t63 * t70 - t50) * pkin(8) + (t63 * t64 + t58) * pkin(5) + t114;
t120 = 0.2e1 * qJD(5) * t63 + t95 * t23 + t92 * t26;
t47 = pkin(5) * t70 - t64 * pkin(8);
t49 = qJDD(3) * t95 - t59 * t92;
t62 = t63 ^ 2;
t18 = -t62 * pkin(5) + t49 * pkin(8) - t70 * t47 + t120;
t97 = sin(qJ(6));
t41 = t100 * t63 - t64 * t97;
t29 = t41 * qJD(6) + t100 * t50 + t49 * t97;
t42 = t100 * t64 + t63 * t97;
t31 = -mrSges(7,1) * t41 + mrSges(7,2) * t42;
t69 = qJD(6) + t70;
t36 = -mrSges(7,2) * t69 + t41 * mrSges(7,3);
t57 = qJDD(6) + t58;
t15 = m(7) * (t100 * t17 - t18 * t97) - t29 * mrSges(7,3) + t57 * mrSges(7,1) - t42 * t31 + t69 * t36;
t28 = -t42 * qJD(6) + t100 * t49 - t50 * t97;
t37 = mrSges(7,1) * t69 - t42 * mrSges(7,3);
t16 = m(7) * (t100 * t18 + t17 * t97) + t28 * mrSges(7,3) - t57 * mrSges(7,2) + t41 * t31 - t69 * t37;
t43 = -mrSges(6,1) * t63 + mrSges(6,2) * t64;
t45 = -mrSges(6,2) * t70 + t63 * mrSges(6,3);
t12 = m(6) * t114 + t58 * mrSges(6,1) - t50 * mrSges(6,3) + t100 * t15 + t97 * t16 - t64 * t43 + t70 * t45;
t46 = mrSges(6,1) * t70 - t64 * mrSges(6,3);
t13 = m(6) * t120 - t58 * mrSges(6,2) + t49 * mrSges(6,3) + t100 * t16 - t97 * t15 + t63 * t43 - t70 * t46;
t65 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t70;
t66 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t71;
t108 = m(5) * t107 + t58 * mrSges(5,1) + t59 * mrSges(5,2) + t95 * t12 + t92 * t13 + t70 * t65 + t71 * t66;
t83 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t125;
t84 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t123;
t129 = -(t101 * t84 - t83 * t98) * qJD(1) + m(4) * (-t104 * pkin(7) + t110) - t81 * mrSges(4,1) + t80 * mrSges(4,2) + t108;
t22 = -qJDD(3) * pkin(4) - t103 * qJ(5) + t71 * t54 + qJDD(5) - t130;
t109 = t28 * mrSges(7,1) + t41 * t36 - m(7) * (-t49 * pkin(5) - t62 * pkin(8) + t64 * t47 + t22) - t29 * mrSges(7,2) - t42 * t37;
t105 = m(6) * t22 - t49 * mrSges(6,1) + t50 * mrSges(6,2) - t63 * t45 + t64 * t46 - t109;
t55 = mrSges(5,1) * t70 + mrSges(5,2) * t71;
t14 = m(5) * t130 + qJDD(3) * mrSges(5,1) - t59 * mrSges(5,3) + qJD(3) * t65 - t71 * t55 - t105;
t78 = (-mrSges(4,1) * t101 + mrSges(4,2) * t98) * qJD(1);
t9 = m(5) * t119 - qJDD(3) * mrSges(5,2) - t58 * mrSges(5,3) - qJD(3) * t66 - t92 * t12 + t95 * t13 - t70 * t55;
t6 = m(4) * t117 + qJDD(3) * mrSges(4,1) - t80 * mrSges(4,3) + qJD(3) * t84 - t78 * t125 + t126 * t14 + t93 * t9;
t7 = m(4) * t128 - qJDD(3) * mrSges(4,2) + t81 * mrSges(4,3) - qJD(3) * t83 + t78 * t123 + t126 * t9 - t93 * t14;
t121 = m(3) * t91 + t101 * t6 + t98 * t7;
t8 = m(3) * t116 + qJDD(1) * mrSges(3,1) - t104 * mrSges(3,2) - t129;
t3 = m(3) * t127 - t104 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t101 * t7 - t98 * t6;
t2 = m(2) * t113 - t104 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t96 * t3 - t94 * t8;
t1 = m(2) * t118 + qJDD(1) * mrSges(2,1) - t104 * mrSges(2,2) + t94 * t3 + t96 * t8;
t4 = [-m(1) * g(1) - t1 * t99 + t102 * t2, t2, t3, t7, t9, t13, t16; -m(1) * g(2) + t1 * t102 + t2 * t99, t1, t8, t6, t14, t12, t15; (-m(1) - m(2)) * g(3) + t121, -m(2) * g(3) + t121, t121, t129, t108, t105, -t109;];
f_new  = t4;
