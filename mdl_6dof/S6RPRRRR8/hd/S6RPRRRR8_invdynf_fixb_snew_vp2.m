% Calculate vector of cutting forces with Newton-Euler
% S6RPRRRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-05-06 04:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRRR8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR8_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR8_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR8_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR8_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR8_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 04:14:19
% EndTime: 2019-05-06 04:14:24
% DurationCPUTime: 2.51s
% Computational Cost: add. (36662->180), mult. (71787->230), div. (0->0), fcn. (49181->10), ass. (0->94)
t102 = cos(qJ(1));
t97 = sin(qJ(1));
t114 = -t102 * g(1) - t97 * g(2);
t112 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t114;
t131 = -m(2) - m(3);
t130 = -pkin(1) - pkin(7);
t129 = (mrSges(2,1) - mrSges(3,2));
t128 = -mrSges(2,2) + mrSges(3,3);
t103 = qJD(1) ^ 2;
t101 = cos(qJ(3));
t123 = qJD(1) * t101;
t122 = qJD(1) * qJD(3);
t96 = sin(qJ(3));
t78 = -t96 * qJDD(1) - t101 * t122;
t82 = (qJD(3) * pkin(3)) - pkin(8) * t123;
t92 = t96 ^ 2;
t106 = -t78 * pkin(3) + t82 * t123 + (-pkin(8) * t92 + t130) * t103 + t112;
t100 = cos(qJ(4));
t95 = sin(qJ(4));
t71 = (t100 * t101 - t95 * t96) * qJD(1);
t118 = t96 * t122;
t79 = t101 * qJDD(1) - t118;
t47 = -t71 * qJD(4) + t100 * t78 - t95 * t79;
t124 = qJD(1) * t96;
t70 = -t100 * t124 - t95 * t123;
t48 = t70 * qJD(4) + t100 * t79 + t95 * t78;
t90 = qJD(3) + qJD(4);
t23 = (-t70 * t90 - t48) * pkin(9) + (t71 * t90 - t47) * pkin(4) + t106;
t119 = t97 * g(1) - t102 * g(2);
t110 = -t103 * qJ(2) + qJDD(2) - t119;
t64 = t130 * qJDD(1) + t110;
t125 = t96 * g(3) + t101 * t64;
t36 = (-t79 - t118) * pkin(8) + (-t101 * t103 * t96 + qJDD(3)) * pkin(3) + t125;
t117 = -t101 * g(3) + t96 * t64;
t39 = -t92 * t103 * pkin(3) + t78 * pkin(8) - qJD(3) * t82 + t117;
t126 = t100 * t39 + t95 * t36;
t54 = -t70 * pkin(4) - t71 * pkin(9);
t88 = t90 ^ 2;
t89 = qJDD(3) + qJDD(4);
t26 = -t88 * pkin(4) + t89 * pkin(9) + t70 * t54 + t126;
t94 = sin(qJ(5));
t99 = cos(qJ(5));
t127 = t94 * t23 + t99 * t26;
t69 = qJD(5) - t70;
t57 = t99 * t71 + t94 * t90;
t29 = -t57 * qJD(5) - t94 * t48 + t99 * t89;
t56 = -t94 * t71 + t99 * t90;
t30 = t56 * qJD(5) + t99 * t48 + t94 * t89;
t93 = sin(qJ(6));
t98 = cos(qJ(6));
t38 = t93 * t56 + t98 * t57;
t19 = -t38 * qJD(6) + t98 * t29 - t93 * t30;
t37 = t98 * t56 - t93 * t57;
t20 = t37 * qJD(6) + t93 * t29 + t98 * t30;
t115 = t100 * t36 - t95 * t39;
t25 = -t89 * pkin(4) - t88 * pkin(9) + t71 * t54 - t115;
t66 = qJD(6) + t69;
t31 = -t66 * mrSges(7,2) + t37 * mrSges(7,3);
t32 = t66 * mrSges(7,1) - t38 * mrSges(7,3);
t51 = t69 * pkin(5) - t57 * pkin(10);
t55 = t56 ^ 2;
t109 = t19 * mrSges(7,1) + t37 * t31 - m(7) * (-t29 * pkin(5) - t55 * pkin(10) + t57 * t51 + t25) - t20 * mrSges(7,2) - t38 * t32;
t49 = -t69 * mrSges(6,2) + t56 * mrSges(6,3);
t50 = t69 * mrSges(6,1) - t57 * mrSges(6,3);
t104 = m(6) * t25 - t29 * mrSges(6,1) + t30 * mrSges(6,2) - t56 * t49 + t57 * t50 - t109;
t53 = -t70 * mrSges(5,1) + t71 * mrSges(5,2);
t61 = -t90 * mrSges(5,2) + t70 * mrSges(5,3);
t11 = m(5) * t115 + t89 * mrSges(5,1) - t48 * mrSges(5,3) - t71 * t53 + t90 * t61 - t104;
t116 = t99 * t23 - t94 * t26;
t46 = qJDD(5) - t47;
t14 = (t56 * t69 - t30) * pkin(10) + (t56 * t57 + t46) * pkin(5) + t116;
t15 = -t55 * pkin(5) + t29 * pkin(10) - t69 * t51 + t127;
t28 = -t37 * mrSges(7,1) + t38 * mrSges(7,2);
t43 = qJDD(6) + t46;
t12 = m(7) * (t98 * t14 - t93 * t15) - t20 * mrSges(7,3) + t43 * mrSges(7,1) - t38 * t28 + t66 * t31;
t13 = m(7) * (t93 * t14 + t98 * t15) + t19 * mrSges(7,3) - t43 * mrSges(7,2) + t37 * t28 - t66 * t32;
t40 = -t56 * mrSges(6,1) + t57 * mrSges(6,2);
t10 = m(6) * t127 - t46 * mrSges(6,2) + t29 * mrSges(6,3) - t93 * t12 + t98 * t13 + t56 * t40 - t69 * t50;
t62 = t90 * mrSges(5,1) - t71 * mrSges(5,3);
t9 = m(6) * t116 + t46 * mrSges(6,1) - t30 * mrSges(6,3) + t98 * t12 + t93 * t13 - t57 * t40 + t69 * t49;
t6 = m(5) * t126 - t89 * mrSges(5,2) + t47 * mrSges(5,3) + t99 * t10 + t70 * t53 - t90 * t62 - t94 * t9;
t77 = (mrSges(4,1) * t96 + mrSges(4,2) * t101) * qJD(1);
t80 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t124;
t3 = m(4) * t125 + qJDD(3) * mrSges(4,1) - t79 * mrSges(4,3) + qJD(3) * t80 + t100 * t11 - t77 * t123 + t95 * t6;
t81 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t123;
t4 = m(4) * t117 - qJDD(3) * mrSges(4,2) + t78 * mrSges(4,3) - qJD(3) * t81 + t100 * t6 - t95 * t11 - t77 * t124;
t120 = t101 * t4 - t96 * t3;
t111 = -m(3) * (-qJDD(1) * pkin(1) + t110) - t101 * t3 - t96 * t4;
t108 = m(5) * t106 - t47 * mrSges(5,1) + t48 * mrSges(5,2) + t94 * t10 - t70 * t61 + t71 * t62 + t99 * t9;
t107 = -t78 * mrSges(4,1) + t108 + m(4) * (t130 * t103 + t112) + t80 * t124 + t81 * t123 + t79 * mrSges(4,2);
t105 = -m(3) * (t103 * pkin(1) - t112) + t107;
t5 = m(2) * t114 + t128 * qJDD(1) - (t129 * t103) + t105;
t1 = m(2) * t119 + t129 * qJDD(1) + t128 * t103 + t111;
t2 = [-m(1) * g(1) - t97 * t1 + t102 * t5, t5, -m(3) * g(3) + t120, t4, t6, t10, t13; -m(1) * g(2) + t102 * t1 + t97 * t5, t1, -(t103 * mrSges(3,2)) - qJDD(1) * mrSges(3,3) - t105, t3, t11, t9, t12; (-m(1) + t131) * g(3) + t120, t131 * g(3) + t120, qJDD(1) * mrSges(3,2) - t103 * mrSges(3,3) - t111, t107, t108, t104, -t109;];
f_new  = t2;
