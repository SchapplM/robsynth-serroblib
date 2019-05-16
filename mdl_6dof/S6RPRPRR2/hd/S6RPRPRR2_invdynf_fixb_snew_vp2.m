% Calculate vector of cutting forces with Newton-Euler
% S6RPRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-05-05 18:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRPRR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:22:23
% EndTime: 2019-05-05 18:22:31
% DurationCPUTime: 3.40s
% Computational Cost: add. (44996->179), mult. (97670->239), div. (0->0), fcn. (66372->12), ass. (0->94)
t102 = cos(qJ(3));
t105 = qJD(1) ^ 2;
t122 = qJD(1) * qJD(3);
t114 = t102 * t122;
t103 = cos(qJ(1));
t99 = sin(qJ(1));
t119 = t99 * g(1) - t103 * g(2);
t77 = qJDD(1) * pkin(1) + t119;
t113 = -t103 * g(1) - t99 * g(2);
t79 = -t105 * pkin(1) + t113;
t93 = sin(pkin(10));
t95 = cos(pkin(10));
t126 = t93 * t77 + t95 * t79;
t52 = -t105 * pkin(2) + qJDD(1) * pkin(7) + t126;
t91 = -g(3) + qJDD(2);
t98 = sin(qJ(3));
t116 = t102 * t91 - t98 * t52;
t80 = t98 * qJDD(1) + t114;
t35 = (-t80 + t114) * qJ(4) + (t102 * t105 * t98 + qJDD(3)) * pkin(3) + t116;
t127 = t102 * t52 + t98 * t91;
t81 = t102 * qJDD(1) - t98 * t122;
t125 = qJD(1) * t98;
t82 = qJD(3) * pkin(3) - qJ(4) * t125;
t90 = t102 ^ 2;
t38 = -t90 * t105 * pkin(3) + t81 * qJ(4) - qJD(3) * t82 + t127;
t92 = sin(pkin(11));
t94 = cos(pkin(11));
t71 = (t102 * t92 + t94 * t98) * qJD(1);
t130 = -0.2e1 * qJD(4) * t71 + t94 * t35 - t92 * t38;
t101 = cos(qJ(5));
t115 = t95 * t77 - t93 * t79;
t111 = -qJDD(1) * pkin(2) - t115;
t108 = -t81 * pkin(3) + qJDD(4) + t82 * t125 + (-qJ(4) * t90 - pkin(7)) * t105 + t111;
t100 = cos(qJ(6));
t104 = qJD(3) ^ 2;
t123 = qJD(1) * t102;
t70 = t94 * t123 - t92 * t125;
t120 = 0.2e1 * qJD(4) * t70 + t92 * t35 + t94 * t38;
t55 = -t70 * pkin(4) - t71 * pkin(8);
t23 = -t104 * pkin(4) + qJDD(3) * pkin(8) + t70 * t55 + t120;
t59 = -t92 * t80 + t94 * t81;
t60 = t94 * t80 + t92 * t81;
t29 = (-qJD(3) * t70 - t60) * pkin(8) + (qJD(3) * t71 - t59) * pkin(4) + t108;
t97 = sin(qJ(5));
t118 = t101 * t29 - t97 * t23;
t62 = t101 * qJD(3) - t97 * t71;
t42 = t62 * qJD(5) + t97 * qJDD(3) + t101 * t60;
t58 = qJDD(5) - t59;
t63 = t97 * qJD(3) + t101 * t71;
t69 = qJD(5) - t70;
t17 = (t62 * t69 - t42) * pkin(9) + (t62 * t63 + t58) * pkin(5) + t118;
t128 = t101 * t23 + t97 * t29;
t41 = -t63 * qJD(5) + t101 * qJDD(3) - t97 * t60;
t49 = t69 * pkin(5) - t63 * pkin(9);
t61 = t62 ^ 2;
t18 = -t61 * pkin(5) + t41 * pkin(9) - t69 * t49 + t128;
t96 = sin(qJ(6));
t43 = t100 * t62 - t96 * t63;
t26 = t43 * qJD(6) + t100 * t42 + t96 * t41;
t44 = t100 * t63 + t96 * t62;
t31 = -t43 * mrSges(7,1) + t44 * mrSges(7,2);
t66 = qJD(6) + t69;
t36 = -t66 * mrSges(7,2) + t43 * mrSges(7,3);
t56 = qJDD(6) + t58;
t15 = m(7) * (t100 * t17 - t96 * t18) - t26 * mrSges(7,3) + t56 * mrSges(7,1) - t44 * t31 + t66 * t36;
t25 = -t44 * qJD(6) + t100 * t41 - t96 * t42;
t37 = t66 * mrSges(7,1) - t44 * mrSges(7,3);
t16 = m(7) * (t100 * t18 + t96 * t17) + t25 * mrSges(7,3) - t56 * mrSges(7,2) + t43 * t31 - t66 * t37;
t45 = -t62 * mrSges(6,1) + t63 * mrSges(6,2);
t47 = -t69 * mrSges(6,2) + t62 * mrSges(6,3);
t12 = m(6) * t118 + t58 * mrSges(6,1) - t42 * mrSges(6,3) + t100 * t15 + t96 * t16 - t63 * t45 + t69 * t47;
t48 = t69 * mrSges(6,1) - t63 * mrSges(6,3);
t13 = m(6) * t128 - t58 * mrSges(6,2) + t41 * mrSges(6,3) + t100 * t16 - t96 * t15 + t62 * t45 - t69 * t48;
t64 = -qJD(3) * mrSges(5,2) + t70 * mrSges(5,3);
t65 = qJD(3) * mrSges(5,1) - t71 * mrSges(5,3);
t109 = -m(5) * t108 + t59 * mrSges(5,1) - t60 * mrSges(5,2) - t101 * t12 - t97 * t13 + t70 * t64 - t71 * t65;
t83 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t125;
t84 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t123;
t129 = -(t102 * t84 - t98 * t83) * qJD(1) + m(4) * (-t105 * pkin(7) + t111) - t81 * mrSges(4,1) + t80 * mrSges(4,2) - t109;
t22 = -qJDD(3) * pkin(4) - t104 * pkin(8) + t71 * t55 - t130;
t110 = t25 * mrSges(7,1) + t43 * t36 - m(7) * (-t41 * pkin(5) - t61 * pkin(9) + t63 * t49 + t22) - t26 * mrSges(7,2) - t44 * t37;
t106 = m(6) * t22 - t41 * mrSges(6,1) + t42 * mrSges(6,2) - t62 * t47 + t63 * t48 - t110;
t54 = -t70 * mrSges(5,1) + t71 * mrSges(5,2);
t14 = m(5) * t130 + qJDD(3) * mrSges(5,1) - t60 * mrSges(5,3) + qJD(3) * t64 - t71 * t54 - t106;
t78 = (-mrSges(4,1) * t102 + mrSges(4,2) * t98) * qJD(1);
t9 = m(5) * t120 - qJDD(3) * mrSges(5,2) + t59 * mrSges(5,3) - qJD(3) * t65 + t101 * t13 - t97 * t12 + t70 * t54;
t6 = m(4) * t116 + qJDD(3) * mrSges(4,1) - t80 * mrSges(4,3) + qJD(3) * t84 - t78 * t125 + t94 * t14 + t92 * t9;
t7 = m(4) * t127 - qJDD(3) * mrSges(4,2) + t81 * mrSges(4,3) - qJD(3) * t83 + t78 * t123 - t92 * t14 + t94 * t9;
t121 = m(3) * t91 + t102 * t6 + t98 * t7;
t8 = m(3) * t115 + qJDD(1) * mrSges(3,1) - t105 * mrSges(3,2) - t129;
t3 = m(3) * t126 - t105 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t102 * t7 - t98 * t6;
t2 = m(2) * t113 - t105 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t95 * t3 - t93 * t8;
t1 = m(2) * t119 + qJDD(1) * mrSges(2,1) - t105 * mrSges(2,2) + t93 * t3 + t95 * t8;
t4 = [-m(1) * g(1) - t99 * t1 + t103 * t2, t2, t3, t7, t9, t13, t16; -m(1) * g(2) + t103 * t1 + t99 * t2, t1, t8, t6, t14, t12, t15; (-m(1) - m(2)) * g(3) + t121, -m(2) * g(3) + t121, t121, t129, -t109, t106, -t110;];
f_new  = t4;
