% Calculate vector of cutting forces with Newton-Euler
% S6RPRRRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 01:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRRP7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP7_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP7_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP7_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP7_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP7_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:39:49
% EndTime: 2019-05-06 01:39:56
% DurationCPUTime: 3.22s
% Computational Cost: add. (40964->191), mult. (96312->240), div. (0->0), fcn. (72300->10), ass. (0->98)
t105 = qJD(1) ^ 2;
t96 = cos(pkin(10));
t94 = t96 ^ 2;
t95 = sin(pkin(10));
t128 = t95 ^ 2 + t94;
t139 = t128 * mrSges(3,3);
t102 = cos(qJ(3));
t117 = -mrSges(3,1) * t96 + mrSges(3,2) * t95;
t114 = mrSges(3,3) * qJDD(1) + t105 * t117;
t104 = qJD(3) ^ 2;
t125 = qJD(1) * qJD(2);
t123 = -t96 * g(3) - 0.2e1 * t95 * t125;
t135 = pkin(2) * t105;
t100 = sin(qJ(1));
t103 = cos(qJ(1));
t116 = -g(1) * t103 - g(2) * t100;
t86 = -pkin(1) * t105 + qJDD(1) * qJ(2) + t116;
t59 = (-pkin(7) * qJDD(1) + t96 * t135 - t86) * t95 + t123;
t119 = -g(3) * t95 + (0.2e1 * t125 + t86) * t96;
t126 = qJDD(1) * t96;
t60 = pkin(7) * t126 - t94 * t135 + t119;
t99 = sin(qJ(3));
t120 = t102 * t59 - t99 * t60;
t134 = t95 * t99;
t84 = (t102 * t96 - t134) * qJD(1);
t115 = t102 * t95 + t96 * t99;
t85 = t115 * qJD(1);
t68 = -pkin(3) * t84 - pkin(8) * t85;
t30 = -qJDD(3) * pkin(3) - t104 * pkin(8) + t85 * t68 - t120;
t101 = cos(qJ(4));
t127 = t84 * qJD(3);
t71 = t115 * qJDD(1) + t127;
t98 = sin(qJ(4));
t75 = qJD(3) * t98 + t101 * t85;
t47 = -qJD(4) * t75 + qJDD(3) * t101 - t71 * t98;
t82 = qJD(4) - t84;
t58 = pkin(4) * t82 - pkin(9) * t75;
t74 = qJD(3) * t101 - t85 * t98;
t73 = t74 ^ 2;
t109 = -t47 * pkin(4) - t73 * pkin(9) + t75 * t58 + t30;
t136 = cos(qJ(5));
t48 = qJD(4) * t74 + qJDD(3) * t98 + t101 * t71;
t97 = sin(qJ(5));
t50 = t136 * t75 + t97 * t74;
t27 = t50 * qJD(5) - t136 * t47 + t97 * t48;
t49 = -t136 * t74 + t97 * t75;
t28 = -t49 * qJD(5) + t136 * t48 + t97 * t47;
t80 = qJD(5) + t82;
t41 = -t49 * mrSges(7,2) + mrSges(7,3) * t80;
t44 = -mrSges(7,1) * t80 + t50 * mrSges(7,2);
t112 = t28 * mrSges(7,3) + t50 * t44 - m(7) * (-0.2e1 * qJD(6) * t50 + (t49 * t80 - t28) * qJ(6) + (t50 * t80 + t27) * pkin(5) + t109) - t27 * mrSges(7,1) - t49 * t41;
t42 = -mrSges(6,2) * t80 - t49 * mrSges(6,3);
t43 = mrSges(6,1) * t80 - t50 * mrSges(6,3);
t108 = m(6) * t109 + t27 * mrSges(6,1) + t28 * mrSges(6,2) + t49 * t42 + t50 * t43 - t112;
t54 = -mrSges(5,2) * t82 + mrSges(5,3) * t74;
t55 = mrSges(5,1) * t82 - mrSges(5,3) * t75;
t106 = m(5) * t30 - t47 * mrSges(5,1) + t48 * mrSges(5,2) - t74 * t54 + t75 * t55 + t108;
t65 = -mrSges(4,1) * t84 + mrSges(4,2) * t85;
t76 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t84;
t12 = m(4) * t120 + qJDD(3) * mrSges(4,1) - t71 * mrSges(4,3) + qJD(3) * t76 - t85 * t65 - t106;
t129 = t102 * t60 + t99 * t59;
t31 = -pkin(3) * t104 + qJDD(3) * pkin(8) + t68 * t84 + t129;
t122 = t100 * g(1) - t103 * g(2);
t118 = qJDD(2) - t122;
t107 = (-pkin(2) * t96 - pkin(1)) * qJDD(1) + (-t128 * pkin(7) - qJ(2)) * t105 + t118;
t81 = t85 * qJD(3);
t70 = -qJDD(1) * t134 + t102 * t126 - t81;
t34 = (-t71 - t127) * pkin(8) + (-t70 + t81) * pkin(3) + t107;
t121 = t101 * t34 - t98 * t31;
t67 = qJDD(4) - t70;
t20 = (t74 * t82 - t48) * pkin(9) + (t74 * t75 + t67) * pkin(4) + t121;
t131 = t101 * t31 + t98 * t34;
t22 = -pkin(4) * t73 + t47 * pkin(9) - t58 * t82 + t131;
t132 = t136 * t22 + t97 * t20;
t37 = pkin(5) * t49 - qJ(6) * t50;
t63 = qJDD(5) + t67;
t78 = t80 ^ 2;
t124 = m(7) * (-pkin(5) * t78 + t63 * qJ(6) + 0.2e1 * qJD(6) * t80 - t49 * t37 + t132) + t80 * t44 + t63 * mrSges(7,3);
t38 = mrSges(7,1) * t49 - mrSges(7,3) * t50;
t130 = -mrSges(6,1) * t49 - mrSges(6,2) * t50 - t38;
t133 = -mrSges(6,3) - mrSges(7,2);
t13 = m(6) * t132 - t63 * mrSges(6,2) + t130 * t49 + t133 * t27 - t80 * t43 + t124;
t113 = t136 * t20 - t97 * t22;
t137 = m(7) * (-t63 * pkin(5) - t78 * qJ(6) + t50 * t37 + qJDD(6) - t113);
t14 = m(6) * t113 - t137 + (t42 + t41) * t80 + (mrSges(6,1) + mrSges(7,1)) * t63 + t130 * t50 + t133 * t28;
t51 = -mrSges(5,1) * t74 + mrSges(5,2) * t75;
t10 = m(5) * t121 + t67 * mrSges(5,1) - t48 * mrSges(5,3) + t97 * t13 + t136 * t14 - t75 * t51 + t82 * t54;
t11 = m(5) * t131 - t67 * mrSges(5,2) + t47 * mrSges(5,3) + t136 * t13 - t97 * t14 + t74 * t51 - t82 * t55;
t77 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t85;
t7 = m(4) * t129 - qJDD(3) * mrSges(4,2) + t70 * mrSges(4,3) - qJD(3) * t77 - t98 * t10 + t101 * t11 + t84 * t65;
t4 = m(3) * t123 + t99 * t7 + t102 * t12 + (-m(3) * t86 - t114) * t95;
t5 = m(3) * t119 + t102 * t7 + t114 * t96 - t99 * t12;
t138 = t96 * t4 + t95 * t5;
t111 = -m(4) * t107 + t70 * mrSges(4,1) - t71 * mrSges(4,2) - t101 * t10 - t98 * t11 + t84 * t76 - t85 * t77;
t110 = m(3) * (-qJDD(1) * pkin(1) - t105 * qJ(2) + t118) - t111;
t6 = m(2) * t122 + (-mrSges(2,2) + t139) * t105 + (mrSges(2,1) - t117) * qJDD(1) - t110;
t1 = m(2) * t116 - t105 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t95 * t4 + t96 * t5;
t2 = [-m(1) * g(1) + t1 * t103 - t100 * t6, t1, t5, t7, t11, t13, -t27 * mrSges(7,2) - t49 * t38 + t124; -m(1) * g(2) + t1 * t100 + t103 * t6, t6, t4, t12, t10, t14, -t112; (-m(1) - m(2)) * g(3) + t138, -m(2) * g(3) + t138, t117 * qJDD(1) - t105 * t139 + t110, -t111, t106, t108, -t63 * mrSges(7,1) + t28 * mrSges(7,2) + t50 * t38 - t80 * t41 + t137;];
f_new  = t2;
