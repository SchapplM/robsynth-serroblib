% Calculate vector of cutting forces with Newton-Euler
% S6PRRPRR6
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-05-05 05:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRPRR6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR6_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR6_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_invdynf_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR6_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR6_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR6_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 05:39:53
% EndTime: 2019-05-05 05:40:06
% DurationCPUTime: 7.19s
% Computational Cost: add. (130466->183), mult. (283280->260), div. (0->0), fcn. (229421->16), ass. (0->107)
t103 = cos(pkin(7));
t113 = qJD(2) ^ 2;
t108 = sin(qJ(2));
t112 = cos(qJ(2));
t100 = sin(pkin(6));
t131 = t100 * t112;
t104 = cos(pkin(6));
t102 = cos(pkin(12));
t98 = sin(pkin(12));
t87 = g(1) * t98 - g(2) * t102;
t134 = t104 * t87;
t88 = -g(1) * t102 - g(2) * t98;
t96 = -g(3) + qJDD(1);
t121 = -t108 * t88 + t112 * t134 + t96 * t131;
t99 = sin(pkin(7));
t59 = pkin(9) * t113 * t99 + qJDD(2) * pkin(2) + t121;
t74 = -t100 * t87 + t104 * t96;
t139 = t103 * t59 + t74 * t99;
t107 = sin(qJ(3));
t111 = cos(qJ(3));
t132 = t100 * t108;
t126 = t108 * t134 + t112 * t88 + t96 * t132;
t129 = qJDD(2) * t99;
t60 = -pkin(2) * t113 + pkin(9) * t129 + t126;
t138 = -t107 * t60 + t139 * t111;
t106 = sin(qJ(5));
t110 = cos(qJ(5));
t101 = cos(pkin(13));
t130 = qJD(2) * t111;
t124 = t99 * t130;
t127 = t139 * t107 + t111 * t60;
t133 = qJD(2) * t99;
t78 = (-pkin(3) * t111 - qJ(4) * t107) * t133;
t94 = qJD(2) * t103 + qJD(3);
t92 = t94 ^ 2;
t93 = qJDD(2) * t103 + qJDD(3);
t33 = -pkin(3) * t92 + qJ(4) * t93 + t78 * t124 + t127;
t70 = t103 * t74;
t80 = (qJD(3) * t130 + qJDD(2) * t107) * t99;
t125 = t107 * t133;
t81 = -qJD(3) * t125 + t111 * t129;
t36 = -pkin(3) * t81 - qJ(4) * t80 + t70 + (-t59 + (pkin(3) * t107 - qJ(4) * t111) * t94 * qJD(2)) * t99;
t97 = sin(pkin(13));
t73 = t101 * t125 + t94 * t97;
t122 = -0.2e1 * qJD(4) * t73 + t101 * t36 - t33 * t97;
t64 = t101 * t80 + t93 * t97;
t72 = t101 * t94 - t97 * t125;
t24 = (-t72 * t124 - t64) * pkin(10) + (t72 * t73 - t81) * pkin(4) + t122;
t128 = 0.2e1 * qJD(4) * t72 + t101 * t33 + t97 * t36;
t63 = t101 * t93 - t80 * t97;
t65 = -pkin(4) * t124 - pkin(10) * t73;
t71 = t72 ^ 2;
t26 = -pkin(4) * t71 + t63 * pkin(10) + t65 * t124 + t128;
t136 = t106 * t24 + t110 * t26;
t105 = sin(qJ(6));
t109 = cos(qJ(6));
t54 = -t106 * t73 + t110 * t72;
t55 = t106 * t72 + t110 * t73;
t46 = -pkin(5) * t54 - pkin(11) * t55;
t75 = qJDD(5) - t81;
t86 = qJD(5) - t124;
t85 = t86 ^ 2;
t21 = -pkin(5) * t85 + pkin(11) * t75 + t54 * t46 + t136;
t32 = -pkin(3) * t93 - qJ(4) * t92 + t78 * t125 + qJDD(4) - t138;
t115 = -t63 * pkin(4) - pkin(10) * t71 + t73 * t65 + t32;
t40 = -t55 * qJD(5) - t106 * t64 + t110 * t63;
t41 = t54 * qJD(5) + t106 * t63 + t110 * t64;
t22 = (-t54 * t86 - t41) * pkin(11) + (t55 * t86 - t40) * pkin(5) + t115;
t47 = -t105 * t55 + t109 * t86;
t30 = t47 * qJD(6) + t105 * t75 + t109 * t41;
t48 = t105 * t86 + t109 * t55;
t37 = -mrSges(7,1) * t47 + mrSges(7,2) * t48;
t39 = qJDD(6) - t40;
t53 = qJD(6) - t54;
t42 = -mrSges(7,2) * t53 + mrSges(7,3) * t47;
t18 = m(7) * (-t105 * t21 + t109 * t22) - t30 * mrSges(7,3) + t39 * mrSges(7,1) - t48 * t37 + t53 * t42;
t29 = -t48 * qJD(6) - t105 * t41 + t109 * t75;
t43 = mrSges(7,1) * t53 - mrSges(7,3) * t48;
t19 = m(7) * (t105 * t22 + t109 * t21) + t29 * mrSges(7,3) - t39 * mrSges(7,2) + t47 * t37 - t53 * t43;
t45 = -mrSges(6,1) * t54 + mrSges(6,2) * t55;
t50 = mrSges(6,1) * t86 - t55 * mrSges(6,3);
t14 = m(6) * t136 - t75 * mrSges(6,2) + t40 * mrSges(6,3) - t105 * t18 + t109 * t19 + t54 * t45 - t86 * t50;
t118 = -t106 * t26 + t110 * t24;
t116 = m(7) * (-pkin(5) * t75 - pkin(11) * t85 + t55 * t46 - t118) - t29 * mrSges(7,1) + t30 * mrSges(7,2) - t47 * t42 + t48 * t43;
t49 = -mrSges(6,2) * t86 + t54 * mrSges(6,3);
t15 = m(6) * t118 + t75 * mrSges(6,1) - t41 * mrSges(6,3) - t55 * t45 + t86 * t49 - t116;
t58 = -mrSges(5,1) * t72 + mrSges(5,2) * t73;
t61 = mrSges(5,2) * t124 + mrSges(5,3) * t72;
t11 = m(5) * t122 - t81 * mrSges(5,1) - t64 * mrSges(5,3) + t106 * t14 + t110 * t15 - t61 * t124 - t73 * t58;
t62 = -mrSges(5,1) * t124 - mrSges(5,3) * t73;
t12 = m(5) * t128 + t81 * mrSges(5,2) + t63 * mrSges(5,3) - t106 * t15 + t110 * t14 + t62 * t124 + t72 * t58;
t76 = mrSges(4,1) * t94 - mrSges(4,3) * t125;
t77 = -mrSges(4,2) * t94 + mrSges(4,3) * t124;
t10 = m(4) * (-t59 * t99 + t70) + t80 * mrSges(4,2) - t81 * mrSges(4,1) + t97 * t12 + t101 * t11 + (t107 * t76 - t111 * t77) * t133;
t117 = -m(6) * t115 + t40 * mrSges(6,1) - t41 * mrSges(6,2) - t105 * t19 - t109 * t18 + t54 * t49 - t55 * t50;
t114 = m(5) * t32 - t63 * mrSges(5,1) + t64 * mrSges(5,2) - t72 * t61 + t73 * t62 - t117;
t79 = (-mrSges(4,1) * t111 + mrSges(4,2) * t107) * t133;
t13 = m(4) * t138 + t93 * mrSges(4,1) - t80 * mrSges(4,3) - t79 * t125 + t94 * t77 - t114;
t9 = m(4) * t127 - t93 * mrSges(4,2) + t81 * mrSges(4,3) + t101 * t12 - t97 * t11 + t79 * t124 - t94 * t76;
t119 = t107 * t9 + t111 * t13;
t4 = m(3) * t121 + qJDD(2) * mrSges(3,1) - t113 * mrSges(3,2) - t99 * t10 + t119 * t103;
t6 = m(3) * t74 + t10 * t103 + t119 * t99;
t8 = m(3) * t126 - t113 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t107 * t13 + t111 * t9;
t123 = m(2) * t96 + t104 * t6 + t4 * t131 + t8 * t132;
t2 = m(2) * t88 - t108 * t4 + t112 * t8;
t1 = m(2) * t87 - t100 * t6 + (t108 * t8 + t112 * t4) * t104;
t3 = [-m(1) * g(1) - t1 * t98 + t102 * t2, t2, t8, t9, t12, t14, t19; -m(1) * g(2) + t1 * t102 + t2 * t98, t1, t4, t13, t11, t15, t18; -m(1) * g(3) + t123, t123, t6, t10, t114, -t117, t116;];
f_new  = t3;
