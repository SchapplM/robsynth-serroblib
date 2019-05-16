% Calculate vector of cutting forces with Newton-Euler
% S6RRPRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
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
% Datum: 2019-05-06 12:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRPP1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:16:02
% EndTime: 2019-05-06 12:16:11
% DurationCPUTime: 3.59s
% Computational Cost: add. (42979->202), mult. (99214->263), div. (0->0), fcn. (70090->10), ass. (0->97)
t141 = -2 * qJD(3);
t100 = sin(qJ(2));
t103 = cos(qJ(2));
t126 = qJD(1) * qJD(2);
t106 = qJD(1) ^ 2;
t101 = sin(qJ(1));
t104 = cos(qJ(1));
t118 = -t104 * g(1) - t101 * g(2);
t86 = -t106 * pkin(1) + qJDD(1) * pkin(7) + t118;
t131 = t100 * t86;
t135 = pkin(2) * t106;
t89 = t100 * qJDD(1) + t103 * t126;
t55 = qJDD(2) * pkin(2) - t89 * qJ(3) - t131 + (qJ(3) * t126 + t100 * t135 - g(3)) * t103;
t121 = -t100 * g(3) + t103 * t86;
t90 = t103 * qJDD(1) - t100 * t126;
t128 = qJD(1) * t100;
t91 = qJD(2) * pkin(2) - qJ(3) * t128;
t95 = t103 ^ 2;
t56 = t90 * qJ(3) - qJD(2) * t91 - t135 * t95 + t121;
t97 = sin(pkin(9));
t98 = cos(pkin(9));
t84 = (t100 * t98 + t103 * t97) * qJD(1);
t140 = t84 * t141 + t98 * t55 - t97 * t56;
t83 = (t100 * t97 - t103 * t98) * qJD(1);
t130 = cos(pkin(10));
t138 = -2 * qJD(5);
t102 = cos(qJ(4));
t105 = qJD(2) ^ 2;
t123 = t83 * t141 + t97 * t55 + t98 * t56;
t66 = t83 * pkin(3) - t84 * pkin(8);
t27 = -t105 * pkin(3) + qJDD(2) * pkin(8) - t83 * t66 + t123;
t122 = t101 * g(1) - t104 * g(2);
t115 = -qJDD(1) * pkin(1) - t122;
t110 = -t90 * pkin(2) + qJDD(3) + t91 * t128 + (-qJ(3) * t95 - pkin(7)) * t106 + t115;
t70 = -t97 * t89 + t98 * t90;
t71 = t98 * t89 + t97 * t90;
t30 = (qJD(2) * t83 - t71) * pkin(8) + (qJD(2) * t84 - t70) * pkin(3) + t110;
t99 = sin(qJ(4));
t120 = t102 * t30 - t99 * t27;
t74 = t102 * qJD(2) - t99 * t84;
t48 = t74 * qJD(4) + t99 * qJDD(2) + t102 * t71;
t69 = qJDD(4) - t70;
t75 = t99 * qJD(2) + t102 * t84;
t82 = qJD(4) + t83;
t20 = (t74 * t82 - t48) * qJ(5) + (t74 * t75 + t69) * pkin(4) + t120;
t133 = t102 * t27 + t99 * t30;
t47 = -t75 * qJD(4) + t102 * qJDD(2) - t99 * t71;
t62 = t82 * pkin(4) - t75 * qJ(5);
t73 = t74 ^ 2;
t22 = -t73 * pkin(4) + t47 * qJ(5) - t82 * t62 + t133;
t96 = sin(pkin(10));
t53 = -t130 * t74 + t96 * t75;
t124 = t130 * t22 + t138 * t53 + t96 * t20;
t54 = t130 * t75 + t96 * t74;
t37 = t53 * pkin(5) - t54 * qJ(6);
t44 = -t82 * mrSges(7,1) + t54 * mrSges(7,2);
t81 = t82 ^ 2;
t125 = m(7) * (-t81 * pkin(5) + t69 * qJ(6) + 0.2e1 * qJD(6) * t82 - t53 * t37 + t124) + t82 * t44 + t69 * mrSges(7,3);
t38 = t53 * mrSges(7,1) - t54 * mrSges(7,3);
t132 = -t53 * mrSges(6,1) - t54 * mrSges(6,2) - t38;
t134 = -mrSges(6,3) - mrSges(7,2);
t33 = -t130 * t47 + t96 * t48;
t43 = t82 * mrSges(6,1) - t54 * mrSges(6,3);
t12 = m(6) * t124 - t69 * mrSges(6,2) + t132 * t53 + t134 * t33 - t82 * t43 + t125;
t114 = t130 * t20 - t96 * t22;
t136 = m(7) * (-t69 * pkin(5) - t81 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t37) * t54 - t114);
t34 = t130 * t48 + t96 * t47;
t41 = -t53 * mrSges(7,2) + t82 * mrSges(7,3);
t42 = -t82 * mrSges(6,2) - t53 * mrSges(6,3);
t13 = m(6) * t114 - t136 + (t42 + t41) * t82 + (mrSges(6,1) + mrSges(7,1)) * t69 + (m(6) * t138 + t132) * t54 + t134 * t34;
t57 = -t74 * mrSges(5,1) + t75 * mrSges(5,2);
t61 = -t82 * mrSges(5,2) + t74 * mrSges(5,3);
t10 = m(5) * t120 + t69 * mrSges(5,1) - t48 * mrSges(5,3) + t96 * t12 + t13 * t130 - t75 * t57 + t82 * t61;
t63 = t82 * mrSges(5,1) - t75 * mrSges(5,3);
t11 = m(5) * t133 - t69 * mrSges(5,2) + t47 * mrSges(5,3) + t12 * t130 - t96 * t13 + t74 * t57 - t82 * t63;
t76 = -qJD(2) * mrSges(4,2) - t83 * mrSges(4,3);
t77 = qJD(2) * mrSges(4,1) - t84 * mrSges(4,3);
t112 = -m(4) * t110 + t70 * mrSges(4,1) - t71 * mrSges(4,2) - t102 * t10 - t99 * t11 - t83 * t76 - t84 * t77;
t92 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t128;
t127 = qJD(1) * t103;
t93 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t127;
t139 = (t100 * t92 - t103 * t93) * qJD(1) + m(3) * (-t106 * pkin(7) + t115) - t90 * mrSges(3,1) + t89 * mrSges(3,2) - t112;
t26 = -qJDD(2) * pkin(3) - t105 * pkin(8) + t84 * t66 - t140;
t108 = -t47 * pkin(4) - t73 * qJ(5) + t75 * t62 + qJDD(5) + t26;
t113 = t34 * mrSges(7,3) + t54 * t44 - m(7) * (-0.2e1 * qJD(6) * t54 + (t53 * t82 - t34) * qJ(6) + (t54 * t82 + t33) * pkin(5) + t108) - t33 * mrSges(7,1) - t53 * t41;
t111 = m(6) * t108 + t33 * mrSges(6,1) + t34 * mrSges(6,2) + t53 * t42 + t54 * t43 - t113;
t107 = m(5) * t26 - t47 * mrSges(5,1) + t48 * mrSges(5,2) - t74 * t61 + t75 * t63 + t111;
t65 = t83 * mrSges(4,1) + t84 * mrSges(4,2);
t14 = m(4) * t140 + qJDD(2) * mrSges(4,1) - t71 * mrSges(4,3) + qJD(2) * t76 - t84 * t65 - t107;
t7 = m(4) * t123 - qJDD(2) * mrSges(4,2) + t70 * mrSges(4,3) - qJD(2) * t77 - t99 * t10 + t102 * t11 - t83 * t65;
t88 = (-mrSges(3,1) * t103 + mrSges(3,2) * t100) * qJD(1);
t4 = m(3) * (-t103 * g(3) - t131) - t89 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t88 * t128 + qJD(2) * t93 + t97 * t7 + t98 * t14;
t5 = m(3) * t121 - qJDD(2) * mrSges(3,2) + t90 * mrSges(3,3) - qJD(2) * t92 + t127 * t88 - t97 * t14 + t98 * t7;
t137 = t100 * t5 + t103 * t4;
t6 = m(2) * t122 + qJDD(1) * mrSges(2,1) - t106 * mrSges(2,2) - t139;
t1 = m(2) * t118 - t106 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t100 * t4 + t103 * t5;
t2 = [-m(1) * g(1) + t104 * t1 - t101 * t6, t1, t5, t7, t11, t12, -t33 * mrSges(7,2) - t53 * t38 + t125; -m(1) * g(2) + t101 * t1 + t104 * t6, t6, t4, t14, t10, t13, -t113; (-m(1) - m(2)) * g(3) + t137, -m(2) * g(3) + t137, t139, -t112, t107, t111, -t69 * mrSges(7,1) + t34 * mrSges(7,2) + t54 * t38 - t82 * t41 + t136;];
f_new  = t2;
