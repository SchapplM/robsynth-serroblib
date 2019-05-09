% Calculate vector of cutting forces with Newton-Euler
% S6RPRRRP4
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
% Datum: 2019-05-06 01:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRRP4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:23:34
% EndTime: 2019-05-06 01:23:42
% DurationCPUTime: 3.57s
% Computational Cost: add. (43980->192), mult. (105491->239), div. (0->0), fcn. (82111->10), ass. (0->96)
t104 = qJD(1) ^ 2;
t95 = cos(pkin(10));
t92 = t95 ^ 2;
t94 = sin(pkin(10));
t129 = t94 ^ 2 + t92;
t138 = t129 * mrSges(3,3);
t101 = cos(qJ(4));
t102 = cos(qJ(3));
t126 = qJD(1) * qJD(2);
t121 = -g(3) * t95 - 0.2e1 * t94 * t126;
t128 = pkin(7) * qJDD(1);
t135 = pkin(2) * t104;
t103 = cos(qJ(1));
t99 = sin(qJ(1));
t115 = -g(1) * t103 - g(2) * t99;
t84 = -pkin(1) * t104 + qJDD(1) * qJ(2) + t115;
t63 = (t95 * t135 - t128 - t84) * t94 + t121;
t117 = -g(3) * t94 + (0.2e1 * t126 + t84) * t95;
t64 = t95 * t128 - t92 * t135 + t117;
t98 = sin(qJ(3));
t119 = t102 * t63 - t64 * t98;
t113 = t102 * t95 - t94 * t98;
t82 = t113 * qJD(1);
t127 = qJD(3) * t82;
t112 = t102 * t94 + t95 * t98;
t75 = t112 * qJDD(1) + t127;
t83 = t112 * qJD(1);
t26 = (-t75 + t127) * pkin(8) + (t82 * t83 + qJDD(3)) * pkin(3) + t119;
t130 = t102 * t64 + t98 * t63;
t74 = -qJD(3) * t83 + t113 * qJDD(1);
t78 = qJD(3) * pkin(3) - pkin(8) * t83;
t81 = t82 ^ 2;
t33 = -pkin(3) * t81 + pkin(8) * t74 - qJD(3) * t78 + t130;
t97 = sin(qJ(4));
t118 = t101 * t26 - t97 * t33;
t66 = t101 * t82 - t83 * t97;
t67 = t101 * t83 + t82 * t97;
t54 = -t66 * pkin(4) - pkin(9) * t67;
t93 = qJD(3) + qJD(4);
t89 = t93 ^ 2;
t90 = qJDD(3) + qJDD(4);
t20 = -pkin(4) * t90 - pkin(9) * t89 + t67 * t54 - t118;
t100 = cos(qJ(5));
t44 = t66 * qJD(4) + t101 * t75 + t74 * t97;
t96 = sin(qJ(5));
t57 = t100 * t67 + t93 * t96;
t29 = -t57 * qJD(5) + t100 * t90 - t44 * t96;
t56 = t100 * t93 - t67 * t96;
t30 = t56 * qJD(5) + t100 * t44 + t90 * t96;
t65 = qJD(5) - t66;
t49 = pkin(5) * t65 - qJ(6) * t57;
t50 = mrSges(7,1) * t65 - mrSges(7,3) * t57;
t55 = t56 ^ 2;
t123 = m(7) * (-t29 * pkin(5) - t55 * qJ(6) + t57 * t49 + qJDD(6) + t20) + t30 * mrSges(7,2) + t57 * t50;
t47 = -mrSges(7,2) * t65 + mrSges(7,3) * t56;
t48 = -mrSges(6,2) * t65 + mrSges(6,3) * t56;
t51 = mrSges(6,1) * t65 - mrSges(6,3) * t57;
t137 = m(6) * t20 + t30 * mrSges(6,2) - (t48 + t47) * t56 - (mrSges(6,1) + mrSges(7,1)) * t29 + t57 * t51 + t123;
t114 = -t95 * mrSges(3,1) + t94 * mrSges(3,2);
t111 = mrSges(3,3) * qJDD(1) + t104 * t114;
t53 = -t66 * mrSges(5,1) + mrSges(5,2) * t67;
t60 = -mrSges(5,2) * t93 + t66 * mrSges(5,3);
t14 = m(5) * t118 + t90 * mrSges(5,1) - t44 * mrSges(5,3) - t67 * t53 + t93 * t60 - t137;
t71 = -mrSges(4,1) * t82 + mrSges(4,2) * t83;
t76 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t82;
t132 = t101 * t33 + t97 * t26;
t21 = -pkin(4) * t89 + pkin(9) * t90 + t66 * t54 + t132;
t122 = g(1) * t99 - t103 * g(2);
t116 = qJDD(2) - t122;
t107 = (-pkin(2) * t95 - pkin(1)) * qJDD(1) + (-t129 * pkin(7) - qJ(2)) * t104 + t116;
t105 = -pkin(3) * t74 - pkin(8) * t81 + t83 * t78 + t107;
t43 = -qJD(4) * t67 + t101 * t74 - t75 * t97;
t24 = (-t66 * t93 - t44) * pkin(9) + (t67 * t93 - t43) * pkin(4) + t105;
t120 = t100 * t24 - t21 * t96;
t42 = qJDD(5) - t43;
t125 = m(7) * (-0.2e1 * qJD(6) * t57 + (t56 * t65 - t30) * qJ(6) + (t56 * t57 + t42) * pkin(5) + t120) + t65 * t47 + t42 * mrSges(7,1);
t45 = -mrSges(7,1) * t56 + mrSges(7,2) * t57;
t46 = -mrSges(6,1) * t56 + mrSges(6,2) * t57;
t11 = m(6) * t120 + t42 * mrSges(6,1) + t65 * t48 + (-t46 - t45) * t57 + (-mrSges(6,3) - mrSges(7,3)) * t30 + t125;
t133 = t100 * t21 + t96 * t24;
t124 = m(7) * (-pkin(5) * t55 + qJ(6) * t29 + 0.2e1 * qJD(6) * t56 - t49 * t65 + t133) + t29 * mrSges(7,3) + t56 * t45;
t13 = m(6) * t133 + t29 * mrSges(6,3) + t56 * t46 + (-t51 - t50) * t65 + (-mrSges(6,2) - mrSges(7,2)) * t42 + t124;
t61 = mrSges(5,1) * t93 - mrSges(5,3) * t67;
t9 = m(5) * t132 - t90 * mrSges(5,2) + t43 * mrSges(5,3) + t100 * t13 - t96 * t11 + t66 * t53 - t93 * t61;
t6 = m(4) * t119 + qJDD(3) * mrSges(4,1) - t75 * mrSges(4,3) + qJD(3) * t76 + t101 * t14 - t83 * t71 + t97 * t9;
t77 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t83;
t7 = m(4) * t130 - qJDD(3) * mrSges(4,2) + t74 * mrSges(4,3) - qJD(3) * t77 + t101 * t9 - t97 * t14 + t82 * t71;
t4 = m(3) * t121 + t98 * t7 + t102 * t6 + (-m(3) * t84 - t111) * t94;
t5 = m(3) * t117 + t102 * t7 + t111 * t95 - t98 * t6;
t136 = t95 * t4 + t94 * t5;
t110 = -m(5) * t105 + t43 * mrSges(5,1) - t44 * mrSges(5,2) - t100 * t11 - t96 * t13 + t66 * t60 - t67 * t61;
t108 = -m(4) * t107 + t74 * mrSges(4,1) - t75 * mrSges(4,2) + t82 * t76 - t83 * t77 + t110;
t106 = m(3) * (-qJDD(1) * pkin(1) - qJ(2) * t104 + t116) - t108;
t8 = (-mrSges(2,2) + t138) * t104 + (mrSges(2,1) - t114) * qJDD(1) + m(2) * t122 - t106;
t1 = m(2) * t115 - t104 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t94 * t4 + t95 * t5;
t2 = [-m(1) * g(1) + t1 * t103 - t8 * t99, t1, t5, t7, t9, t13, -t42 * mrSges(7,2) - t65 * t50 + t124; -m(1) * g(2) + t1 * t99 + t103 * t8, t8, t4, t6, t14, t11, -t30 * mrSges(7,3) - t57 * t45 + t125; (-m(1) - m(2)) * g(3) + t136, -m(2) * g(3) + t136, t114 * qJDD(1) - t104 * t138 + t106, -t108, -t110, t137, -t29 * mrSges(7,1) - t56 * t47 + t123;];
f_new  = t2;
