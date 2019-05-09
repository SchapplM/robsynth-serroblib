% Calculate vector of cutting forces with Newton-Euler
% S6RPRRRP5
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
% Datum: 2019-05-06 01:32
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRRP5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:28:13
% EndTime: 2019-05-06 01:28:21
% DurationCPUTime: 3.51s
% Computational Cost: add. (43562->192), mult. (104548->239), div. (0->0), fcn. (81239->10), ass. (0->97)
t102 = qJD(1) ^ 2;
t94 = cos(pkin(10));
t91 = t94 ^ 2;
t93 = sin(pkin(10));
t126 = t93 ^ 2 + t91;
t139 = t126 * mrSges(3,3);
t100 = cos(qJ(3));
t123 = qJD(1) * qJD(2);
t119 = -t94 * g(3) - 0.2e1 * t93 * t123;
t125 = pkin(7) * qJDD(1);
t134 = pkin(2) * t102;
t101 = cos(qJ(1));
t98 = sin(qJ(1));
t114 = -g(1) * t101 - g(2) * t98;
t81 = -pkin(1) * t102 + qJDD(1) * qJ(2) + t114;
t59 = (t94 * t134 - t125 - t81) * t93 + t119;
t116 = -g(3) * t93 + (0.2e1 * t123 + t81) * t94;
t60 = t94 * t125 - t91 * t134 + t116;
t97 = sin(qJ(3));
t117 = t100 * t59 - t97 * t60;
t112 = t100 * t94 - t93 * t97;
t79 = t112 * qJD(1);
t124 = t79 * qJD(3);
t111 = t100 * t93 + t94 * t97;
t72 = t111 * qJDD(1) + t124;
t80 = t111 * qJD(1);
t25 = (-t72 + t124) * pkin(8) + (t79 * t80 + qJDD(3)) * pkin(3) + t117;
t127 = t100 * t60 + t97 * t59;
t71 = -t80 * qJD(3) + t112 * qJDD(1);
t75 = qJD(3) * pkin(3) - pkin(8) * t80;
t78 = t79 ^ 2;
t31 = -pkin(3) * t78 + pkin(8) * t71 - qJD(3) * t75 + t127;
t96 = sin(qJ(4));
t99 = cos(qJ(4));
t118 = t99 * t25 - t96 * t31;
t63 = t79 * t99 - t80 * t96;
t64 = t79 * t96 + t80 * t99;
t51 = -pkin(4) * t63 - pkin(9) * t64;
t92 = qJD(3) + qJD(4);
t88 = t92 ^ 2;
t89 = qJDD(3) + qJDD(4);
t20 = -t89 * pkin(4) - t88 * pkin(9) + t64 * t51 - t118;
t135 = cos(qJ(5));
t41 = t63 * qJD(4) + t71 * t96 + t72 * t99;
t95 = sin(qJ(5));
t53 = t135 * t64 + t95 * t92;
t27 = t53 * qJD(5) - t135 * t89 + t41 * t95;
t52 = -t135 * t92 + t64 * t95;
t28 = -t52 * qJD(5) + t135 * t41 + t95 * t89;
t62 = qJD(5) - t63;
t45 = -mrSges(7,2) * t52 + mrSges(7,3) * t62;
t121 = m(7) * (-0.2e1 * qJD(6) * t53 + (t52 * t62 - t28) * qJ(6) + (t53 * t62 + t27) * pkin(5) + t20) + t27 * mrSges(7,1) + t52 * t45;
t46 = -mrSges(6,2) * t62 - mrSges(6,3) * t52;
t47 = mrSges(6,1) * t62 - mrSges(6,3) * t53;
t48 = -mrSges(7,1) * t62 + mrSges(7,2) * t53;
t138 = m(6) * t20 + t27 * mrSges(6,1) + (t47 - t48) * t53 + (mrSges(6,2) - mrSges(7,3)) * t28 + t52 * t46 + t121;
t113 = -mrSges(3,1) * t94 + mrSges(3,2) * t93;
t110 = mrSges(3,3) * qJDD(1) + t102 * t113;
t50 = -mrSges(5,1) * t63 + mrSges(5,2) * t64;
t56 = -mrSges(5,2) * t92 + t63 * mrSges(5,3);
t10 = m(5) * t118 + t89 * mrSges(5,1) - t41 * mrSges(5,3) - t64 * t50 + t92 * t56 - t138;
t68 = -mrSges(4,1) * t79 + mrSges(4,2) * t80;
t73 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t79;
t130 = t96 * t25 + t99 * t31;
t21 = -pkin(4) * t88 + pkin(9) * t89 + t63 * t51 + t130;
t120 = t98 * g(1) - t101 * g(2);
t115 = qJDD(2) - t120;
t105 = (-pkin(2) * t94 - pkin(1)) * qJDD(1) + (-t126 * pkin(7) - qJ(2)) * t102 + t115;
t103 = -t71 * pkin(3) - t78 * pkin(8) + t80 * t75 + t105;
t40 = -t64 * qJD(4) + t71 * t99 - t72 * t96;
t23 = (-t63 * t92 - t41) * pkin(9) + (t64 * t92 - t40) * pkin(4) + t103;
t131 = t135 * t21 + t95 * t23;
t39 = qJDD(5) - t40;
t42 = pkin(5) * t52 - qJ(6) * t53;
t61 = t62 ^ 2;
t122 = m(7) * (-pkin(5) * t61 + qJ(6) * t39 + 0.2e1 * qJD(6) * t62 - t42 * t52 + t131) + t62 * t48 + t39 * mrSges(7,3);
t43 = mrSges(7,1) * t52 - mrSges(7,3) * t53;
t129 = -mrSges(6,1) * t52 - mrSges(6,2) * t53 - t43;
t132 = -mrSges(6,3) - mrSges(7,2);
t12 = m(6) * t131 - t39 * mrSges(6,2) + t129 * t52 + t132 * t27 - t62 * t47 + t122;
t109 = t135 * t23 - t95 * t21;
t136 = m(7) * (-t39 * pkin(5) - t61 * qJ(6) + t53 * t42 + qJDD(6) - t109);
t14 = m(6) * t109 - t136 + (t46 + t45) * t62 + t129 * t53 + (mrSges(6,1) + mrSges(7,1)) * t39 + t132 * t28;
t57 = mrSges(5,1) * t92 - t64 * mrSges(5,3);
t9 = m(5) * t130 - t89 * mrSges(5,2) + t40 * mrSges(5,3) + t135 * t12 - t95 * t14 + t63 * t50 - t92 * t57;
t6 = m(4) * t117 + qJDD(3) * mrSges(4,1) - t72 * mrSges(4,3) + qJD(3) * t73 + t99 * t10 - t80 * t68 + t96 * t9;
t74 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t80;
t7 = m(4) * t127 - qJDD(3) * mrSges(4,2) + t71 * mrSges(4,3) - qJD(3) * t74 - t96 * t10 + t79 * t68 + t99 * t9;
t4 = m(3) * t119 + t97 * t7 + t100 * t6 + (-m(3) * t81 - t110) * t93;
t5 = m(3) * t116 + t100 * t7 + t110 * t94 - t97 * t6;
t137 = t94 * t4 + t93 * t5;
t108 = -m(5) * t103 + t40 * mrSges(5,1) - t41 * mrSges(5,2) - t95 * t12 - t135 * t14 + t63 * t56 - t64 * t57;
t106 = -m(4) * t105 + t71 * mrSges(4,1) - t72 * mrSges(4,2) + t79 * t73 - t80 * t74 + t108;
t104 = m(3) * (-qJDD(1) * pkin(1) - t102 * qJ(2) + t115) - t106;
t8 = (-mrSges(2,2) + t139) * t102 + (mrSges(2,1) - t113) * qJDD(1) + m(2) * t120 - t104;
t1 = m(2) * t114 - t102 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t93 * t4 + t94 * t5;
t2 = [-m(1) * g(1) + t1 * t101 - t8 * t98, t1, t5, t7, t9, t12, -t27 * mrSges(7,2) - t52 * t43 + t122; -m(1) * g(2) + t1 * t98 + t101 * t8, t8, t4, t6, t10, t14, -t28 * mrSges(7,3) - t53 * t48 + t121; (-m(1) - m(2)) * g(3) + t137, -m(2) * g(3) + t137, t113 * qJDD(1) - t102 * t139 + t104, -t106, -t108, t138, -t39 * mrSges(7,1) + t28 * mrSges(7,2) + t53 * t43 - t62 * t45 + t136;];
f_new  = t2;
