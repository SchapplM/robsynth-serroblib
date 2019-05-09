% Calculate vector of cutting forces with Newton-Euler
% S6RRRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-05-07 07:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPRP5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:51:01
% EndTime: 2019-05-07 07:51:11
% DurationCPUTime: 4.44s
% Computational Cost: add. (60266->200), mult. (125677->256), div. (0->0), fcn. (88874->10), ass. (0->98)
t101 = sin(qJ(3));
t102 = sin(qJ(2));
t104 = cos(qJ(3));
t105 = cos(qJ(2));
t126 = qJD(1) * qJD(2);
t121 = t105 * t126;
t108 = qJD(1) ^ 2;
t103 = sin(qJ(1));
t106 = cos(qJ(1));
t123 = t103 * g(1) - t106 * g(2);
t79 = -qJDD(1) * pkin(1) - t108 * pkin(7) - t123;
t88 = t102 * qJDD(1) + t121;
t95 = t102 * t126;
t89 = t105 * qJDD(1) - t95;
t55 = (-t88 - t121) * pkin(8) + (-t89 + t95) * pkin(2) + t79;
t107 = qJD(2) ^ 2;
t118 = -t106 * g(1) - t103 * g(2);
t80 = -t108 * pkin(1) + qJDD(1) * pkin(7) + t118;
t122 = -t102 * g(3) + t105 * t80;
t127 = t105 * qJD(1);
t87 = (-pkin(2) * t105 - pkin(8) * t102) * qJD(1);
t58 = -t107 * pkin(2) + qJDD(2) * pkin(8) + t87 * t127 + t122;
t120 = -t101 * t58 + t104 * t55;
t128 = qJD(1) * t102;
t84 = t104 * qJD(2) - t101 * t128;
t67 = t84 * qJD(3) + t101 * qJDD(2) + t104 * t88;
t100 = sin(qJ(5));
t134 = cos(qJ(5));
t83 = qJDD(3) - t89;
t85 = t101 * qJD(2) + t104 * t128;
t94 = qJD(3) - t127;
t29 = (t84 * t94 - t67) * qJ(4) + (t84 * t85 + t83) * pkin(3) + t120;
t130 = t101 * t55 + t104 * t58;
t66 = -t85 * qJD(3) + t104 * qJDD(2) - t101 * t88;
t73 = t94 * pkin(3) - t85 * qJ(4);
t82 = t84 ^ 2;
t31 = -t82 * pkin(3) + t66 * qJ(4) - t94 * t73 + t130;
t98 = sin(pkin(10));
t99 = cos(pkin(10));
t70 = t98 * t84 + t99 * t85;
t119 = -0.2e1 * qJD(4) * t70 + t99 * t29 - t98 * t31;
t47 = t98 * t66 + t99 * t67;
t69 = t99 * t84 - t98 * t85;
t18 = (t69 * t94 - t47) * pkin(9) + (t69 * t70 + t83) * pkin(4) + t119;
t124 = 0.2e1 * qJD(4) * t69 + t98 * t29 + t99 * t31;
t46 = t99 * t66 - t98 * t67;
t61 = t94 * pkin(4) - t70 * pkin(9);
t68 = t69 ^ 2;
t20 = -t68 * pkin(4) + t46 * pkin(9) - t94 * t61 + t124;
t116 = -t100 * t20 + t134 * t18;
t50 = t100 * t70 - t134 * t69;
t51 = t100 * t69 + t134 * t70;
t35 = t50 * mrSges(7,1) - t51 * mrSges(7,3);
t131 = -t50 * mrSges(6,1) - t51 * mrSges(6,2) - t35;
t133 = -mrSges(6,3) - mrSges(7,2);
t34 = t50 * pkin(5) - t51 * qJ(6);
t81 = qJDD(5) + t83;
t93 = qJD(5) + t94;
t92 = t93 ^ 2;
t136 = m(7) * (-t81 * pkin(5) - t92 * qJ(6) + t51 * t34 + qJDD(6) - t116);
t26 = -t50 * qJD(5) + t100 * t46 + t134 * t47;
t40 = -t50 * mrSges(7,2) + t93 * mrSges(7,3);
t41 = -t93 * mrSges(6,2) - t50 * mrSges(6,3);
t10 = m(6) * t116 - t136 + (t41 + t40) * t93 + (mrSges(6,1) + mrSges(7,1)) * t81 + t131 * t51 + t133 * t26;
t52 = -t69 * mrSges(5,1) + t70 * mrSges(5,2);
t59 = -t94 * mrSges(5,2) + t69 * mrSges(5,3);
t132 = t100 * t18 + t134 * t20;
t43 = -t93 * mrSges(7,1) + t51 * mrSges(7,2);
t125 = m(7) * (-t92 * pkin(5) + t81 * qJ(6) + 0.2e1 * qJD(6) * t93 - t50 * t34 + t132) + t93 * t43 + t81 * mrSges(7,3);
t25 = t51 * qJD(5) + t100 * t47 - t134 * t46;
t42 = t93 * mrSges(6,1) - t51 * mrSges(6,3);
t9 = m(6) * t132 - t81 * mrSges(6,2) + t131 * t50 + t133 * t25 - t93 * t42 + t125;
t7 = m(5) * t119 + t83 * mrSges(5,1) - t47 * mrSges(5,3) + t134 * t10 + t100 * t9 - t70 * t52 + t94 * t59;
t71 = -t84 * mrSges(4,1) + t85 * mrSges(4,2);
t72 = -t94 * mrSges(4,2) + t84 * mrSges(4,3);
t60 = t94 * mrSges(5,1) - t70 * mrSges(5,3);
t8 = m(5) * t124 - t83 * mrSges(5,2) + t46 * mrSges(5,3) - t100 * t10 + t134 * t9 + t69 * t52 - t94 * t60;
t5 = m(4) * t120 + t83 * mrSges(4,1) - t67 * mrSges(4,3) + t99 * t7 - t85 * t71 + t94 * t72 + t98 * t8;
t74 = t94 * mrSges(4,1) - t85 * mrSges(4,3);
t6 = m(4) * t130 - t83 * mrSges(4,2) + t66 * mrSges(4,3) - t98 * t7 + t84 * t71 - t94 * t74 + t99 * t8;
t90 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t128;
t91 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t127;
t137 = m(3) * t79 - t89 * mrSges(3,1) + t88 * mrSges(3,2) + t101 * t6 + t104 * t5 + (t102 * t90 - t105 * t91) * qJD(1);
t129 = -t105 * g(3) - t102 * t80;
t57 = -qJDD(2) * pkin(2) - t107 * pkin(8) + t87 * t128 - t129;
t113 = -t66 * pkin(3) - t82 * qJ(4) + t85 * t73 + qJDD(4) + t57;
t111 = -t46 * pkin(4) - t68 * pkin(9) + t70 * t61 + t113;
t115 = t26 * mrSges(7,3) + t51 * t43 - m(7) * (t111 + (t50 * t93 - t26) * qJ(6) + (t51 * t93 + t25) * pkin(5) - 0.2e1 * qJD(6) * t51) - t25 * mrSges(7,1) - t50 * t40;
t112 = m(6) * t111 + t25 * mrSges(6,1) + t26 * mrSges(6,2) + t50 * t41 + t51 * t42 - t115;
t110 = -m(5) * t113 + t46 * mrSges(5,1) - t47 * mrSges(5,2) + t69 * t59 - t70 * t60 - t112;
t109 = m(4) * t57 - t66 * mrSges(4,1) + t67 * mrSges(4,2) - t84 * t72 + t85 * t74 - t110;
t86 = (-mrSges(3,1) * t105 + mrSges(3,2) * t102) * qJD(1);
t12 = m(3) * t129 + qJDD(2) * mrSges(3,1) - t88 * mrSges(3,3) + qJD(2) * t91 - t86 * t128 - t109;
t4 = m(3) * t122 - qJDD(2) * mrSges(3,2) + t89 * mrSges(3,3) - qJD(2) * t90 - t101 * t5 + t104 * t6 + t86 * t127;
t135 = t102 * t4 + t105 * t12;
t2 = m(2) * t123 + qJDD(1) * mrSges(2,1) - t108 * mrSges(2,2) - t137;
t1 = m(2) * t118 - t108 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t102 * t12 + t105 * t4;
t3 = [-m(1) * g(1) + t106 * t1 - t103 * t2, t1, t4, t6, t8, t9, -t25 * mrSges(7,2) - t50 * t35 + t125; -m(1) * g(2) + t103 * t1 + t106 * t2, t2, t12, t5, t7, t10, -t115; (-m(1) - m(2)) * g(3) + t135, -m(2) * g(3) + t135, t137, t109, -t110, t112, -t81 * mrSges(7,1) + t26 * mrSges(7,2) + t51 * t35 - t93 * t40 + t136;];
f_new  = t3;
