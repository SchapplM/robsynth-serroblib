% Calculate vector of cutting forces with Newton-Euler
% S6RRRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 04:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRRP2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 04:29:20
% EndTime: 2019-05-08 04:29:33
% DurationCPUTime: 4.13s
% Computational Cost: add. (54416->205), mult. (117803->260), div. (0->0), fcn. (87409->10), ass. (0->99)
t101 = sin(qJ(2));
t105 = cos(qJ(2));
t107 = qJD(1) ^ 2;
t102 = sin(qJ(1));
t106 = cos(qJ(1));
t121 = t102 * g(1) - t106 * g(2);
t115 = -qJDD(1) * pkin(1) - t121;
t126 = qJD(1) * t101;
t124 = qJD(1) * qJD(2);
t85 = t105 * qJDD(1) - t101 * t124;
t88 = qJD(2) * pkin(2) - pkin(8) * t126;
t97 = t105 ^ 2;
t111 = -t85 * pkin(2) + t88 * t126 + (-pkin(8) * t97 - pkin(7)) * t107 + t115;
t100 = sin(qJ(3));
t104 = cos(qJ(3));
t79 = (t100 * t105 + t101 * t104) * qJD(1);
t84 = t101 * qJDD(1) + t105 * t124;
t61 = -t79 * qJD(3) - t100 * t84 + t104 * t85;
t96 = qJD(2) + qJD(3);
t75 = t96 * pkin(3) - t79 * pkin(9);
t78 = (-t100 * t101 + t104 * t105) * qJD(1);
t77 = t78 ^ 2;
t109 = -t61 * pkin(3) - t77 * pkin(9) + t79 * t75 + t111;
t136 = cos(qJ(5));
t103 = cos(qJ(4));
t117 = -t106 * g(1) - t102 * g(2);
t81 = -t107 * pkin(1) + qJDD(1) * pkin(7) + t117;
t127 = t101 * t81;
t135 = pkin(2) * t107;
t55 = qJDD(2) * pkin(2) - t84 * pkin(8) - t127 + (pkin(8) * t124 + t101 * t135 - g(3)) * t105;
t120 = -t101 * g(3) + t105 * t81;
t56 = t85 * pkin(8) - qJD(2) * t88 - t97 * t135 + t120;
t119 = -t100 * t56 + t104 * t55;
t62 = t78 * qJD(3) + t100 * t85 + t104 * t84;
t95 = qJDD(2) + qJDD(3);
t25 = (t78 * t96 - t62) * pkin(9) + (t78 * t79 + t95) * pkin(3) + t119;
t128 = t100 * t55 + t104 * t56;
t29 = -t77 * pkin(3) + t61 * pkin(9) - t96 * t75 + t128;
t99 = sin(qJ(4));
t131 = t103 * t29 + t99 * t25;
t70 = t103 * t78 - t99 * t79;
t71 = t103 * t79 + t99 * t78;
t51 = -t70 * pkin(4) - t71 * pkin(10);
t93 = qJD(4) + t96;
t91 = t93 ^ 2;
t92 = qJDD(4) + t95;
t21 = -t91 * pkin(4) + t92 * pkin(10) + t70 * t51 + t131;
t38 = -t71 * qJD(4) + t103 * t61 - t99 * t62;
t39 = t70 * qJD(4) + t103 * t62 + t99 * t61;
t23 = (-t70 * t93 - t39) * pkin(10) + (t71 * t93 - t38) * pkin(4) + t109;
t98 = sin(qJ(5));
t132 = t136 * t21 + t98 * t23;
t36 = qJDD(5) - t38;
t57 = -t136 * t93 + t98 * t71;
t58 = t136 * t71 + t98 * t93;
t42 = t57 * pkin(5) - t58 * qJ(6);
t67 = qJD(5) - t70;
t48 = -t67 * mrSges(7,1) + t58 * mrSges(7,2);
t66 = t67 ^ 2;
t123 = m(7) * (-t66 * pkin(5) + t36 * qJ(6) + 0.2e1 * qJD(6) * t67 - t57 * t42 + t132) + t36 * mrSges(7,3) + t67 * t48;
t43 = t57 * mrSges(7,1) - t58 * mrSges(7,3);
t130 = -t57 * mrSges(6,1) - t58 * mrSges(6,2) - t43;
t133 = -mrSges(6,3) - mrSges(7,2);
t30 = t58 * qJD(5) - t136 * t92 + t98 * t39;
t47 = t67 * mrSges(6,1) - t58 * mrSges(6,3);
t12 = m(6) * t132 - t36 * mrSges(6,2) + t130 * t57 + t133 * t30 - t67 * t47 + t123;
t114 = t136 * t23 - t98 * t21;
t137 = m(7) * (-t36 * pkin(5) - t66 * qJ(6) + t58 * t42 + qJDD(6) - t114);
t31 = -t57 * qJD(5) + t136 * t39 + t98 * t92;
t45 = -t57 * mrSges(7,2) + t67 * mrSges(7,3);
t46 = -t67 * mrSges(6,2) - t57 * mrSges(6,3);
t14 = m(6) * t114 - t137 + (t46 + t45) * t67 + t130 * t58 + (mrSges(6,1) + mrSges(7,1)) * t36 + t133 * t31;
t64 = -t93 * mrSges(5,2) + t70 * mrSges(5,3);
t65 = t93 * mrSges(5,1) - t71 * mrSges(5,3);
t113 = -m(5) * t109 + t38 * mrSges(5,1) - t39 * mrSges(5,2) - t98 * t12 - t136 * t14 + t70 * t64 - t71 * t65;
t73 = -t96 * mrSges(4,2) + t78 * mrSges(4,3);
t74 = t96 * mrSges(4,1) - t79 * mrSges(4,3);
t110 = -m(4) * t111 + t61 * mrSges(4,1) - t62 * mrSges(4,2) + t78 * t73 - t79 * t74 + t113;
t86 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t126;
t125 = qJD(1) * t105;
t87 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t125;
t140 = (t101 * t86 - t105 * t87) * qJD(1) + m(3) * (-t107 * pkin(7) + t115) - t85 * mrSges(3,1) + t84 * mrSges(3,2) - t110;
t118 = t103 * t25 - t99 * t29;
t20 = -t92 * pkin(4) - t91 * pkin(10) + t71 * t51 - t118;
t122 = m(7) * (-0.2e1 * qJD(6) * t58 + (t57 * t67 - t31) * qJ(6) + (t58 * t67 + t30) * pkin(5) + t20) + t30 * mrSges(7,1) + t57 * t45;
t139 = m(6) * t20 + t30 * mrSges(6,1) + (t47 - t48) * t58 + (mrSges(6,2) - mrSges(7,3)) * t31 + t57 * t46 + t122;
t50 = -t70 * mrSges(5,1) + t71 * mrSges(5,2);
t10 = m(5) * t118 + t92 * mrSges(5,1) - t39 * mrSges(5,3) - t71 * t50 + t93 * t64 - t139;
t72 = -t78 * mrSges(4,1) + t79 * mrSges(4,2);
t9 = m(5) * t131 - t92 * mrSges(5,2) + t38 * mrSges(5,3) + t136 * t12 - t98 * t14 + t70 * t50 - t93 * t65;
t6 = m(4) * t119 + t95 * mrSges(4,1) - t62 * mrSges(4,3) + t103 * t10 - t79 * t72 + t96 * t73 + t99 * t9;
t7 = m(4) * t128 - t95 * mrSges(4,2) + t61 * mrSges(4,3) - t99 * t10 + t103 * t9 + t78 * t72 - t96 * t74;
t83 = (-mrSges(3,1) * t105 + mrSges(3,2) * t101) * qJD(1);
t4 = m(3) * (-t105 * g(3) - t127) - t84 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t83 * t126 + qJD(2) * t87 + t100 * t7 + t104 * t6;
t5 = m(3) * t120 - qJDD(2) * mrSges(3,2) + t85 * mrSges(3,3) - qJD(2) * t86 - t100 * t6 + t104 * t7 + t83 * t125;
t138 = t101 * t5 + t105 * t4;
t8 = m(2) * t121 + qJDD(1) * mrSges(2,1) - t107 * mrSges(2,2) - t140;
t1 = m(2) * t117 - t107 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t101 * t4 + t105 * t5;
t2 = [-m(1) * g(1) + t106 * t1 - t102 * t8, t1, t5, t7, t9, t12, -t30 * mrSges(7,2) - t57 * t43 + t123; -m(1) * g(2) + t102 * t1 + t106 * t8, t8, t4, t6, t10, t14, -t31 * mrSges(7,3) - t58 * t48 + t122; (-m(1) - m(2)) * g(3) + t138, -m(2) * g(3) + t138, t140, -t110, -t113, t139, -t36 * mrSges(7,1) + t31 * mrSges(7,2) + t58 * t43 - t67 * t45 + t137;];
f_new  = t2;
