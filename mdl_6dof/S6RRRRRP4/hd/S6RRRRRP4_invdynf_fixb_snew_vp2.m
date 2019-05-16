% Calculate vector of cutting forces with Newton-Euler
% S6RRRRRP4
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
% Datum: 2019-05-08 04:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRRP4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 04:46:55
% EndTime: 2019-05-08 04:47:07
% DurationCPUTime: 3.84s
% Computational Cost: add. (54372->203), mult. (108391->259), div. (0->0), fcn. (78086->10), ass. (0->99)
t104 = sin(qJ(2));
t108 = cos(qJ(2));
t110 = qJD(1) ^ 2;
t101 = sin(qJ(5));
t137 = cos(qJ(5));
t102 = sin(qJ(4));
t106 = cos(qJ(4));
t100 = t108 ^ 2;
t105 = sin(qJ(1));
t109 = cos(qJ(1));
t125 = t105 * g(1) - t109 * g(2);
t119 = -qJDD(1) * pkin(1) - t125;
t129 = qJD(1) * t104;
t127 = qJD(1) * qJD(2);
t91 = t108 * qJDD(1) - t104 * t127;
t94 = qJD(2) * pkin(2) - pkin(8) * t129;
t114 = -t91 * pkin(2) + t94 * t129 + (-pkin(8) * t100 - pkin(7)) * t110 + t119;
t103 = sin(qJ(3));
t107 = cos(qJ(3));
t84 = (t103 * t108 + t104 * t107) * qJD(1);
t90 = t104 * qJDD(1) + t108 * t127;
t63 = -t84 * qJD(3) - t103 * t90 + t107 * t91;
t128 = qJD(1) * t108;
t83 = -t103 * t129 + t107 * t128;
t64 = t83 * qJD(3) + t103 * t91 + t107 * t90;
t99 = qJD(2) + qJD(3);
t31 = (-t83 * t99 - t64) * pkin(9) + (t84 * t99 - t63) * pkin(3) + t114;
t121 = -t109 * g(1) - t105 * g(2);
t86 = -t110 * pkin(1) + qJDD(1) * pkin(7) + t121;
t130 = t104 * t86;
t136 = pkin(2) * t110;
t55 = qJDD(2) * pkin(2) - t90 * pkin(8) - t130 + (pkin(8) * t127 + t104 * t136 - g(3)) * t108;
t124 = -t104 * g(3) + t108 * t86;
t56 = t91 * pkin(8) - qJD(2) * t94 - t100 * t136 + t124;
t131 = t103 * t55 + t107 * t56;
t71 = -t83 * pkin(3) - t84 * pkin(9);
t97 = t99 ^ 2;
t98 = qJDD(2) + qJDD(3);
t34 = -t97 * pkin(3) + t98 * pkin(9) + t83 * t71 + t131;
t123 = -t102 * t34 + t106 * t31;
t74 = -t102 * t84 + t106 * t99;
t43 = t74 * qJD(4) + t102 * t98 + t106 * t64;
t61 = qJDD(4) - t63;
t75 = t102 * t99 + t106 * t84;
t82 = qJD(4) - t83;
t20 = (t74 * t82 - t43) * pkin(10) + (t74 * t75 + t61) * pkin(4) + t123;
t133 = t102 * t31 + t106 * t34;
t42 = -t75 * qJD(4) - t102 * t64 + t106 * t98;
t68 = t82 * pkin(4) - t75 * pkin(10);
t73 = t74 ^ 2;
t22 = -t73 * pkin(4) + t42 * pkin(10) - t82 * t68 + t133;
t134 = t101 * t20 + t137 * t22;
t49 = t101 * t75 - t137 * t74;
t50 = t101 * t74 + t137 * t75;
t37 = t49 * pkin(5) - t50 * qJ(6);
t80 = qJD(5) + t82;
t47 = -t80 * mrSges(7,1) + t50 * mrSges(7,2);
t58 = qJDD(5) + t61;
t78 = t80 ^ 2;
t126 = m(7) * (-t78 * pkin(5) + t58 * qJ(6) + 0.2e1 * qJD(6) * t80 - t49 * t37 + t134) + t80 * t47 + t58 * mrSges(7,3);
t38 = t49 * mrSges(7,1) - t50 * mrSges(7,3);
t132 = -t49 * mrSges(6,1) - t50 * mrSges(6,2) - t38;
t135 = -mrSges(6,3) - mrSges(7,2);
t27 = t50 * qJD(5) + t101 * t43 - t137 * t42;
t46 = t80 * mrSges(6,1) - t50 * mrSges(6,3);
t12 = m(6) * t134 - t58 * mrSges(6,2) + t132 * t49 + t135 * t27 - t80 * t46 + t126;
t118 = -t101 * t22 + t137 * t20;
t138 = m(7) * (-t58 * pkin(5) - t78 * qJ(6) + t50 * t37 + qJDD(6) - t118);
t28 = -t49 * qJD(5) + t101 * t42 + t137 * t43;
t44 = -t49 * mrSges(7,2) + t80 * mrSges(7,3);
t45 = -t80 * mrSges(6,2) - t49 * mrSges(6,3);
t13 = m(6) * t118 - t138 + (t45 + t44) * t80 + (mrSges(6,1) + mrSges(7,1)) * t58 + t132 * t50 + t135 * t28;
t54 = -t74 * mrSges(5,1) + t75 * mrSges(5,2);
t66 = -t82 * mrSges(5,2) + t74 * mrSges(5,3);
t10 = m(5) * t123 + t61 * mrSges(5,1) - t43 * mrSges(5,3) + t101 * t12 + t137 * t13 - t75 * t54 + t82 * t66;
t67 = t82 * mrSges(5,1) - t75 * mrSges(5,3);
t11 = m(5) * t133 - t61 * mrSges(5,2) + t42 * mrSges(5,3) - t101 * t13 + t137 * t12 + t74 * t54 - t82 * t67;
t76 = -t99 * mrSges(4,2) + t83 * mrSges(4,3);
t77 = t99 * mrSges(4,1) - t84 * mrSges(4,3);
t116 = -m(4) * t114 + t63 * mrSges(4,1) - t64 * mrSges(4,2) - t106 * t10 - t102 * t11 + t83 * t76 - t84 * t77;
t92 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t129;
t93 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t128;
t140 = (t104 * t92 - t108 * t93) * qJD(1) + m(3) * (-t110 * pkin(7) + t119) - t91 * mrSges(3,1) + t90 * mrSges(3,2) - t116;
t122 = -t103 * t56 + t107 * t55;
t33 = -t98 * pkin(3) - t97 * pkin(9) + t84 * t71 - t122;
t115 = -t42 * pkin(4) - t73 * pkin(10) + t75 * t68 + t33;
t117 = t28 * mrSges(7,3) + t50 * t47 - m(7) * (-0.2e1 * qJD(6) * t50 + (t49 * t80 - t28) * qJ(6) + (t50 * t80 + t27) * pkin(5) + t115) - t27 * mrSges(7,1) - t49 * t44;
t113 = m(6) * t115 + t27 * mrSges(6,1) + t28 * mrSges(6,2) + t49 * t45 + t50 * t46 - t117;
t111 = m(5) * t33 - t42 * mrSges(5,1) + t43 * mrSges(5,2) - t74 * t66 + t75 * t67 + t113;
t70 = -t83 * mrSges(4,1) + t84 * mrSges(4,2);
t14 = m(4) * t122 + t98 * mrSges(4,1) - t64 * mrSges(4,3) - t84 * t70 + t99 * t76 - t111;
t7 = m(4) * t131 - t98 * mrSges(4,2) + t63 * mrSges(4,3) - t102 * t10 + t106 * t11 + t83 * t70 - t99 * t77;
t89 = (-mrSges(3,1) * t108 + mrSges(3,2) * t104) * qJD(1);
t4 = m(3) * (-t108 * g(3) - t130) - t90 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t89 * t129 + qJD(2) * t93 + t103 * t7 + t107 * t14;
t5 = m(3) * t124 - qJDD(2) * mrSges(3,2) + t91 * mrSges(3,3) - qJD(2) * t92 - t103 * t14 + t107 * t7 + t89 * t128;
t139 = t104 * t5 + t108 * t4;
t6 = m(2) * t125 + qJDD(1) * mrSges(2,1) - t110 * mrSges(2,2) - t140;
t1 = m(2) * t121 - t110 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t104 * t4 + t108 * t5;
t2 = [-m(1) * g(1) + t109 * t1 - t105 * t6, t1, t5, t7, t11, t12, -t27 * mrSges(7,2) - t49 * t38 + t126; -m(1) * g(2) + t105 * t1 + t109 * t6, t6, t4, t14, t10, t13, -t117; (-m(1) - m(2)) * g(3) + t139, -m(2) * g(3) + t139, t140, -t116, t111, t113, -t58 * mrSges(7,1) + t28 * mrSges(7,2) + t50 * t38 - t80 * t44 + t138;];
f_new  = t2;
