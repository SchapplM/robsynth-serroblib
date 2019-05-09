% Calculate vector of cutting forces with Newton-Euler
% S6RRRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-05-07 19:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRPR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 19:54:39
% EndTime: 2019-05-07 19:54:47
% DurationCPUTime: 3.47s
% Computational Cost: add. (42593->205), mult. (93346->259), div. (0->0), fcn. (68202->10), ass. (0->102)
t101 = sin(qJ(2));
t105 = cos(qJ(2));
t107 = qJD(1) ^ 2;
t102 = sin(qJ(1));
t106 = cos(qJ(1));
t125 = t102 * g(1) - t106 * g(2);
t119 = -qJDD(1) * pkin(1) - t125;
t128 = qJD(1) * t101;
t126 = qJD(1) * qJD(2);
t87 = qJDD(1) * t105 - t101 * t126;
t90 = qJD(2) * pkin(2) - pkin(8) * t128;
t97 = t105 ^ 2;
t112 = -t87 * pkin(2) + t90 * t128 + (-pkin(8) * t97 - pkin(7)) * t107 + t119;
t100 = sin(qJ(3));
t104 = cos(qJ(3));
t81 = (t100 * t105 + t101 * t104) * qJD(1);
t86 = qJDD(1) * t101 + t105 * t126;
t57 = -qJD(3) * t81 - t100 * t86 + t104 * t87;
t96 = qJD(2) + qJD(3);
t76 = pkin(3) * t96 - pkin(9) * t81;
t80 = (-t100 * t101 + t104 * t105) * qJD(1);
t79 = t80 ^ 2;
t110 = -t57 * pkin(3) - pkin(9) * t79 + t81 * t76 + t112;
t103 = cos(qJ(6));
t142 = cos(qJ(4));
t99 = sin(qJ(4));
t71 = -t142 * t80 + t81 * t99;
t93 = qJD(4) + t96;
t139 = t71 * t93;
t145 = -2 * qJD(5);
t58 = qJD(3) * t80 + t100 * t87 + t104 * t86;
t37 = -t71 * qJD(4) + t142 * t58 + t99 * t57;
t72 = t142 * t81 + t99 * t80;
t108 = (-t37 + t139) * qJ(5) + t110 + (t93 * pkin(4) + t145) * t72;
t121 = -g(1) * t106 - g(2) * t102;
t83 = -pkin(1) * t107 + qJDD(1) * pkin(7) + t121;
t129 = t101 * t83;
t141 = pkin(2) * t107;
t50 = qJDD(2) * pkin(2) - t86 * pkin(8) - t129 + (pkin(8) * t126 + t101 * t141 - g(3)) * t105;
t124 = -g(3) * t101 + t105 * t83;
t51 = pkin(8) * t87 - qJD(2) * t90 - t97 * t141 + t124;
t123 = -t100 * t51 + t104 * t50;
t95 = qJDD(2) + qJDD(3);
t23 = (t80 * t96 - t58) * pkin(9) + (t80 * t81 + t95) * pkin(3) + t123;
t131 = t100 * t50 + t104 * t51;
t27 = -pkin(3) * t79 + t57 * pkin(9) - t76 * t96 + t131;
t122 = t142 * t23 - t99 * t27;
t44 = pkin(4) * t71 - qJ(5) * t72;
t91 = t93 ^ 2;
t92 = qJDD(4) + t95;
t20 = -t92 * pkin(4) - t91 * qJ(5) + t72 * t44 + qJDD(5) - t122;
t14 = (t71 * t72 - t92) * pkin(10) + (t37 + t139) * pkin(5) + t20;
t36 = t72 * qJD(4) - t142 * t57 + t58 * t99;
t64 = t72 * pkin(5) - pkin(10) * t93;
t70 = t71 ^ 2;
t15 = -t72 * t64 - t70 * pkin(5) + t108 + (pkin(4) + pkin(10)) * t36;
t98 = sin(qJ(6));
t53 = t103 * t71 - t93 * t98;
t29 = t53 * qJD(6) + t103 * t92 + t36 * t98;
t34 = qJDD(6) + t37;
t54 = t103 * t93 + t71 * t98;
t39 = -mrSges(7,1) * t53 + mrSges(7,2) * t54;
t67 = qJD(6) + t72;
t40 = -mrSges(7,2) * t67 + mrSges(7,3) * t53;
t12 = m(7) * (t103 * t14 - t15 * t98) - t29 * mrSges(7,3) + t34 * mrSges(7,1) - t54 * t39 + t67 * t40;
t28 = -t54 * qJD(6) + t103 * t36 - t92 * t98;
t41 = mrSges(7,1) * t67 - mrSges(7,3) * t54;
t13 = m(7) * (t103 * t15 + t14 * t98) + t28 * mrSges(7,3) - t34 * mrSges(7,2) + t53 * t39 - t67 * t41;
t61 = t72 * mrSges(6,1) + mrSges(6,2) * t93;
t118 = -t103 * t13 + t98 * t12 - m(6) * (t36 * pkin(4) + t108) + t37 * mrSges(6,3) + t72 * t61;
t60 = t71 * mrSges(6,1) - mrSges(6,3) * t93;
t130 = -mrSges(5,2) * t93 - t71 * mrSges(5,3) - t60;
t135 = mrSges(6,2) - mrSges(5,1);
t63 = mrSges(5,1) * t93 - t72 * mrSges(5,3);
t111 = m(5) * t110 + t37 * mrSges(5,2) + t130 * t71 - t135 * t36 + t72 * t63 - t118;
t74 = -mrSges(4,2) * t96 + mrSges(4,3) * t80;
t75 = mrSges(4,1) * t96 - mrSges(4,3) * t81;
t109 = m(4) * t112 - t57 * mrSges(4,1) + t58 * mrSges(4,2) - t80 * t74 + t81 * t75 + t111;
t88 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t128;
t127 = qJD(1) * t105;
t89 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t127;
t147 = t109 + (t101 * t88 - t105 * t89) * qJD(1) - t87 * mrSges(3,1) + t86 * mrSges(3,2) + m(3) * (-pkin(7) * t107 + t119);
t133 = t142 * t27 + t99 * t23;
t114 = -t91 * pkin(4) + t92 * qJ(5) - t71 * t44 + t133;
t116 = -t28 * mrSges(7,1) - t53 * t40 + m(7) * (-t36 * pkin(5) - t70 * pkin(10) + ((2 * qJD(5)) + t64) * t93 + t114) + t29 * mrSges(7,2) + t54 * t41;
t113 = -m(6) * (t93 * t145 - t114) + t116;
t46 = -mrSges(6,2) * t71 - mrSges(6,3) * t72;
t132 = -mrSges(5,1) * t71 - mrSges(5,2) * t72 - t46;
t134 = -mrSges(5,3) - mrSges(6,1);
t10 = m(5) * t133 + (-t63 + t61) * t93 + (-mrSges(5,2) + mrSges(6,3)) * t92 + t132 * t71 + t134 * t36 + t113;
t73 = -mrSges(4,1) * t80 + mrSges(4,2) * t81;
t117 = -m(6) * t20 - t103 * t12 - t98 * t13;
t9 = m(5) * t122 + t130 * t93 + t132 * t72 + t134 * t37 - t135 * t92 + t117;
t6 = m(4) * t123 + t95 * mrSges(4,1) - t58 * mrSges(4,3) + t99 * t10 + t142 * t9 - t81 * t73 + t96 * t74;
t7 = m(4) * t131 - t95 * mrSges(4,2) + t57 * mrSges(4,3) + t142 * t10 + t80 * t73 - t96 * t75 - t99 * t9;
t85 = (-mrSges(3,1) * t105 + mrSges(3,2) * t101) * qJD(1);
t4 = m(3) * (-t105 * g(3) - t129) - t86 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t85 * t128 + qJD(2) * t89 + t100 * t7 + t104 * t6;
t5 = m(3) * t124 - qJDD(2) * mrSges(3,2) + t87 * mrSges(3,3) - qJD(2) * t88 - t100 * t6 + t104 * t7 + t85 * t127;
t144 = t101 * t5 + t105 * t4;
t8 = m(2) * t125 + qJDD(1) * mrSges(2,1) - t107 * mrSges(2,2) - t147;
t1 = m(2) * t121 - t107 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t101 * t4 + t105 * t5;
t2 = [-m(1) * g(1) + t1 * t106 - t102 * t8, t1, t5, t7, t10, -t36 * mrSges(6,2) - t71 * t60 - t118, t13; -m(1) * g(2) + t1 * t102 + t106 * t8, t8, t4, t6, t9, t36 * mrSges(6,1) - t92 * mrSges(6,3) + t71 * t46 - t93 * t61 - t113, t12; (-m(1) - m(2)) * g(3) + t144, -m(2) * g(3) + t144, t147, t109, t111, t37 * mrSges(6,1) + t92 * mrSges(6,2) + t72 * t46 + t93 * t60 - t117, t116;];
f_new  = t2;
