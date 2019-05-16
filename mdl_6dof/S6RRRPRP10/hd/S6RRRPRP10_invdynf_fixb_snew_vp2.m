% Calculate vector of cutting forces with Newton-Euler
% S6RRRPRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-05-07 09:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPRP10_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP10_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP10_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP10_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP10_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP10_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP10_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 08:50:13
% EndTime: 2019-05-07 08:50:28
% DurationCPUTime: 5.24s
% Computational Cost: add. (92239->209), mult. (201297->282), div. (0->0), fcn. (158957->12), ass. (0->107)
t103 = sin(pkin(6));
t108 = sin(qJ(2));
t111 = cos(qJ(2));
t130 = qJD(1) * qJD(2);
t92 = (-qJDD(1) * t111 + t108 * t130) * t103;
t106 = sin(qJ(5));
t144 = cos(qJ(5));
t102 = sin(pkin(11));
t104 = cos(pkin(11));
t107 = sin(qJ(3));
t110 = cos(qJ(3));
t132 = qJD(1) * t111;
t105 = cos(pkin(6));
t113 = qJD(1) ^ 2;
t109 = sin(qJ(1));
t112 = cos(qJ(1));
t125 = t109 * g(1) - t112 * g(2);
t143 = pkin(8) * t103;
t87 = qJDD(1) * pkin(1) + t113 * t143 + t125;
t136 = t105 * t87;
t122 = -t112 * g(1) - t109 * g(2);
t88 = -t113 * pkin(1) + qJDD(1) * t143 + t122;
t137 = t108 * t136 + t111 * t88;
t133 = qJD(1) * t103;
t90 = (-pkin(2) * t111 - pkin(9) * t108) * t133;
t99 = t105 * qJD(1) + qJD(2);
t97 = t99 ^ 2;
t98 = t105 * qJDD(1) + qJDD(2);
t50 = -t97 * pkin(2) + t98 * pkin(9) + (-g(3) * t108 + t90 * t132) * t103 + t137;
t142 = t105 * g(3);
t91 = (qJDD(1) * t108 + t111 * t130) * t103;
t51 = t92 * pkin(2) - t91 * pkin(9) - t142 + (-t87 + (pkin(2) * t108 - pkin(9) * t111) * t99 * qJD(1)) * t103;
t138 = t107 * t51 + t110 * t50;
t127 = t108 * t133;
t80 = t107 * t127 - t110 * t99;
t81 = t107 * t99 + t110 * t127;
t66 = t80 * pkin(3) - t81 * qJ(4);
t84 = qJDD(3) + t92;
t126 = t103 * t132;
t96 = qJD(3) - t126;
t95 = t96 ^ 2;
t27 = -t95 * pkin(3) + t84 * qJ(4) - t80 * t66 + t138;
t134 = t103 * t111;
t121 = -g(3) * t134 - t108 * t88 + t111 * t136;
t49 = -t98 * pkin(2) - t97 * pkin(9) + t90 * t127 - t121;
t64 = t81 * qJD(3) + t107 * t91 - t110 * t98;
t65 = -t80 * qJD(3) + t107 * t98 + t110 * t91;
t30 = (t80 * t96 - t65) * qJ(4) + (t81 * t96 + t64) * pkin(3) + t49;
t73 = t102 * t96 + t104 * t81;
t123 = -0.2e1 * qJD(4) * t73 - t102 * t27 + t104 * t30;
t56 = t102 * t84 + t104 * t65;
t72 = -t102 * t81 + t104 * t96;
t20 = (t72 * t80 - t56) * pkin(10) + (t72 * t73 + t64) * pkin(4) + t123;
t128 = 0.2e1 * qJD(4) * t72 + t102 * t30 + t104 * t27;
t55 = -t102 * t65 + t104 * t84;
t61 = t80 * pkin(4) - t73 * pkin(10);
t71 = t72 ^ 2;
t22 = -t71 * pkin(4) + t55 * pkin(10) - t80 * t61 + t128;
t120 = -t106 * t22 + t144 * t20;
t53 = t106 * t73 - t144 * t72;
t54 = t106 * t72 + t144 * t73;
t37 = t53 * pkin(5) - t54 * qJ(6);
t63 = qJDD(5) + t64;
t79 = qJD(5) + t80;
t78 = t79 ^ 2;
t145 = m(7) * (-t63 * pkin(5) - t78 * qJ(6) + t54 * t37 + qJDD(6) - t120);
t141 = -mrSges(6,3) - mrSges(7,2);
t140 = t106 * t20 + t144 * t22;
t38 = t53 * mrSges(7,1) - t54 * mrSges(7,3);
t139 = -t53 * mrSges(6,1) - t54 * mrSges(6,2) - t38;
t135 = t103 * t108;
t124 = -t107 * t50 + t110 * t51;
t26 = -t84 * pkin(3) - t95 * qJ(4) + t81 * t66 + qJDD(4) - t124;
t116 = -t55 * pkin(4) - t71 * pkin(10) + t73 * t61 + t26;
t33 = t54 * qJD(5) + t106 * t56 - t144 * t55;
t34 = -t53 * qJD(5) + t106 * t55 + t144 * t56;
t41 = -t53 * mrSges(7,2) + t79 * mrSges(7,3);
t44 = -t79 * mrSges(7,1) + t54 * mrSges(7,2);
t118 = t34 * mrSges(7,3) + t54 * t44 - m(7) * (-0.2e1 * qJD(6) * t54 + (t53 * t79 - t34) * qJ(6) + (t54 * t79 + t33) * pkin(5) + t116) - t33 * mrSges(7,1) - t53 * t41;
t42 = -t79 * mrSges(6,2) - t53 * mrSges(6,3);
t43 = t79 * mrSges(6,1) - t54 * mrSges(6,3);
t117 = m(6) * t116 + t33 * mrSges(6,1) + t34 * mrSges(6,2) + t53 * t42 + t54 * t43 - t118;
t59 = -t80 * mrSges(5,2) + t72 * mrSges(5,3);
t60 = t80 * mrSges(5,1) - t73 * mrSges(5,3);
t114 = m(5) * t26 - t55 * mrSges(5,1) + t56 * mrSges(5,2) - t72 * t59 + t73 * t60 + t117;
t67 = t80 * mrSges(4,1) + t81 * mrSges(4,2);
t74 = -t96 * mrSges(4,2) - t80 * mrSges(4,3);
t14 = m(4) * t124 + t84 * mrSges(4,1) - t65 * mrSges(4,3) - t81 * t67 + t96 * t74 - t114;
t85 = t99 * mrSges(3,1) - mrSges(3,3) * t127;
t89 = (-mrSges(3,1) * t111 + mrSges(3,2) * t108) * t133;
t129 = m(7) * (-t78 * pkin(5) + t63 * qJ(6) + 0.2e1 * qJD(6) * t79 - t53 * t37 + t140) + t79 * t44 + t63 * mrSges(7,3);
t12 = m(6) * t140 - t63 * mrSges(6,2) + t139 * t53 + t141 * t33 - t79 * t43 + t129;
t13 = m(6) * t120 - t145 + (t42 + t41) * t79 + (mrSges(6,1) + mrSges(7,1)) * t63 + t139 * t54 + t141 * t34;
t57 = -t72 * mrSges(5,1) + t73 * mrSges(5,2);
t10 = m(5) * t123 + t64 * mrSges(5,1) - t56 * mrSges(5,3) + t106 * t12 + t144 * t13 - t73 * t57 + t80 * t59;
t11 = m(5) * t128 - t64 * mrSges(5,2) + t55 * mrSges(5,3) - t106 * t13 + t144 * t12 + t72 * t57 - t80 * t60;
t75 = t96 * mrSges(4,1) - t81 * mrSges(4,3);
t9 = m(4) * t138 - t84 * mrSges(4,2) - t64 * mrSges(4,3) - t102 * t10 + t104 * t11 - t80 * t67 - t96 * t75;
t4 = m(3) * (-g(3) * t135 + t137) - t92 * mrSges(3,3) - t98 * mrSges(3,2) + t89 * t126 - t99 * t85 + t110 * t9 - t107 * t14;
t86 = -t99 * mrSges(3,2) + mrSges(3,3) * t126;
t6 = m(3) * (-t103 * t87 - t142) + t91 * mrSges(3,2) + t92 * mrSges(3,1) + t107 * t9 + t110 * t14 + (t108 * t85 - t111 * t86) * t133;
t115 = m(4) * t49 + t64 * mrSges(4,1) + t65 * mrSges(4,2) + t104 * t10 + t102 * t11 + t80 * t74 + t81 * t75;
t8 = m(3) * t121 + t98 * mrSges(3,1) - t91 * mrSges(3,3) - t89 * t127 + t99 * t86 - t115;
t131 = t105 * t6 + t8 * t134 + t4 * t135;
t2 = m(2) * t122 - t113 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t108 * t8 + t111 * t4;
t1 = m(2) * t125 + qJDD(1) * mrSges(2,1) - t113 * mrSges(2,2) - t103 * t6 + (t108 * t4 + t111 * t8) * t105;
t3 = [-m(1) * g(1) - t109 * t1 + t112 * t2, t2, t4, t9, t11, t12, -t33 * mrSges(7,2) - t53 * t38 + t129; -m(1) * g(2) + t112 * t1 + t109 * t2, t1, t8, t14, t10, t13, -t118; (-m(1) - m(2)) * g(3) + t131, -m(2) * g(3) + t131, t6, t115, t114, t117, -t63 * mrSges(7,1) + t34 * mrSges(7,2) + t54 * t38 - t79 * t41 + t145;];
f_new  = t3;
