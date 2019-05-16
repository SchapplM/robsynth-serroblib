% Calculate vector of cutting forces with Newton-Euler
% S6RRRRPR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
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
% Datum: 2019-05-08 01:41
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRPR13_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR13_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR13_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR13_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR13_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR13_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR13_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 01:19:06
% EndTime: 2019-05-08 01:19:23
% DurationCPUTime: 5.02s
% Computational Cost: add. (70768->214), mult. (149808->280), div. (0->0), fcn. (117957->12), ass. (0->111)
t107 = sin(pkin(6));
t112 = sin(qJ(2));
t116 = cos(qJ(2));
t135 = qJD(1) * qJD(2);
t96 = (-qJDD(1) * t116 + t112 * t135) * t107;
t137 = qJD(1) * t116;
t133 = t107 * t137;
t101 = qJD(3) - t133;
t100 = t101 ^ 2;
t111 = sin(qJ(3));
t115 = cos(qJ(3));
t108 = cos(pkin(6));
t104 = t108 * qJD(1) + qJD(2);
t102 = t104 ^ 2;
t103 = t108 * qJDD(1) + qJDD(2);
t118 = qJD(1) ^ 2;
t113 = sin(qJ(1));
t117 = cos(qJ(1));
t131 = t113 * g(1) - t117 * g(2);
t152 = pkin(8) * t107;
t91 = qJDD(1) * pkin(1) + t118 * t152 + t131;
t141 = t108 * t91;
t128 = -t117 * g(1) - t113 * g(2);
t92 = -t118 * pkin(1) + qJDD(1) * t152 + t128;
t142 = t112 * t141 + t116 * t92;
t138 = qJD(1) * t107;
t94 = (-pkin(2) * t116 - pkin(9) * t112) * t138;
t50 = -t102 * pkin(2) + t103 * pkin(9) + (-g(3) * t112 + t94 * t137) * t107 + t142;
t151 = t108 * g(3);
t95 = (qJDD(1) * t112 + t116 * t135) * t107;
t51 = t96 * pkin(2) - t95 * pkin(9) - t151 + (-t91 + (pkin(2) * t112 - pkin(9) * t116) * t104 * qJD(1)) * t107;
t146 = -t111 * t50 + t115 * t51;
t134 = t112 * t138;
t84 = t115 * t104 - t111 * t134;
t85 = t111 * t104 + t115 * t134;
t71 = -t84 * pkin(3) - t85 * pkin(10);
t88 = qJDD(3) + t96;
t121 = t88 * pkin(3) + t100 * pkin(10) - t85 * t71 + t146;
t110 = sin(qJ(4));
t153 = cos(qJ(4));
t73 = -t153 * t101 + t110 * t85;
t82 = qJD(4) - t84;
t150 = t73 * t82;
t69 = t84 * qJD(3) + t111 * t103 + t115 * t95;
t39 = -t73 * qJD(4) + t110 * t88 + t153 * t69;
t156 = (-t39 + t150) * qJ(5) - t121;
t154 = 2 * qJD(5);
t109 = sin(qJ(6));
t114 = cos(qJ(6));
t74 = t110 * t101 + t153 * t85;
t38 = t74 * qJD(4) + t110 * t69 - t153 * t88;
t54 = t109 * t73 + t114 * t74;
t26 = -t54 * qJD(6) - t109 * t39 + t114 * t38;
t53 = -t109 * t74 + t114 * t73;
t27 = t53 * qJD(6) + t109 * t38 + t114 * t39;
t80 = qJD(6) - t82;
t42 = -t80 * mrSges(7,2) + t53 * mrSges(7,3);
t43 = t80 * mrSges(7,1) - t54 * mrSges(7,3);
t63 = -t82 * pkin(5) - t74 * pkin(11);
t72 = t73 ^ 2;
t130 = m(7) * (-t72 * pkin(11) + (-pkin(4) - pkin(5)) * t38 + (-pkin(4) * t82 + t154 + t63) * t74 - t156) + t27 * mrSges(7,2) - t26 * mrSges(7,1) + t54 * t43 - t53 * t42;
t59 = -t73 * mrSges(6,2) + t82 * mrSges(6,3);
t123 = m(6) * (-0.2e1 * qJD(5) * t74 + (t74 * t82 + t38) * pkin(4) + t156) + t38 * mrSges(6,1) + t73 * t59 - t130;
t60 = -t82 * mrSges(5,2) - t73 * mrSges(5,3);
t61 = t82 * mrSges(5,1) - t74 * mrSges(5,3);
t62 = -t82 * mrSges(6,1) + t74 * mrSges(6,2);
t155 = -m(5) * t121 + t38 * mrSges(5,1) + (t61 - t62) * t74 + (mrSges(5,2) - mrSges(6,3)) * t39 + t73 * t60 + t123;
t148 = -mrSges(5,3) - mrSges(6,2);
t145 = t111 * t51 + t115 * t50;
t31 = -t100 * pkin(3) + t88 * pkin(10) + t84 * t71 + t145;
t139 = t107 * t116;
t127 = -g(3) * t139 - t112 * t92 + t116 * t141;
t49 = -t103 * pkin(2) - t102 * pkin(9) + t94 * t134 - t127;
t68 = -t85 * qJD(3) + t115 * t103 - t111 * t95;
t33 = (-t101 * t84 - t69) * pkin(10) + (t101 * t85 - t68) * pkin(3) + t49;
t147 = t110 * t33 + t153 * t31;
t56 = t73 * mrSges(6,1) - t74 * mrSges(6,3);
t144 = -t73 * mrSges(5,1) - t74 * mrSges(5,2) - t56;
t140 = t107 * t112;
t70 = -t84 * mrSges(4,1) + t85 * mrSges(4,2);
t75 = -t101 * mrSges(4,2) + t84 * mrSges(4,3);
t12 = m(4) * t146 + t88 * mrSges(4,1) - t69 * mrSges(4,3) + t101 * t75 - t85 * t70 - t155;
t89 = t104 * mrSges(3,1) - mrSges(3,3) * t134;
t55 = t73 * pkin(4) - t74 * qJ(5);
t67 = qJDD(4) - t68;
t81 = t82 ^ 2;
t125 = -t81 * pkin(4) + t67 * qJ(5) + t82 * t154 - t73 * t55 + t147;
t129 = -t110 * t31 + t153 * t33;
t21 = -t67 * pkin(4) - t81 * qJ(5) + t74 * t55 + qJDD(5) - t129;
t16 = (-t39 - t150) * pkin(11) + (t73 * t74 - t67) * pkin(5) + t21;
t17 = -t72 * pkin(5) + t38 * pkin(11) + t82 * t63 + t125;
t36 = -t53 * mrSges(7,1) + t54 * mrSges(7,2);
t66 = qJDD(6) - t67;
t14 = m(7) * (-t109 * t17 + t114 * t16) - t27 * mrSges(7,3) + t66 * mrSges(7,1) - t54 * t36 + t80 * t42;
t15 = m(7) * (t109 * t16 + t114 * t17) + t26 * mrSges(7,3) - t66 * mrSges(7,2) + t53 * t36 - t80 * t43;
t126 = m(6) * t125 + t67 * mrSges(6,3) - t109 * t14 + t114 * t15 + t82 * t62;
t10 = m(5) * t147 - t67 * mrSges(5,2) + t144 * t73 + t148 * t38 - t82 * t61 + t126;
t122 = -m(6) * t21 - t109 * t15 - t114 * t14;
t11 = m(5) * t129 + (t60 + t59) * t82 + t144 * t74 + (mrSges(5,1) + mrSges(6,1)) * t67 + t148 * t39 + t122;
t76 = t101 * mrSges(4,1) - t85 * mrSges(4,3);
t9 = m(4) * t145 - t88 * mrSges(4,2) + t68 * mrSges(4,3) + t153 * t10 - t101 * t76 - t110 * t11 + t84 * t70;
t93 = (-mrSges(3,1) * t116 + mrSges(3,2) * t112) * t138;
t4 = m(3) * (-g(3) * t140 + t142) - t96 * mrSges(3,3) - t103 * mrSges(3,2) + t93 * t133 - t104 * t89 + t115 * t9 - t111 * t12;
t90 = -t104 * mrSges(3,2) + mrSges(3,3) * t133;
t6 = m(3) * (-t107 * t91 - t151) + t95 * mrSges(3,2) + t96 * mrSges(3,1) + t111 * t9 + t115 * t12 + (t112 * t89 - t116 * t90) * t138;
t119 = m(4) * t49 - t68 * mrSges(4,1) + t69 * mrSges(4,2) + t110 * t10 + t153 * t11 - t84 * t75 + t85 * t76;
t8 = m(3) * t127 + t103 * mrSges(3,1) - t95 * mrSges(3,3) + t104 * t90 - t93 * t134 - t119;
t136 = t108 * t6 + t8 * t139 + t4 * t140;
t2 = m(2) * t128 - t118 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t112 * t8 + t116 * t4;
t1 = m(2) * t131 + qJDD(1) * mrSges(2,1) - t118 * mrSges(2,2) - t107 * t6 + (t112 * t4 + t116 * t8) * t108;
t3 = [-m(1) * g(1) - t113 * t1 + t117 * t2, t2, t4, t9, t10, -t38 * mrSges(6,2) - t73 * t56 + t126, t15; -m(1) * g(2) + t117 * t1 + t113 * t2, t1, t8, t12, t11, -t39 * mrSges(6,3) - t74 * t62 + t123, t14; (-m(1) - m(2)) * g(3) + t136, -m(2) * g(3) + t136, t6, t119, t155, -t67 * mrSges(6,1) + t39 * mrSges(6,2) + t74 * t56 - t82 * t59 - t122, t130;];
f_new  = t3;
