% Calculate vector of cutting forces with Newton-Euler
% S6RRPRRR14
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-09 12:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRRR14_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR14_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_invdynf_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR14_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR14_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-09 11:05:33
% EndTime: 2019-05-09 11:07:26
% DurationCPUTime: 51.02s
% Computational Cost: add. (881902->237), mult. (2391917->349), div. (0->0), fcn. (2062186->18), ass. (0->141)
t107 = sin(pkin(8));
t111 = cos(pkin(8));
t106 = sin(pkin(14));
t110 = cos(pkin(14));
t108 = sin(pkin(7));
t154 = t108 * t110;
t112 = cos(pkin(7));
t109 = sin(pkin(6));
t117 = sin(qJ(2));
t122 = cos(qJ(2));
t146 = qJD(1) * qJD(2);
t101 = (qJDD(1) * t117 + t122 * t146) * t109;
t113 = cos(pkin(6));
t103 = t113 * qJDD(1) + qJDD(2);
t104 = t113 * qJD(1) + qJD(2);
t124 = qJD(1) ^ 2;
t118 = sin(qJ(1));
t123 = cos(qJ(1));
t141 = t118 * g(1) - t123 * g(2);
t168 = pkin(10) * t109;
t98 = qJDD(1) * pkin(1) + t124 * t168 + t141;
t159 = t113 * t98;
t137 = -t123 * g(1) - t118 * g(2);
t99 = -t124 * pkin(1) + qJDD(1) * t168 + t137;
t139 = -t117 * t99 + t122 * t159;
t149 = qJD(1) * t117;
t156 = qJ(3) * t112;
t148 = qJD(1) * t122;
t143 = t109 * t148;
t89 = (t104 * t108 + t112 * t143) * qJ(3);
t150 = qJD(1) * t109;
t157 = qJ(3) * t108;
t95 = (-pkin(2) * t122 - t117 * t157) * t150;
t63 = -t101 * t156 + t103 * pkin(2) + t104 * t89 + (-g(3) * t122 - t95 * t149) * t109 + t139;
t160 = t112 * t63;
t102 = (qJDD(1) * t122 - t117 * t146) * t109;
t127 = t102 * t112 + t103 * t108;
t163 = t117 * t159 + t122 * t99;
t144 = t109 * t149;
t94 = t104 * pkin(2) - t144 * t156;
t64 = -t104 * t94 + (-g(3) * t117 + t95 * t148) * t109 + t127 * qJ(3) + t163;
t165 = t113 * g(3);
t68 = -t101 * t157 - t102 * pkin(2) - t165 + (-t98 + (t117 * t94 - t122 * t89) * qJD(1)) * t109;
t151 = t112 * t122;
t155 = t106 * t108;
t87 = t104 * t155 + (t106 * t151 + t110 * t117) * t150;
t128 = -0.2e1 * qJD(3) * t87 - t106 * t64 + t110 * t160 + t68 * t154;
t166 = pkin(11) * t111;
t167 = pkin(11) * t107;
t86 = t104 * t154 + (-t106 * t117 + t110 * t151) * t150;
t74 = -t86 * pkin(3) - t87 * t167;
t92 = t112 * t104 - t108 * t143;
t132 = t107 * t92 + t111 * t86;
t77 = t132 * pkin(11);
t83 = t110 * t101 + t127 * t106;
t88 = -t108 * t102 + t112 * t103;
t32 = t88 * pkin(3) - t83 * t166 - t87 * t74 + t92 * t77 + t128;
t136 = -t108 * t63 + t112 * t68 + qJDD(3);
t79 = t92 * pkin(3) - t87 * t166;
t82 = -t106 * t101 + t127 * t110;
t37 = -t82 * pkin(3) - t83 * t167 - t86 * t77 + t87 * t79 + t136;
t170 = t107 * t37 + t111 * t32;
t116 = sin(qJ(4));
t121 = cos(qJ(4));
t133 = t107 * t88 + t111 * t82;
t138 = 0.2e1 * qJD(3) * t86 + t106 * t160 + t110 * t64 + t68 * t155;
t33 = t133 * pkin(11) + t86 * t74 - t92 * t79 + t138;
t169 = -t116 * t33 + t170 * t121;
t72 = t132 * t116 + t121 * t87;
t51 = -t72 * qJD(4) - t116 * t83 + t133 * t121;
t71 = -t116 * t87 + t132 * t121;
t115 = sin(qJ(5));
t120 = cos(qJ(5));
t145 = t170 * t116 + t121 * t33;
t54 = -t71 * pkin(4) - t72 * pkin(12);
t73 = -t107 * t82 + t111 * t88 + qJDD(4);
t78 = -t107 * t86 + t111 * t92 + qJD(4);
t76 = t78 ^ 2;
t24 = -t76 * pkin(4) + t73 * pkin(12) + t71 * t54 + t145;
t140 = -t107 * t32 + t111 * t37;
t52 = t71 * qJD(4) + t133 * t116 + t121 * t83;
t26 = (-t71 * t78 - t52) * pkin(12) + (t72 * t78 - t51) * pkin(4) + t140;
t164 = t115 * t26 + t120 * t24;
t153 = t109 * t117;
t152 = t109 * t122;
t114 = sin(qJ(6));
t119 = cos(qJ(6));
t56 = -t115 * t72 + t120 * t78;
t57 = t115 * t78 + t120 * t72;
t44 = -t56 * pkin(5) - t57 * pkin(13);
t50 = qJDD(5) - t51;
t70 = qJD(5) - t71;
t69 = t70 ^ 2;
t20 = -t69 * pkin(5) + t50 * pkin(13) + t56 * t44 + t164;
t23 = -t73 * pkin(4) - t76 * pkin(12) + t72 * t54 - t169;
t39 = -t57 * qJD(5) - t115 * t52 + t120 * t73;
t40 = t56 * qJD(5) + t115 * t73 + t120 * t52;
t21 = (-t56 * t70 - t40) * pkin(13) + (t57 * t70 - t39) * pkin(5) + t23;
t46 = -t114 * t57 + t119 * t70;
t28 = t46 * qJD(6) + t114 * t50 + t119 * t40;
t47 = t114 * t70 + t119 * t57;
t34 = -t46 * mrSges(7,1) + t47 * mrSges(7,2);
t38 = qJDD(6) - t39;
t55 = qJD(6) - t56;
t41 = -t55 * mrSges(7,2) + t46 * mrSges(7,3);
t17 = m(7) * (-t114 * t20 + t119 * t21) - t28 * mrSges(7,3) + t38 * mrSges(7,1) - t47 * t34 + t55 * t41;
t27 = -t47 * qJD(6) - t114 * t40 + t119 * t50;
t42 = t55 * mrSges(7,1) - t47 * mrSges(7,3);
t18 = m(7) * (t114 * t21 + t119 * t20) + t27 * mrSges(7,3) - t38 * mrSges(7,2) + t46 * t34 - t55 * t42;
t43 = -t56 * mrSges(6,1) + t57 * mrSges(6,2);
t49 = t70 * mrSges(6,1) - t57 * mrSges(6,3);
t15 = m(6) * t164 - t50 * mrSges(6,2) + t39 * mrSges(6,3) - t114 * t17 + t119 * t18 + t56 * t43 - t70 * t49;
t131 = -t115 * t24 + t120 * t26;
t126 = m(7) * (-t50 * pkin(5) - t69 * pkin(13) + t57 * t44 - t131) - t27 * mrSges(7,1) + t28 * mrSges(7,2) - t46 * t41 + t47 * t42;
t48 = -t70 * mrSges(6,2) + t56 * mrSges(6,3);
t16 = m(6) * t131 + t50 * mrSges(6,1) - t40 * mrSges(6,3) - t57 * t43 + t70 * t48 - t126;
t58 = -t78 * mrSges(5,2) + t71 * mrSges(5,3);
t59 = t78 * mrSges(5,1) - t72 * mrSges(5,3);
t13 = m(5) * t140 - t51 * mrSges(5,1) + t52 * mrSges(5,2) + t115 * t15 + t120 * t16 - t71 * t58 + t72 * t59;
t53 = -t71 * mrSges(5,1) + t72 * mrSges(5,2);
t12 = m(5) * t145 - t73 * mrSges(5,2) + t51 * mrSges(5,3) - t115 * t16 + t120 * t15 + t71 * t53 - t78 * t59;
t125 = m(6) * t23 - t39 * mrSges(6,1) + t40 * mrSges(6,2) + t114 * t18 + t119 * t17 - t56 * t48 + t57 * t49;
t14 = m(5) * t169 + t73 * mrSges(5,1) - t52 * mrSges(5,3) - t72 * t53 + t78 * t58 - t125;
t130 = t116 * t12 + t121 * t14;
t80 = -t92 * mrSges(4,2) + t86 * mrSges(4,3);
t81 = t92 * mrSges(4,1) - t87 * mrSges(4,3);
t10 = m(4) * t136 - t82 * mrSges(4,1) + t83 * mrSges(4,2) + t130 * t107 + t111 * t13 - t86 * t80 + t87 * t81;
t75 = -t86 * mrSges(4,1) + t87 * mrSges(4,2);
t11 = m(4) * t138 - t88 * mrSges(4,2) + t82 * mrSges(4,3) - t116 * t14 + t121 * t12 + t86 * t75 - t92 * t81;
t9 = m(4) * t128 + t88 * mrSges(4,1) - t83 * mrSges(4,3) - t107 * t13 + t130 * t111 - t87 * t75 + t92 * t80;
t135 = t106 * t11 + t110 * t9;
t142 = (-mrSges(3,1) * t122 + mrSges(3,2) * t117) * t150 ^ 2;
t97 = -t104 * mrSges(3,2) + mrSges(3,3) * t143;
t4 = m(3) * (-g(3) * t152 + t139) - t101 * mrSges(3,3) + t103 * mrSges(3,1) - t117 * t142 + t104 * t97 - t108 * t10 + t135 * t112;
t96 = t104 * mrSges(3,1) - mrSges(3,3) * t144;
t6 = m(3) * (-t109 * t98 - t165) + t101 * mrSges(3,2) - t102 * mrSges(3,1) + t112 * t10 + t135 * t108 + (t117 * t96 - t122 * t97) * t150;
t8 = m(3) * (-g(3) * t153 + t163) + t102 * mrSges(3,3) - t103 * mrSges(3,2) + t122 * t142 - t104 * t96 + t110 * t11 - t106 * t9;
t147 = t113 * t6 + t4 * t152 + t8 * t153;
t2 = m(2) * t137 - t124 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t117 * t4 + t122 * t8;
t1 = m(2) * t141 + qJDD(1) * mrSges(2,1) - t124 * mrSges(2,2) - t109 * t6 + (t117 * t8 + t122 * t4) * t113;
t3 = [-m(1) * g(1) - t118 * t1 + t123 * t2, t2, t8, t11, t12, t15, t18; -m(1) * g(2) + t123 * t1 + t118 * t2, t1, t4, t9, t14, t16, t17; (-m(1) - m(2)) * g(3) + t147, -m(2) * g(3) + t147, t6, t10, t13, t125, t126;];
f_new  = t3;
