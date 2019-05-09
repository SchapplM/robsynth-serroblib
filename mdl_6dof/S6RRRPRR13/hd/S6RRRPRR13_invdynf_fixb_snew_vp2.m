% Calculate vector of cutting forces with Newton-Euler
% S6RRRPRR13
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-05-07 15:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPRR13_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR13_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR13_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR13_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_invdynf_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR13_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR13_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR13_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 15:32:40
% EndTime: 2019-05-07 15:33:23
% DurationCPUTime: 21.56s
% Computational Cost: add. (397242->229), mult. (994588->328), div. (0->0), fcn. (842384->16), ass. (0->130)
t110 = sin(pkin(7));
t113 = cos(pkin(7));
t117 = sin(qJ(3));
t122 = cos(qJ(3));
t111 = sin(pkin(6));
t118 = sin(qJ(2));
t123 = cos(qJ(2));
t145 = qJD(1) * qJD(2);
t101 = (qJDD(1) * t118 + t123 * t145) * t111;
t114 = cos(pkin(6));
t106 = t114 * qJDD(1) + qJDD(2);
t107 = t114 * qJD(1) + qJD(2);
t125 = qJD(1) ^ 2;
t119 = sin(qJ(1));
t124 = cos(qJ(1));
t139 = t119 * g(1) - t124 * g(2);
t162 = pkin(9) * t111;
t98 = qJDD(1) * pkin(1) + t125 * t162 + t139;
t156 = t114 * t98;
t134 = -t124 * g(1) - t119 * g(2);
t99 = -t125 * pkin(1) + qJDD(1) * t162 + t134;
t137 = -t118 * t99 + t123 * t156;
t148 = qJD(1) * t118;
t160 = pkin(10) * t113;
t147 = qJD(1) * t123;
t141 = t111 * t147;
t136 = t113 * t141;
t90 = (t107 * t110 + t136) * pkin(10);
t149 = qJD(1) * t111;
t161 = pkin(10) * t110;
t95 = (-pkin(2) * t123 - t118 * t161) * t149;
t58 = -t101 * t160 + t106 * pkin(2) + t107 * t90 + (-g(3) * t123 - t95 * t148) * t111 + t137;
t102 = (qJDD(1) * t123 - t118 * t145) * t111;
t155 = t102 * t113;
t130 = t106 * t110 + t155;
t157 = t118 * t156 + t123 * t99;
t142 = t111 * t148;
t94 = t107 * pkin(2) - t142 * t160;
t59 = -t107 * t94 + (-g(3) * t118 + t95 * t147) * t111 + t130 * pkin(10) + t157;
t159 = t114 * g(3);
t65 = -t101 * t161 - t102 * pkin(2) - t159 + (-t98 + (t118 * t94 - t123 * t90) * qJD(1)) * t111;
t163 = -t117 * t59 + (t110 * t65 + t113 * t58) * t122;
t116 = sin(qJ(5));
t121 = cos(qJ(5));
t109 = sin(pkin(13));
t112 = cos(pkin(13));
t150 = t113 * t117;
t154 = t110 * t117;
t143 = t122 * t59 + t58 * t150 + t65 * t154;
t153 = t110 * t122;
t84 = -t107 * t153 + t117 * t142 - t122 * t136;
t85 = t107 * t154 + (t118 * t122 + t123 * t150) * t149;
t73 = t84 * pkin(3) - t85 * qJ(4);
t86 = -t110 * t102 + t113 * t106 + qJDD(3);
t92 = t113 * t107 - t110 * t141 + qJD(3);
t89 = t92 ^ 2;
t33 = -t89 * pkin(3) + t86 * qJ(4) - t84 * t73 + t143;
t138 = -t110 * t58 + t113 * t65;
t71 = t85 * qJD(3) + t117 * t101 - t106 * t153 - t122 * t155;
t72 = -t84 * qJD(3) + t122 * t101 + t130 * t117;
t36 = (t84 * t92 - t72) * qJ(4) + (t85 * t92 + t71) * pkin(3) + t138;
t79 = t109 * t92 + t112 * t85;
t135 = -0.2e1 * qJD(4) * t79 - t109 * t33 + t112 * t36;
t62 = t109 * t86 + t112 * t72;
t78 = -t109 * t85 + t112 * t92;
t24 = (t78 * t84 - t62) * pkin(11) + (t78 * t79 + t71) * pkin(4) + t135;
t144 = 0.2e1 * qJD(4) * t78 + t109 * t36 + t112 * t33;
t61 = -t109 * t72 + t112 * t86;
t69 = t84 * pkin(4) - t79 * pkin(11);
t77 = t78 ^ 2;
t26 = -t77 * pkin(4) + t61 * pkin(11) - t84 * t69 + t144;
t158 = t116 * t24 + t121 * t26;
t152 = t111 * t118;
t151 = t111 * t123;
t115 = sin(qJ(6));
t120 = cos(qJ(6));
t56 = -t116 * t79 + t121 * t78;
t57 = t116 * t78 + t121 * t79;
t46 = -t56 * pkin(5) - t57 * pkin(12);
t70 = qJDD(5) + t71;
t83 = qJD(5) + t84;
t82 = t83 ^ 2;
t21 = -t82 * pkin(5) + t70 * pkin(12) + t56 * t46 + t158;
t32 = -t86 * pkin(3) - t89 * qJ(4) + t85 * t73 + qJDD(4) - t163;
t127 = -t61 * pkin(4) - t77 * pkin(11) + t79 * t69 + t32;
t39 = -t57 * qJD(5) - t116 * t62 + t121 * t61;
t40 = t56 * qJD(5) + t116 * t61 + t121 * t62;
t22 = (-t56 * t83 - t40) * pkin(12) + (t57 * t83 - t39) * pkin(5) + t127;
t47 = -t115 * t57 + t120 * t83;
t30 = t47 * qJD(6) + t115 * t70 + t120 * t40;
t38 = qJDD(6) - t39;
t48 = t115 * t83 + t120 * t57;
t41 = -t47 * mrSges(7,1) + t48 * mrSges(7,2);
t55 = qJD(6) - t56;
t42 = -t55 * mrSges(7,2) + t47 * mrSges(7,3);
t18 = m(7) * (-t115 * t21 + t120 * t22) - t30 * mrSges(7,3) + t38 * mrSges(7,1) - t48 * t41 + t55 * t42;
t29 = -t48 * qJD(6) - t115 * t40 + t120 * t70;
t43 = t55 * mrSges(7,1) - t48 * mrSges(7,3);
t19 = m(7) * (t115 * t22 + t120 * t21) + t29 * mrSges(7,3) - t38 * mrSges(7,2) + t47 * t41 - t55 * t43;
t45 = -t56 * mrSges(6,1) + t57 * mrSges(6,2);
t50 = t83 * mrSges(6,1) - t57 * mrSges(6,3);
t13 = m(6) * t158 - t70 * mrSges(6,2) + t39 * mrSges(6,3) - t115 * t18 + t120 * t19 + t56 * t45 - t83 * t50;
t131 = -t116 * t26 + t121 * t24;
t128 = m(7) * (-t70 * pkin(5) - t82 * pkin(12) + t57 * t46 - t131) - t29 * mrSges(7,1) + t30 * mrSges(7,2) - t47 * t42 + t48 * t43;
t49 = -t83 * mrSges(6,2) + t56 * mrSges(6,3);
t15 = m(6) * t131 + t70 * mrSges(6,1) - t40 * mrSges(6,3) - t57 * t45 + t83 * t49 - t128;
t60 = -t78 * mrSges(5,1) + t79 * mrSges(5,2);
t67 = -t84 * mrSges(5,2) + t78 * mrSges(5,3);
t11 = m(5) * t135 + t71 * mrSges(5,1) - t62 * mrSges(5,3) + t116 * t13 + t121 * t15 - t79 * t60 + t84 * t67;
t68 = t84 * mrSges(5,1) - t79 * mrSges(5,3);
t12 = m(5) * t144 - t71 * mrSges(5,2) + t61 * mrSges(5,3) - t116 * t15 + t121 * t13 + t78 * t60 - t84 * t68;
t80 = -t92 * mrSges(4,2) - t84 * mrSges(4,3);
t81 = t92 * mrSges(4,1) - t85 * mrSges(4,3);
t10 = m(4) * t138 + t71 * mrSges(4,1) + t72 * mrSges(4,2) + t109 * t12 + t112 * t11 + t84 * t80 + t85 * t81;
t129 = -m(6) * t127 + t39 * mrSges(6,1) - t40 * mrSges(6,2) - t115 * t19 - t120 * t18 + t56 * t49 - t57 * t50;
t126 = m(5) * t32 - t61 * mrSges(5,1) + t62 * mrSges(5,2) - t78 * t67 + t79 * t68 - t129;
t74 = t84 * mrSges(4,1) + t85 * mrSges(4,2);
t14 = m(4) * t163 + t86 * mrSges(4,1) - t72 * mrSges(4,3) - t85 * t74 + t92 * t80 - t126;
t9 = m(4) * t143 - t86 * mrSges(4,2) - t71 * mrSges(4,3) - t109 * t11 + t112 * t12 - t84 * t74 - t92 * t81;
t133 = t117 * t9 + t122 * t14;
t140 = (-mrSges(3,1) * t123 + mrSges(3,2) * t118) * t149 ^ 2;
t97 = -t107 * mrSges(3,2) + mrSges(3,3) * t141;
t4 = m(3) * (-g(3) * t151 + t137) - t101 * mrSges(3,3) + t106 * mrSges(3,1) - t118 * t140 + t107 * t97 - t110 * t10 + t133 * t113;
t96 = t107 * mrSges(3,1) - mrSges(3,3) * t142;
t6 = m(3) * (-t111 * t98 - t159) + t101 * mrSges(3,2) - t102 * mrSges(3,1) + t113 * t10 + t133 * t110 + (t118 * t96 - t123 * t97) * t149;
t8 = m(3) * (-g(3) * t152 + t157) + t102 * mrSges(3,3) - t106 * mrSges(3,2) + t123 * t140 - t107 * t96 + t122 * t9 - t117 * t14;
t146 = t114 * t6 + t4 * t151 + t8 * t152;
t2 = m(2) * t134 - t125 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t118 * t4 + t123 * t8;
t1 = m(2) * t139 + qJDD(1) * mrSges(2,1) - t125 * mrSges(2,2) - t111 * t6 + (t118 * t8 + t123 * t4) * t114;
t3 = [-m(1) * g(1) - t119 * t1 + t124 * t2, t2, t8, t9, t12, t13, t19; -m(1) * g(2) + t124 * t1 + t119 * t2, t1, t4, t14, t11, t15, t18; (-m(1) - m(2)) * g(3) + t146, -m(2) * g(3) + t146, t6, t10, t126, -t129, t128;];
f_new  = t3;
