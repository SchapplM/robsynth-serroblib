% Calculate vector of cutting forces with Newton-Euler
% S6RRRRRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 14:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRRR8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR8_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR8_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_invdynf_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR8_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR8_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR8_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 13:32:18
% EndTime: 2019-05-08 13:33:53
% DurationCPUTime: 23.28s
% Computational Cost: add. (426372->228), mult. (1045830->323), div. (0->0), fcn. (891769->16), ass. (0->130)
t108 = sin(pkin(7));
t110 = cos(pkin(7));
t115 = sin(qJ(3));
t121 = cos(qJ(3));
t109 = sin(pkin(6));
t116 = sin(qJ(2));
t122 = cos(qJ(2));
t143 = qJD(1) * qJD(2);
t102 = (qJDD(1) * t116 + t122 * t143) * t109;
t111 = cos(pkin(6));
t105 = t111 * qJDD(1) + qJDD(2);
t106 = t111 * qJD(1) + qJD(2);
t124 = qJD(1) ^ 2;
t117 = sin(qJ(1));
t123 = cos(qJ(1));
t134 = -t123 * g(1) - t117 * g(2);
t160 = pkin(9) * t109;
t100 = -t124 * pkin(1) + qJDD(1) * t160 + t134;
t138 = t117 * g(1) - t123 * g(2);
t99 = qJDD(1) * pkin(1) + t124 * t160 + t138;
t152 = t111 * t99;
t135 = -t116 * t100 + t122 * t152;
t146 = qJD(1) * t116;
t158 = pkin(10) * t110;
t145 = qJD(1) * t122;
t140 = t109 * t145;
t129 = t106 * t108 + t110 * t140;
t91 = t129 * pkin(10);
t147 = qJD(1) * t109;
t159 = pkin(10) * t108;
t96 = (-pkin(2) * t122 - t116 * t159) * t147;
t58 = -t102 * t158 + t105 * pkin(2) + t106 * t91 + (-g(3) * t122 - t96 * t146) * t109 + t135;
t103 = (qJDD(1) * t122 - t116 * t143) * t109;
t130 = t103 * t110 + t105 * t108;
t153 = t122 * t100 + t116 * t152;
t141 = t109 * t146;
t94 = t106 * pkin(2) - t141 * t158;
t59 = -t106 * t94 + (-g(3) * t116 + t96 * t145) * t109 + t130 * pkin(10) + t153;
t157 = t111 * g(3);
t65 = -t102 * t159 - t103 * pkin(2) - t157 + (-t99 + (t116 * t94 - t122 * t91) * qJD(1)) * t109;
t161 = -t115 * t59 + (t108 * t65 + t110 * t58) * t121;
t85 = -t115 * t141 + t129 * t121;
t148 = t110 * t115;
t151 = t108 * t115;
t86 = t106 * t151 + (t116 * t121 + t122 * t148) * t147;
t72 = -t86 * qJD(3) - t115 * t102 + t130 * t121;
t113 = sin(qJ(5));
t119 = cos(qJ(5));
t114 = sin(qJ(4));
t120 = cos(qJ(4));
t142 = t121 * t59 + t58 * t148 + t65 * t151;
t75 = -t85 * pkin(3) - t86 * pkin(11);
t87 = -t108 * t103 + t110 * t105 + qJDD(3);
t92 = t110 * t106 - t108 * t140 + qJD(3);
t90 = t92 ^ 2;
t33 = -t90 * pkin(3) + t87 * pkin(11) + t85 * t75 + t142;
t137 = -t108 * t58 + t110 * t65;
t73 = t85 * qJD(3) + t121 * t102 + t130 * t115;
t36 = (-t85 * t92 - t73) * pkin(11) + (t86 * t92 - t72) * pkin(3) + t137;
t136 = -t114 * t33 + t120 * t36;
t77 = -t114 * t86 + t120 * t92;
t48 = t77 * qJD(4) + t114 * t87 + t120 * t73;
t71 = qJDD(4) - t72;
t78 = t114 * t92 + t120 * t86;
t84 = qJD(4) - t85;
t24 = (t77 * t84 - t48) * pkin(12) + (t77 * t78 + t71) * pkin(4) + t136;
t155 = t114 * t36 + t120 * t33;
t47 = -t78 * qJD(4) - t114 * t73 + t120 * t87;
t69 = t84 * pkin(4) - t78 * pkin(12);
t76 = t77 ^ 2;
t26 = -t76 * pkin(4) + t47 * pkin(12) - t84 * t69 + t155;
t156 = t113 * t24 + t119 * t26;
t150 = t109 * t116;
t149 = t109 * t122;
t112 = sin(qJ(6));
t118 = cos(qJ(6));
t60 = -t113 * t78 + t119 * t77;
t61 = t113 * t77 + t119 * t78;
t46 = -t60 * pkin(5) - t61 * pkin(13);
t70 = qJDD(5) + t71;
t82 = qJD(5) + t84;
t81 = t82 ^ 2;
t21 = -t81 * pkin(5) + t70 * pkin(13) + t60 * t46 + t156;
t32 = -t87 * pkin(3) - t90 * pkin(11) + t86 * t75 - t161;
t126 = -t47 * pkin(4) - t76 * pkin(12) + t78 * t69 + t32;
t39 = -t61 * qJD(5) - t113 * t48 + t119 * t47;
t40 = t60 * qJD(5) + t113 * t47 + t119 * t48;
t22 = (-t60 * t82 - t40) * pkin(13) + (t61 * t82 - t39) * pkin(5) + t126;
t49 = -t112 * t61 + t118 * t82;
t30 = t49 * qJD(6) + t112 * t70 + t118 * t40;
t38 = qJDD(6) - t39;
t50 = t112 * t82 + t118 * t61;
t41 = -t49 * mrSges(7,1) + t50 * mrSges(7,2);
t57 = qJD(6) - t60;
t42 = -t57 * mrSges(7,2) + t49 * mrSges(7,3);
t18 = m(7) * (-t112 * t21 + t118 * t22) - t30 * mrSges(7,3) + t38 * mrSges(7,1) - t50 * t41 + t57 * t42;
t29 = -t50 * qJD(6) - t112 * t40 + t118 * t70;
t43 = t57 * mrSges(7,1) - t50 * mrSges(7,3);
t19 = m(7) * (t112 * t22 + t118 * t21) + t29 * mrSges(7,3) - t38 * mrSges(7,2) + t49 * t41 - t57 * t43;
t45 = -t60 * mrSges(6,1) + t61 * mrSges(6,2);
t52 = t82 * mrSges(6,1) - t61 * mrSges(6,3);
t14 = m(6) * t156 - t70 * mrSges(6,2) + t39 * mrSges(6,3) - t112 * t18 + t118 * t19 + t60 * t45 - t82 * t52;
t131 = -t113 * t26 + t119 * t24;
t127 = m(7) * (-t70 * pkin(5) - t81 * pkin(13) + t61 * t46 - t131) - t29 * mrSges(7,1) + t30 * mrSges(7,2) - t49 * t42 + t50 * t43;
t51 = -t82 * mrSges(6,2) + t60 * mrSges(6,3);
t15 = m(6) * t131 + t70 * mrSges(6,1) - t40 * mrSges(6,3) - t61 * t45 + t82 * t51 - t127;
t62 = -t77 * mrSges(5,1) + t78 * mrSges(5,2);
t67 = -t84 * mrSges(5,2) + t77 * mrSges(5,3);
t11 = m(5) * t136 + t71 * mrSges(5,1) - t48 * mrSges(5,3) + t113 * t14 + t119 * t15 - t78 * t62 + t84 * t67;
t68 = t84 * mrSges(5,1) - t78 * mrSges(5,3);
t12 = m(5) * t155 - t71 * mrSges(5,2) + t47 * mrSges(5,3) - t113 * t15 + t119 * t14 + t77 * t62 - t84 * t68;
t79 = -t92 * mrSges(4,2) + t85 * mrSges(4,3);
t80 = t92 * mrSges(4,1) - t86 * mrSges(4,3);
t10 = m(4) * t137 - t72 * mrSges(4,1) + t73 * mrSges(4,2) + t120 * t11 + t114 * t12 - t85 * t79 + t86 * t80;
t128 = -m(6) * t126 + t39 * mrSges(6,1) - t40 * mrSges(6,2) - t112 * t19 - t118 * t18 + t60 * t51 - t61 * t52;
t125 = m(5) * t32 - t47 * mrSges(5,1) + t48 * mrSges(5,2) - t77 * t67 + t78 * t68 - t128;
t74 = -t85 * mrSges(4,1) + t86 * mrSges(4,2);
t13 = m(4) * t161 + t87 * mrSges(4,1) - t73 * mrSges(4,3) - t86 * t74 + t92 * t79 - t125;
t9 = m(4) * t142 - t87 * mrSges(4,2) + t72 * mrSges(4,3) - t114 * t11 + t120 * t12 + t85 * t74 - t92 * t80;
t133 = t115 * t9 + t121 * t13;
t139 = (-mrSges(3,1) * t122 + mrSges(3,2) * t116) * t147 ^ 2;
t98 = -t106 * mrSges(3,2) + mrSges(3,3) * t140;
t4 = m(3) * (-g(3) * t149 + t135) - t102 * mrSges(3,3) + t105 * mrSges(3,1) - t116 * t139 + t106 * t98 - t108 * t10 + t133 * t110;
t97 = t106 * mrSges(3,1) - mrSges(3,3) * t141;
t6 = m(3) * (-t109 * t99 - t157) + t102 * mrSges(3,2) - t103 * mrSges(3,1) + t110 * t10 + t133 * t108 + (t116 * t97 - t122 * t98) * t147;
t8 = m(3) * (-g(3) * t150 + t153) + t103 * mrSges(3,3) - t105 * mrSges(3,2) + t122 * t139 - t106 * t97 + t121 * t9 - t115 * t13;
t144 = t111 * t6 + t4 * t149 + t8 * t150;
t2 = m(2) * t134 - t124 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t116 * t4 + t122 * t8;
t1 = m(2) * t138 + qJDD(1) * mrSges(2,1) - t124 * mrSges(2,2) - t109 * t6 + (t116 * t8 + t122 * t4) * t111;
t3 = [-m(1) * g(1) - t117 * t1 + t123 * t2, t2, t8, t9, t12, t14, t19; -m(1) * g(2) + t123 * t1 + t117 * t2, t1, t4, t13, t11, t15, t18; (-m(1) - m(2)) * g(3) + t144, -m(2) * g(3) + t144, t6, t10, t125, -t128, t127;];
f_new  = t3;
