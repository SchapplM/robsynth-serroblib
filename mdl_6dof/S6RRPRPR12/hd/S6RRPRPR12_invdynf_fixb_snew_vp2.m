% Calculate vector of cutting forces with Newton-Euler
% S6RRPRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-05-06 16:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRPR12_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR12_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR12_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR12_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR12_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR12_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR12_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 16:13:57
% EndTime: 2019-05-06 16:14:08
% DurationCPUTime: 4.70s
% Computational Cost: add. (58788->216), mult. (135127->285), div. (0->0), fcn. (100146->12), ass. (0->113)
t158 = -2 * qJD(3);
t110 = sin(qJ(2));
t105 = sin(pkin(6));
t141 = qJD(1) * t105;
t133 = t110 * t141;
t107 = cos(pkin(6));
t99 = qJD(1) * t107 + qJD(2);
t157 = (pkin(2) * t99 + t158) * t133;
t114 = cos(qJ(2));
t143 = t105 * t110;
t116 = qJD(1) ^ 2;
t111 = sin(qJ(1));
t115 = cos(qJ(1));
t132 = g(1) * t111 - g(2) * t115;
t82 = pkin(8) * t105 * t116 + qJDD(1) * pkin(1) + t132;
t147 = t107 * t82;
t130 = -g(1) * t115 - g(2) * t111;
t139 = qJDD(1) * t105;
t83 = -pkin(1) * t116 + pkin(8) * t139 + t130;
t127 = -g(3) * t143 + t110 * t147 + t114 * t83;
t140 = qJD(1) * t114;
t134 = t105 * t140;
t84 = (-pkin(2) * t114 - qJ(3) * t110) * t141;
t97 = t99 ^ 2;
t98 = qJDD(1) * t107 + qJDD(2);
t156 = t97 * pkin(2) - qJ(3) * t98 - t134 * t84 + t158 * t99 - t127;
t155 = 2 * qJD(5);
t154 = -pkin(2) - pkin(9);
t153 = t107 * g(3);
t152 = mrSges(3,1) - mrSges(4,2);
t151 = mrSges(3,3) + mrSges(4,1);
t109 = sin(qJ(4));
t113 = cos(qJ(4));
t144 = t105 ^ 2 * t116;
t135 = t114 ^ 2 * t144;
t87 = pkin(3) * t133 - pkin(9) * t99;
t88 = (qJD(2) * t140 + qJDD(1) * t110) * t105;
t89 = -qJD(2) * t133 + t114 * t139;
t32 = -pkin(3) * t135 - t153 - t88 * qJ(3) + t154 * t89 + (-t82 + (-qJ(3) * t114 * t99 - t110 * t87) * qJD(1)) * t105 + t157;
t142 = t105 * t114;
t146 = g(3) * t142 + t110 * t83;
t125 = -t97 * qJ(3) + t133 * t84 + qJDD(3) + t146;
t35 = t88 * pkin(3) + t154 * t98 + (-pkin(3) * t141 * t99 - pkin(9) * t110 * t144 - t147) * t114 + t125;
t150 = t109 * t35 + t113 * t32;
t81 = mrSges(4,1) * t133 + mrSges(4,2) * t99;
t149 = mrSges(3,1) * t99 - mrSges(3,3) * t133 - t81;
t85 = (mrSges(4,2) * t114 - mrSges(4,3) * t110) * t141;
t148 = t85 + (-mrSges(3,1) * t114 + mrSges(3,2) * t110) * t141;
t119 = t89 * pkin(3) - pkin(9) * t135 + t87 * t99 - t156;
t108 = sin(qJ(6));
t112 = cos(qJ(6));
t73 = -t109 * t134 + t113 * t99;
t56 = -qJD(4) * t73 - t109 * t98 - t113 * t89;
t93 = qJD(4) + t133;
t64 = pkin(4) * t93 - qJ(5) * t73;
t72 = -t109 * t99 - t113 * t134;
t71 = t72 ^ 2;
t117 = -t56 * pkin(4) - t71 * qJ(5) + t64 * t73 + qJDD(5) + t119;
t104 = sin(pkin(11));
t106 = cos(pkin(11));
t131 = -t109 * t32 + t113 * t35;
t57 = qJD(4) * t72 - t109 * t89 + t113 * t98;
t77 = qJDD(4) + t88;
t22 = (t72 * t93 - t57) * qJ(5) + (t72 * t73 + t77) * pkin(4) + t131;
t24 = -pkin(4) * t71 + qJ(5) * t56 - t64 * t93 + t150;
t60 = -t104 * t73 + t106 * t72;
t137 = t104 * t22 + t106 * t24 + t155 * t60;
t61 = t104 * t72 + t106 * t73;
t47 = -pkin(5) * t60 - pkin(10) * t61;
t91 = t93 ^ 2;
t19 = -pkin(5) * t91 + pkin(10) * t77 + t47 * t60 + t137;
t41 = -t104 * t57 + t106 * t56;
t42 = t104 * t56 + t106 * t57;
t20 = t117 + (-t60 * t93 - t42) * pkin(10) + (t61 * t93 - t41) * pkin(5);
t49 = -t108 * t61 + t112 * t93;
t28 = qJD(6) * t49 + t108 * t77 + t112 * t42;
t50 = t108 * t93 + t112 * t61;
t36 = -mrSges(7,1) * t49 + mrSges(7,2) * t50;
t59 = qJD(6) - t60;
t37 = -mrSges(7,2) * t59 + mrSges(7,3) * t49;
t40 = qJDD(6) - t41;
t16 = m(7) * (-t108 * t19 + t112 * t20) - t28 * mrSges(7,3) + t40 * mrSges(7,1) - t50 * t36 + t59 * t37;
t27 = -qJD(6) * t50 - t108 * t42 + t112 * t77;
t38 = mrSges(7,1) * t59 - mrSges(7,3) * t50;
t17 = m(7) * (t108 * t20 + t112 * t19) + t27 * mrSges(7,3) - t40 * mrSges(7,2) + t49 * t36 - t59 * t38;
t51 = -mrSges(6,2) * t93 + mrSges(6,3) * t60;
t52 = mrSges(6,1) * t93 - mrSges(6,3) * t61;
t122 = m(6) * t117 - t41 * mrSges(6,1) + mrSges(6,2) * t42 + t108 * t17 + t112 * t16 - t60 * t51 + t52 * t61;
t63 = -mrSges(5,2) * t93 + mrSges(5,3) * t72;
t65 = mrSges(5,1) * t93 - mrSges(5,3) * t73;
t120 = m(5) * t119 - t56 * mrSges(5,1) + mrSges(5,2) * t57 - t72 * t63 + t65 * t73 + t122;
t118 = -m(4) * t156 + t120;
t11 = t118 + t148 * t134 + m(3) * t127 - t149 * t99 + (-mrSges(3,2) + mrSges(4,3)) * t98 + t151 * t89;
t129 = -t105 * t82 - t153;
t46 = -mrSges(6,1) * t60 + mrSges(6,2) * t61;
t12 = m(6) * t137 - mrSges(6,2) * t77 + mrSges(6,3) * t41 - t108 * t16 + t112 * t17 + t46 * t60 - t52 * t93;
t128 = t104 * t24 - t106 * t22;
t121 = m(7) * (-t77 * pkin(5) - t91 * pkin(10) + (t155 + t47) * t61 + t128) - t27 * mrSges(7,1) + t28 * mrSges(7,2) - t49 * t37 + t50 * t38;
t13 = m(6) * (-0.2e1 * qJD(5) * t61 - t128) - t42 * mrSges(6,3) + t77 * mrSges(6,1) - t61 * t46 + t93 * t51 - t121;
t62 = -mrSges(5,1) * t72 + mrSges(5,2) * t73;
t8 = m(5) * t131 + mrSges(5,1) * t77 - mrSges(5,3) * t57 + t104 * t12 + t106 * t13 - t62 * t73 + t63 * t93;
t80 = -mrSges(4,1) * t134 - mrSges(4,3) * t99;
t9 = m(5) * t150 - mrSges(5,2) * t77 + mrSges(5,3) * t56 - t104 * t13 + t106 * t12 + t62 * t72 - t65 * t93;
t126 = -t109 * t8 + t113 * t9 + m(4) * (-t89 * pkin(2) + (-t134 * t99 - t88) * qJ(3) + t129 + t157) + t80 * t134 - t88 * mrSges(4,3);
t79 = -mrSges(3,2) * t99 + mrSges(3,3) * t134;
t5 = m(3) * t129 + t88 * mrSges(3,2) - t152 * t89 + (t110 * t149 - t114 * t79) * t141 + t126;
t136 = t114 * t147;
t124 = -m(4) * (-t98 * pkin(2) + t125 - t136) - t109 * t9 - t113 * t8;
t6 = m(3) * (t136 - t146) + (t79 - t80) * t99 + t152 * t98 - t151 * t88 - t148 * t133 + t124;
t138 = t107 * t5 + t11 * t143 + t142 * t6;
t2 = m(2) * t130 - mrSges(2,1) * t116 - qJDD(1) * mrSges(2,2) + t11 * t114 - t110 * t6;
t1 = m(2) * t132 + qJDD(1) * mrSges(2,1) - t116 * mrSges(2,2) - t105 * t5 + (t11 * t110 + t114 * t6) * t107;
t3 = [-m(1) * g(1) - t1 * t111 + t115 * t2, t2, t11, t89 * mrSges(4,2) - t133 * t81 + t126, t9, t12, t17; -m(1) * g(2) + t1 * t115 + t111 * t2, t1, t6, -t89 * mrSges(4,1) - t98 * mrSges(4,3) - t134 * t85 - t99 * t81 - t118, t8, t13, t16; (-m(1) - m(2)) * g(3) + t138, -m(2) * g(3) + t138, t5, t88 * mrSges(4,1) + t98 * mrSges(4,2) + t133 * t85 + t99 * t80 - t124, t120, t122, t121;];
f_new  = t3;
