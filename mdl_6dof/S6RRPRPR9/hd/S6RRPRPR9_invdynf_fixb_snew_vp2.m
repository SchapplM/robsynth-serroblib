% Calculate vector of cutting forces with Newton-Euler
% S6RRPRPR9
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-05-06 15:23
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRPR9_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR9_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR9_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR9_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR9_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR9_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR9_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 15:12:19
% EndTime: 2019-05-06 15:12:40
% DurationCPUTime: 9.26s
% Computational Cost: add. (166135->213), mult. (380525->297), div. (0->0), fcn. (309261->14), ass. (0->113)
t149 = cos(qJ(4));
t115 = cos(pkin(6));
t148 = t115 * g(3);
t117 = sin(qJ(4));
t121 = cos(qJ(2));
t118 = sin(qJ(2));
t112 = sin(pkin(6));
t142 = qJD(1) * t112;
t136 = t118 * t142;
t139 = qJDD(1) * t112;
t101 = -qJD(2) * t136 + t121 * t139;
t111 = sin(pkin(11));
t114 = cos(pkin(11));
t107 = qJD(1) * t115 + qJD(2);
t105 = t107 ^ 2;
t106 = qJDD(1) * t115 + qJDD(2);
t141 = qJD(1) * t121;
t123 = qJD(1) ^ 2;
t119 = sin(qJ(1));
t122 = cos(qJ(1));
t134 = g(1) * t119 - g(2) * t122;
t96 = pkin(8) * t112 * t123 + qJDD(1) * pkin(1) + t134;
t145 = t115 * t96;
t130 = -g(1) * t122 - g(2) * t119;
t97 = -pkin(1) * t123 + pkin(8) * t139 + t130;
t146 = t118 * t145 + t121 * t97;
t98 = (-pkin(2) * t121 - qJ(3) * t118) * t142;
t63 = -t105 * pkin(2) + t106 * qJ(3) + (-g(3) * t118 + t141 * t98) * t112 + t146;
t100 = (qJD(2) * t141 + qJDD(1) * t118) * t112;
t64 = -t101 * pkin(2) - t148 - t100 * qJ(3) + (-t96 + (pkin(2) * t118 - qJ(3) * t121) * t107 * qJD(1)) * t112;
t90 = t107 * t111 + t114 * t136;
t132 = -0.2e1 * qJD(3) * t90 - t111 * t63 + t114 * t64;
t135 = t112 * t141;
t81 = t100 * t114 + t106 * t111;
t89 = t107 * t114 - t111 * t136;
t33 = (-t135 * t89 - t81) * pkin(9) + (t89 * t90 - t101) * pkin(3) + t132;
t137 = 0.2e1 * qJD(3) * t89 + t111 * t64 + t114 * t63;
t80 = -t100 * t111 + t106 * t114;
t82 = -pkin(3) * t135 - pkin(9) * t90;
t88 = t89 ^ 2;
t37 = -pkin(3) * t88 + pkin(9) * t80 + t135 * t82 + t137;
t147 = t117 * t33 + t149 * t37;
t144 = t112 * t118;
t143 = t112 * t121;
t110 = sin(pkin(12));
t113 = cos(pkin(12));
t129 = -g(3) * t143 - t118 * t97 + t121 * t145;
t62 = -t106 * pkin(2) - t105 * qJ(3) + t136 * t98 + qJDD(3) - t129;
t126 = -t80 * pkin(3) - t88 * pkin(9) + t82 * t90 + t62;
t116 = sin(qJ(6));
t120 = cos(qJ(6));
t103 = qJD(4) - t135;
t102 = t103 ^ 2;
t74 = t117 * t90 - t149 * t89;
t75 = t117 * t89 + t149 * t90;
t57 = pkin(4) * t74 - qJ(5) * t75;
t93 = qJDD(4) - t101;
t25 = -pkin(4) * t102 + qJ(5) * t93 - t57 * t74 + t147;
t51 = qJD(4) * t75 + t117 * t81 - t149 * t80;
t52 = -qJD(4) * t74 + t117 * t80 + t149 * t81;
t28 = (t103 * t74 - t52) * qJ(5) + (t103 * t75 + t51) * pkin(4) + t126;
t69 = t103 * t110 + t113 * t75;
t133 = -0.2e1 * qJD(5) * t69 - t110 * t25 + t113 * t28;
t45 = t110 * t93 + t113 * t52;
t68 = t103 * t113 - t110 * t75;
t19 = (t68 * t74 - t45) * pkin(10) + (t68 * t69 + t51) * pkin(5) + t133;
t138 = 0.2e1 * qJD(5) * t68 + t110 * t28 + t113 * t25;
t44 = -t110 * t52 + t113 * t93;
t55 = pkin(5) * t74 - pkin(10) * t69;
t67 = t68 ^ 2;
t20 = -t67 * pkin(5) + t44 * pkin(10) - t55 * t74 + t138;
t46 = -t116 * t69 + t120 * t68;
t31 = t46 * qJD(6) + t116 * t44 + t120 * t45;
t47 = t116 * t68 + t120 * t69;
t38 = -mrSges(7,1) * t46 + mrSges(7,2) * t47;
t73 = qJD(6) + t74;
t41 = -mrSges(7,2) * t73 + t46 * mrSges(7,3);
t50 = qJDD(6) + t51;
t17 = m(7) * (-t116 * t20 + t120 * t19) - t31 * mrSges(7,3) + t50 * mrSges(7,1) - t47 * t38 + t73 * t41;
t30 = -t47 * qJD(6) - t116 * t45 + t120 * t44;
t42 = mrSges(7,1) * t73 - t47 * mrSges(7,3);
t18 = m(7) * (t116 * t19 + t120 * t20) + t30 * mrSges(7,3) - t50 * mrSges(7,2) + t46 * t38 - t73 * t42;
t48 = -mrSges(6,1) * t68 + mrSges(6,2) * t69;
t53 = -mrSges(6,2) * t74 + mrSges(6,3) * t68;
t14 = m(6) * t133 + t51 * mrSges(6,1) - t45 * mrSges(6,3) + t116 * t18 + t120 * t17 - t48 * t69 + t53 * t74;
t54 = mrSges(6,1) * t74 - mrSges(6,3) * t69;
t15 = m(6) * t138 - t51 * mrSges(6,2) + t44 * mrSges(6,3) - t116 * t17 + t120 * t18 + t48 * t68 - t54 * t74;
t70 = -mrSges(5,2) * t103 - mrSges(5,3) * t74;
t71 = mrSges(5,1) * t103 - mrSges(5,3) * t75;
t127 = m(5) * t126 + t51 * mrSges(5,1) + t52 * mrSges(5,2) + t110 * t15 + t113 * t14 + t70 * t74 + t71 * t75;
t78 = mrSges(4,2) * t135 + mrSges(4,3) * t89;
t79 = -mrSges(4,1) * t135 - mrSges(4,3) * t90;
t124 = m(4) * t62 - mrSges(4,1) * t80 + mrSges(4,2) * t81 - t78 * t89 + t79 * t90 + t127;
t95 = -mrSges(3,2) * t107 + mrSges(3,3) * t135;
t99 = (-mrSges(3,1) * t121 + mrSges(3,2) * t118) * t142;
t10 = m(3) * t129 + mrSges(3,1) * t106 - mrSges(3,3) * t100 + t107 * t95 - t136 * t99 - t124;
t58 = mrSges(5,1) * t74 + mrSges(5,2) * t75;
t11 = m(5) * t147 - mrSges(5,2) * t93 - t51 * mrSges(5,3) - t103 * t71 - t110 * t14 + t113 * t15 - t58 * t74;
t131 = -t117 * t37 + t149 * t33;
t24 = -pkin(4) * t93 - qJ(5) * t102 + t57 * t75 + qJDD(5) - t131;
t128 = t30 * mrSges(7,1) + t46 * t41 - m(7) * (-t44 * pkin(5) - t67 * pkin(10) + t55 * t69 + t24) - t31 * mrSges(7,2) - t47 * t42;
t125 = m(6) * t24 - t44 * mrSges(6,1) + t45 * mrSges(6,2) - t53 * t68 + t54 * t69 - t128;
t16 = m(5) * t131 + mrSges(5,1) * t93 - t52 * mrSges(5,3) + t103 * t70 - t58 * t75 - t125;
t76 = -mrSges(4,1) * t89 + mrSges(4,2) * t90;
t7 = m(4) * t132 - mrSges(4,1) * t101 - mrSges(4,3) * t81 + t11 * t117 - t135 * t78 + t149 * t16 - t76 * t90;
t8 = m(4) * t137 + mrSges(4,2) * t101 + mrSges(4,3) * t80 + t11 * t149 - t117 * t16 + t135 * t79 + t76 * t89;
t94 = mrSges(3,1) * t107 - mrSges(3,3) * t136;
t4 = m(3) * (-g(3) * t144 + t146) + t101 * mrSges(3,3) - t106 * mrSges(3,2) + t99 * t135 - t107 * t94 + t114 * t8 - t111 * t7;
t6 = m(3) * (-t112 * t96 - t148) + t100 * mrSges(3,2) - t101 * mrSges(3,1) + t111 * t8 + t114 * t7 + (t118 * t94 - t121 * t95) * t142;
t140 = t10 * t143 + t115 * t6 + t144 * t4;
t2 = m(2) * t130 - mrSges(2,1) * t123 - qJDD(1) * mrSges(2,2) - t10 * t118 + t121 * t4;
t1 = m(2) * t134 + qJDD(1) * mrSges(2,1) - t123 * mrSges(2,2) - t112 * t6 + (t10 * t121 + t118 * t4) * t115;
t3 = [-m(1) * g(1) - t1 * t119 + t122 * t2, t2, t4, t8, t11, t15, t18; -m(1) * g(2) + t1 * t122 + t119 * t2, t1, t10, t7, t16, t14, t17; (-m(1) - m(2)) * g(3) + t140, -m(2) * g(3) + t140, t6, t124, t127, t125, -t128;];
f_new  = t3;
