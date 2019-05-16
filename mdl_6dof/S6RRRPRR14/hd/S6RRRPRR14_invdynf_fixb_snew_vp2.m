% Calculate vector of cutting forces with Newton-Euler
% S6RRRPRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2019-05-07 16:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPRR14_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR14_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR14_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR14_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR14_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR14_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR14_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 16:33:09
% EndTime: 2019-05-07 16:33:25
% DurationCPUTime: 5.47s
% Computational Cost: add. (69269->212), mult. (148420->279), div. (0->0), fcn. (115591->12), ass. (0->111)
t106 = sin(pkin(6));
t111 = sin(qJ(2));
t115 = cos(qJ(2));
t135 = qJD(1) * qJD(2);
t94 = (-qJDD(1) * t115 + t111 * t135) * t106;
t110 = sin(qJ(3));
t153 = cos(qJ(3));
t107 = cos(pkin(6));
t103 = t107 * qJD(1) + qJD(2);
t101 = t103 ^ 2;
t102 = t107 * qJDD(1) + qJDD(2);
t137 = qJD(1) * t115;
t117 = qJD(1) ^ 2;
t112 = sin(qJ(1));
t116 = cos(qJ(1));
t132 = t112 * g(1) - t116 * g(2);
t152 = pkin(8) * t106;
t89 = qJDD(1) * pkin(1) + t117 * t152 + t132;
t142 = t107 * t89;
t129 = -t116 * g(1) - t112 * g(2);
t90 = -t117 * pkin(1) + qJDD(1) * t152 + t129;
t143 = t111 * t142 + t115 * t90;
t138 = qJD(1) * t106;
t92 = (-pkin(2) * t115 - pkin(9) * t111) * t138;
t46 = -t101 * pkin(2) + t102 * pkin(9) + (-g(3) * t111 + t92 * t137) * t106 + t143;
t151 = t107 * g(3);
t93 = (qJDD(1) * t111 + t115 * t135) * t106;
t47 = t94 * pkin(2) - t93 * pkin(9) - t151 + (-t89 + (pkin(2) * t111 - pkin(9) * t115) * t103 * qJD(1)) * t106;
t146 = t110 * t47 + t153 * t46;
t134 = t111 * t138;
t81 = -t153 * t103 + t110 * t134;
t82 = t110 * t103 + t153 * t134;
t62 = t81 * pkin(3) - t82 * qJ(4);
t86 = qJDD(3) + t94;
t133 = t106 * t137;
t99 = -qJD(3) + t133;
t98 = t99 ^ 2;
t155 = t98 * pkin(3) - t86 * qJ(4) + 0.2e1 * qJD(4) * t99 + t81 * t62 - t146;
t109 = sin(qJ(5));
t108 = sin(qJ(6));
t113 = cos(qJ(6));
t114 = cos(qJ(5));
t150 = t81 * t99;
t130 = -t110 * t46 + t153 * t47;
t31 = -t86 * pkin(3) - t98 * qJ(4) + t82 * t62 + qJDD(4) - t130;
t61 = -t81 * qJD(3) + t110 * t102 + t153 * t93;
t22 = (t81 * t82 - t86) * pkin(10) + (t61 - t150) * pkin(4) + t31;
t139 = t106 * t115;
t128 = -g(3) * t139 - t111 * t90 + t115 * t142;
t45 = -t102 * pkin(2) - t101 * pkin(9) + t92 * t134 - t128;
t118 = (-t61 - t150) * qJ(4) + t45 + (-t99 * pkin(3) - 0.2e1 * qJD(4)) * t82;
t60 = t82 * qJD(3) - t153 * t102 + t110 * t93;
t73 = t82 * pkin(4) + t99 * pkin(10);
t80 = t81 ^ 2;
t26 = -t80 * pkin(4) - t82 * t73 + (pkin(3) + pkin(10)) * t60 + t118;
t131 = -t109 * t26 + t114 * t22;
t67 = t109 * t99 + t114 * t81;
t38 = t67 * qJD(5) + t109 * t60 + t114 * t86;
t58 = qJDD(5) + t61;
t68 = t109 * t81 - t114 * t99;
t79 = qJD(5) + t82;
t16 = (t67 * t79 - t38) * pkin(11) + (t67 * t68 + t58) * pkin(5) + t131;
t147 = t109 * t22 + t114 * t26;
t37 = -t68 * qJD(5) - t109 * t86 + t114 * t60;
t55 = t79 * pkin(5) - t68 * pkin(11);
t66 = t67 ^ 2;
t17 = -t66 * pkin(5) + t37 * pkin(11) - t79 * t55 + t147;
t48 = -t108 * t68 + t113 * t67;
t29 = t48 * qJD(6) + t108 * t37 + t113 * t38;
t49 = t108 * t67 + t113 * t68;
t35 = -t48 * mrSges(7,1) + t49 * mrSges(7,2);
t77 = qJD(6) + t79;
t39 = -t77 * mrSges(7,2) + t48 * mrSges(7,3);
t56 = qJDD(6) + t58;
t14 = m(7) * (-t108 * t17 + t113 * t16) - t29 * mrSges(7,3) + t56 * mrSges(7,1) - t49 * t35 + t77 * t39;
t28 = -t49 * qJD(6) - t108 * t38 + t113 * t37;
t40 = t77 * mrSges(7,1) - t49 * mrSges(7,3);
t15 = m(7) * (t108 * t16 + t113 * t17) + t28 * mrSges(7,3) - t56 * mrSges(7,2) + t48 * t35 - t77 * t40;
t50 = -t67 * mrSges(6,1) + t68 * mrSges(6,2);
t53 = -t79 * mrSges(6,2) + t67 * mrSges(6,3);
t11 = m(6) * t131 + t58 * mrSges(6,1) - t38 * mrSges(6,3) + t108 * t15 + t113 * t14 - t68 * t50 + t79 * t53;
t54 = t79 * mrSges(6,1) - t68 * mrSges(6,3);
t12 = m(6) * t147 - t58 * mrSges(6,2) + t37 * mrSges(6,3) - t108 * t14 + t113 * t15 + t67 * t50 - t79 * t54;
t72 = t82 * mrSges(5,1) - t99 * mrSges(5,2);
t126 = t109 * t11 - t114 * t12 - m(5) * (t60 * pkin(3) + t118) + t61 * mrSges(5,3) + t82 * t72;
t71 = t81 * mrSges(5,1) + t99 * mrSges(5,3);
t144 = t99 * mrSges(4,2) - t81 * mrSges(4,3) - t71;
t149 = mrSges(4,1) - mrSges(5,2);
t70 = -t99 * mrSges(4,1) - t82 * mrSges(4,3);
t154 = m(4) * t45 + t61 * mrSges(4,2) + t144 * t81 + t149 * t60 + t82 * t70 - t126;
t148 = -mrSges(4,3) - mrSges(5,1);
t64 = -t81 * mrSges(5,2) - t82 * mrSges(5,3);
t145 = -t81 * mrSges(4,1) - t82 * mrSges(4,2) - t64;
t140 = t106 * t111;
t122 = -t60 * pkin(4) - t80 * pkin(10) - t99 * t73 - t155;
t123 = -t28 * mrSges(7,1) - t48 * t39 + m(7) * (-t37 * pkin(5) - t66 * pkin(11) + t68 * t55 + t122) + t29 * mrSges(7,2) + t49 * t40;
t120 = m(6) * t122 - t37 * mrSges(6,1) + t38 * mrSges(6,2) - t67 * t53 + t68 * t54 + t123;
t119 = -m(5) * t155 + t120;
t13 = (t70 - t72) * t99 + (-mrSges(4,2) + mrSges(5,3)) * t86 + t145 * t81 + t148 * t60 + t119 + m(4) * t146;
t87 = t103 * mrSges(3,1) - mrSges(3,3) * t134;
t124 = -m(5) * t31 - t109 * t12 - t114 * t11;
t9 = m(4) * t130 - t144 * t99 + t145 * t82 + t148 * t61 + t149 * t86 + t124;
t91 = (-mrSges(3,1) * t115 + mrSges(3,2) * t111) * t138;
t4 = m(3) * (-g(3) * t140 + t143) - t94 * mrSges(3,3) - t102 * mrSges(3,2) + t91 * t133 - t103 * t87 + t153 * t13 - t110 * t9;
t88 = -t103 * mrSges(3,2) + mrSges(3,3) * t133;
t6 = m(3) * (-t106 * t89 - t151) + t93 * mrSges(3,2) + t94 * mrSges(3,1) + t110 * t13 + t153 * t9 + (t111 * t87 - t115 * t88) * t138;
t8 = m(3) * t128 + t102 * mrSges(3,1) - t93 * mrSges(3,3) + t103 * t88 - t91 * t134 - t154;
t136 = t107 * t6 + t8 * t139 + t4 * t140;
t2 = m(2) * t129 - t117 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t111 * t8 + t115 * t4;
t1 = m(2) * t132 + qJDD(1) * mrSges(2,1) - t117 * mrSges(2,2) - t106 * t6 + (t111 * t4 + t115 * t8) * t107;
t3 = [-m(1) * g(1) - t112 * t1 + t116 * t2, t2, t4, t13, -t60 * mrSges(5,2) - t81 * t71 - t126, t12, t15; -m(1) * g(2) + t116 * t1 + t112 * t2, t1, t8, t9, t60 * mrSges(5,1) - t86 * mrSges(5,3) + t81 * t64 + t99 * t72 - t119, t11, t14; (-m(1) - m(2)) * g(3) + t136, -m(2) * g(3) + t136, t6, t154, t61 * mrSges(5,1) + t86 * mrSges(5,2) + t82 * t64 - t99 * t71 - t124, t120, t123;];
f_new  = t3;
