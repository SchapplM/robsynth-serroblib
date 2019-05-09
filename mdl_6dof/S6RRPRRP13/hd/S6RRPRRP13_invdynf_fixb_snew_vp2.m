% Calculate vector of cutting forces with Newton-Euler
% S6RRPRRP13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-05-06 19:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRRP13_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP13_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP13_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP13_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP13_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP13_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP13_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 19:06:10
% EndTime: 2019-05-06 19:06:18
% DurationCPUTime: 2.68s
% Computational Cost: add. (30134->213), mult. (67232->269), div. (0->0), fcn. (48449->10), ass. (0->104)
t157 = -2 * qJD(3);
t106 = sin(qJ(2));
t102 = sin(pkin(6));
t137 = qJD(1) * t102;
t127 = t106 * t137;
t103 = cos(pkin(6));
t97 = t103 * qJD(1) + qJD(2);
t156 = (pkin(2) * t97 + t157) * t127;
t110 = cos(qJ(2));
t139 = t102 * t106;
t112 = qJD(1) ^ 2;
t107 = sin(qJ(1));
t111 = cos(qJ(1));
t126 = t107 * g(1) - t111 * g(2);
t80 = t112 * t102 * pkin(8) + qJDD(1) * pkin(1) + t126;
t142 = t103 * t80;
t123 = -t111 * g(1) - t107 * g(2);
t134 = qJDD(1) * t102;
t81 = -t112 * pkin(1) + pkin(8) * t134 + t123;
t121 = -g(3) * t139 + t106 * t142 + t110 * t81;
t136 = qJD(1) * t110;
t128 = t102 * t136;
t82 = (-pkin(2) * t110 - qJ(3) * t106) * t137;
t95 = t97 ^ 2;
t96 = t103 * qJDD(1) + qJDD(2);
t155 = t95 * pkin(2) - t96 * qJ(3) - t82 * t128 + t97 * t157 - t121;
t105 = sin(qJ(4));
t109 = cos(qJ(4));
t140 = t102 ^ 2 * t112;
t129 = t110 ^ 2 * t140;
t152 = t103 * g(3);
t153 = -pkin(2) - pkin(9);
t85 = pkin(3) * t127 - t97 * pkin(9);
t86 = (qJD(2) * t136 + qJDD(1) * t106) * t102;
t87 = -qJD(2) * t127 + t110 * t134;
t30 = -pkin(3) * t129 - t152 - t86 * qJ(3) + t153 * t87 + (-t80 + (-qJ(3) * t110 * t97 - t106 * t85) * qJD(1)) * t102 + t156;
t138 = t102 * t110;
t145 = g(3) * t138 + t106 * t81;
t119 = -t95 * qJ(3) + t82 * t127 + qJDD(3) + t145;
t32 = t86 * pkin(3) + t153 * t96 + (-pkin(3) * t97 * t137 - pkin(9) * t106 * t140 - t142) * t110 + t119;
t124 = -t105 * t30 + t109 * t32;
t70 = -t105 * t97 - t109 * t128;
t71 = -t105 * t128 + t109 * t97;
t58 = -t70 * pkin(4) - t71 * pkin(10);
t75 = qJDD(4) + t86;
t91 = qJD(4) + t127;
t89 = t91 ^ 2;
t21 = -t75 * pkin(4) - t89 * pkin(10) + t71 * t58 - t124;
t104 = sin(qJ(5));
t108 = cos(qJ(5));
t56 = t70 * qJD(4) - t105 * t87 + t109 * t96;
t61 = t104 * t91 + t108 * t71;
t35 = -t61 * qJD(5) - t104 * t56 + t108 * t75;
t60 = -t104 * t71 + t108 * t91;
t36 = t60 * qJD(5) + t104 * t75 + t108 * t56;
t69 = qJD(5) - t70;
t48 = t69 * pkin(5) - t61 * qJ(6);
t49 = t69 * mrSges(7,1) - t61 * mrSges(7,3);
t59 = t60 ^ 2;
t131 = m(7) * (-t35 * pkin(5) - t59 * qJ(6) + t61 * t48 + qJDD(6) + t21) + t36 * mrSges(7,2) + t61 * t49;
t46 = -t69 * mrSges(7,2) + t60 * mrSges(7,3);
t47 = -t69 * mrSges(6,2) + t60 * mrSges(6,3);
t50 = t69 * mrSges(6,1) - t61 * mrSges(6,3);
t154 = m(6) * t21 + t36 * mrSges(6,2) - (t47 + t46) * t60 - (mrSges(6,1) + mrSges(7,1)) * t35 + t61 * t50 + t131;
t151 = mrSges(3,1) - mrSges(4,2);
t149 = mrSges(3,3) + mrSges(4,1);
t147 = t105 * t32 + t109 * t30;
t22 = -t89 * pkin(4) + t75 * pkin(10) + t70 * t58 + t147;
t113 = t87 * pkin(3) - pkin(9) * t129 + t97 * t85 - t155;
t55 = -t71 * qJD(4) - t105 * t96 - t109 * t87;
t25 = (-t70 * t91 - t56) * pkin(10) + (t71 * t91 - t55) * pkin(4) + t113;
t148 = t104 * t25 + t108 * t22;
t79 = mrSges(4,1) * t127 + t97 * mrSges(4,2);
t144 = t97 * mrSges(3,1) - mrSges(3,3) * t127 - t79;
t83 = (mrSges(4,2) * t110 - mrSges(4,3) * t106) * t137;
t143 = t83 + (-mrSges(3,1) * t110 + mrSges(3,2) * t106) * t137;
t125 = -t104 * t22 + t108 * t25;
t53 = qJDD(5) - t55;
t133 = m(7) * (-0.2e1 * qJD(6) * t61 + (t60 * t69 - t36) * qJ(6) + (t60 * t61 + t53) * pkin(5) + t125) + t69 * t46 + t53 * mrSges(7,1);
t43 = -t60 * mrSges(7,1) + t61 * mrSges(7,2);
t44 = -t60 * mrSges(6,1) + t61 * mrSges(6,2);
t12 = m(6) * t125 + t53 * mrSges(6,1) + t69 * t47 + (-t44 - t43) * t61 + (-mrSges(6,3) - mrSges(7,3)) * t36 + t133;
t132 = m(7) * (-t59 * pkin(5) + t35 * qJ(6) + 0.2e1 * qJD(6) * t60 - t69 * t48 + t148) + t35 * mrSges(7,3) + t60 * t43;
t14 = m(6) * t148 + t35 * mrSges(6,3) + t60 * t44 + (-t50 - t49) * t69 + (-mrSges(6,2) - mrSges(7,2)) * t53 + t132;
t57 = -t70 * mrSges(5,1) + t71 * mrSges(5,2);
t63 = t91 * mrSges(5,1) - t71 * mrSges(5,3);
t10 = m(5) * t147 - t75 * mrSges(5,2) + t55 * mrSges(5,3) - t104 * t12 + t108 * t14 + t70 * t57 - t91 * t63;
t122 = -t102 * t80 - t152;
t62 = -t91 * mrSges(5,2) + t70 * mrSges(5,3);
t15 = m(5) * t124 + t75 * mrSges(5,1) - t56 * mrSges(5,3) - t71 * t57 + t91 * t62 - t154;
t78 = -mrSges(4,1) * t128 - t97 * mrSges(4,3);
t120 = t109 * t10 - t105 * t15 + m(4) * (-t87 * pkin(2) + (-t97 * t128 - t86) * qJ(3) + t122 + t156) + t78 * t128 - t86 * mrSges(4,3);
t77 = -t97 * mrSges(3,2) + mrSges(3,3) * t128;
t5 = m(3) * t122 + t86 * mrSges(3,2) - t151 * t87 + (t144 * t106 - t110 * t77) * t137 + t120;
t130 = t110 * t142;
t118 = -m(4) * (-t96 * pkin(2) + t119 - t130) - t105 * t10 - t109 * t15;
t6 = m(3) * (t130 - t145) + (t77 - t78) * t97 + t151 * t96 - t149 * t86 - t143 * t127 + t118;
t116 = m(5) * t113 - t55 * mrSges(5,1) + t56 * mrSges(5,2) + t104 * t14 + t108 * t12 - t70 * t62 + t71 * t63;
t114 = -m(4) * t155 + t116;
t8 = m(3) * t121 - t144 * t97 + (-mrSges(3,2) + mrSges(4,3)) * t96 + t149 * t87 + t143 * t128 + t114;
t135 = t103 * t5 + t6 * t138 + t8 * t139;
t2 = m(2) * t123 - t112 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t106 * t6 + t110 * t8;
t1 = m(2) * t126 + qJDD(1) * mrSges(2,1) - t112 * mrSges(2,2) - t102 * t5 + (t106 * t8 + t110 * t6) * t103;
t3 = [-m(1) * g(1) - t107 * t1 + t111 * t2, t2, t8, t87 * mrSges(4,2) - t79 * t127 + t120, t10, t14, -t53 * mrSges(7,2) - t69 * t49 + t132; -m(1) * g(2) + t111 * t1 + t107 * t2, t1, t6, -t87 * mrSges(4,1) - t96 * mrSges(4,3) - t83 * t128 - t97 * t79 - t114, t15, t12, -t36 * mrSges(7,3) - t61 * t43 + t133; (-m(1) - m(2)) * g(3) + t135, -m(2) * g(3) + t135, t5, t86 * mrSges(4,1) + t96 * mrSges(4,2) + t83 * t127 + t97 * t78 - t118, t116, t154, -t35 * mrSges(7,1) - t60 * t46 + t131;];
f_new  = t3;
