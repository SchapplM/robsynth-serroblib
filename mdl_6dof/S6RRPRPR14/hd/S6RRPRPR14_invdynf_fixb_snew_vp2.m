% Calculate vector of cutting forces with Newton-Euler
% S6RRPRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
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
% Datum: 2019-05-06 17:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRPR14_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR14_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR14_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR14_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR14_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR14_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR14_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 16:54:39
% EndTime: 2019-05-06 16:54:45
% DurationCPUTime: 2.24s
% Computational Cost: add. (22967->215), mult. (51707->268), div. (0->0), fcn. (36462->10), ass. (0->108)
t161 = -2 * qJD(3);
t109 = cos(qJ(2));
t101 = sin(pkin(6));
t111 = qJD(1) ^ 2;
t141 = t111 * t101 ^ 2;
t133 = t109 ^ 2 * t141;
t105 = sin(qJ(2));
t139 = t101 * t105;
t102 = cos(pkin(6));
t106 = sin(qJ(1));
t110 = cos(qJ(1));
t129 = t106 * g(1) - g(2) * t110;
t78 = pkin(8) * t101 * t111 + qJDD(1) * pkin(1) + t129;
t142 = t102 * t78;
t127 = -g(1) * t110 - g(2) * t106;
t134 = qJDD(1) * t101;
t79 = -pkin(1) * t111 + pkin(8) * t134 + t127;
t125 = -g(3) * t139 + t105 * t142 + t109 * t79;
t136 = qJD(1) * t109;
t130 = t101 * t136;
t137 = qJD(1) * t101;
t80 = (-pkin(2) * t109 - qJ(3) * t105) * t137;
t96 = qJD(1) * t102 + qJD(2);
t94 = t96 ^ 2;
t95 = qJDD(1) * t102 + qJDD(2);
t158 = pkin(2) * t94 - t95 * qJ(3) - t80 * t130 + t96 * t161 - t125;
t131 = t105 * t137;
t83 = pkin(3) * t131 - pkin(9) * t96;
t85 = -qJD(2) * t131 + t109 * t134;
t113 = pkin(3) * t85 - pkin(9) * t133 + t96 * t83 - t158;
t103 = sin(qJ(6));
t107 = cos(qJ(6));
t104 = sin(qJ(4));
t108 = cos(qJ(4));
t67 = t104 * t96 + t108 * t130;
t90 = qJD(4) + t131;
t153 = t67 * t90;
t157 = -2 * qJD(5);
t47 = -t67 * qJD(4) - t104 * t85 + t108 * t95;
t68 = -t104 * t130 + t108 * t96;
t112 = (-t47 + t153) * qJ(5) + t113 + (t90 * pkin(4) + t157) * t68;
t154 = g(3) * t102;
t156 = -pkin(2) - pkin(9);
t159 = (pkin(2) * t96 + t161) * t131;
t84 = (qJD(2) * t136 + qJDD(1) * t105) * t101;
t27 = -pkin(3) * t133 - t154 - qJ(3) * t84 + t156 * t85 + (-t78 + (-qJ(3) * t109 * t96 - t105 * t83) * qJD(1)) * t101 + t159;
t138 = t101 * t109;
t145 = g(3) * t138 + t105 * t79;
t122 = -qJ(3) * t94 + t80 * t131 + qJDD(3) + t145;
t29 = pkin(3) * t84 + t156 * t95 + (-pkin(3) * t96 * t137 - pkin(9) * t105 * t141 - t142) * t109 + t122;
t128 = -t104 * t27 + t108 * t29;
t48 = t67 * pkin(4) - qJ(5) * t68;
t73 = qJDD(4) + t84;
t87 = t90 ^ 2;
t20 = -pkin(4) * t73 - qJ(5) * t87 + t68 * t48 + qJDD(5) - t128;
t15 = (t67 * t68 - t73) * pkin(10) + (t47 + t153) * pkin(5) + t20;
t46 = qJD(4) * t68 + t104 * t95 + t108 * t85;
t57 = pkin(5) * t68 - pkin(10) * t90;
t66 = t67 ^ 2;
t18 = -t66 * pkin(5) + t112 + (pkin(4) + pkin(10)) * t46 - t57 * t68;
t51 = -t103 * t90 + t107 * t67;
t32 = t51 * qJD(6) + t103 * t46 + t107 * t73;
t52 = t103 * t67 + t107 * t90;
t37 = -mrSges(7,1) * t51 + mrSges(7,2) * t52;
t65 = qJD(6) + t68;
t39 = -mrSges(7,2) * t65 + mrSges(7,3) * t51;
t43 = qJDD(6) + t47;
t13 = m(7) * (-t103 * t18 + t107 * t15) - t32 * mrSges(7,3) + t43 * mrSges(7,1) - t52 * t37 + t65 * t39;
t31 = -t52 * qJD(6) - t103 * t73 + t107 * t46;
t40 = mrSges(7,1) * t65 - mrSges(7,3) * t52;
t14 = m(7) * (t103 * t15 + t107 * t18) + t31 * mrSges(7,3) - t43 * mrSges(7,2) + t51 * t37 - t65 * t40;
t54 = mrSges(6,1) * t68 + mrSges(6,2) * t90;
t123 = t103 * t13 - t107 * t14 - m(6) * (t46 * pkin(4) + t112) + t47 * mrSges(6,3) + t68 * t54;
t53 = t67 * mrSges(6,1) - mrSges(6,3) * t90;
t146 = -mrSges(5,2) * t90 - t67 * mrSges(5,3) - t53;
t151 = mrSges(5,1) - mrSges(6,2);
t56 = mrSges(5,1) * t90 - mrSges(5,3) * t68;
t114 = m(5) * t113 + t47 * mrSges(5,2) + t146 * t67 + t151 * t46 + t68 * t56 - t123;
t160 = -m(4) * t158 + t114;
t152 = mrSges(3,1) - mrSges(4,2);
t150 = mrSges(3,3) + mrSges(4,1);
t149 = -mrSges(5,3) - mrSges(6,1);
t148 = t104 * t29 + t108 * t27;
t50 = -t67 * mrSges(6,2) - mrSges(6,3) * t68;
t147 = -t67 * mrSges(5,1) - mrSges(5,2) * t68 - t50;
t77 = mrSges(4,1) * t131 + mrSges(4,2) * t96;
t144 = mrSges(3,1) * t96 - mrSges(3,3) * t131 - t77;
t81 = (mrSges(4,2) * t109 - mrSges(4,3) * t105) * t137;
t143 = t81 + (-mrSges(3,1) * t109 + mrSges(3,2) * t105) * t137;
t117 = -pkin(4) * t87 + qJ(5) * t73 - t67 * t48 + t148;
t119 = -t31 * mrSges(7,1) - t51 * t39 + m(7) * (-t46 * pkin(5) - t66 * pkin(10) + ((2 * qJD(5)) + t57) * t90 + t117) + t32 * mrSges(7,2) + t52 * t40;
t115 = -m(6) * (t90 * t157 - t117) + t119;
t11 = m(5) * t148 + (-t56 + t54) * t90 + (-mrSges(5,2) + mrSges(6,3)) * t73 + t147 * t67 + t149 * t46 + t115;
t126 = -t101 * t78 - t154;
t76 = -mrSges(4,1) * t130 - mrSges(4,3) * t96;
t120 = -m(6) * t20 - t103 * t14 - t107 * t13;
t9 = m(5) * t128 + t146 * t90 + t147 * t68 + t149 * t47 + t151 * t73 + t120;
t124 = -t104 * t9 + t108 * t11 + m(4) * (-pkin(2) * t85 + (-t96 * t130 - t84) * qJ(3) + t126 + t159) + t76 * t130 - t84 * mrSges(4,3);
t75 = -mrSges(3,2) * t96 + mrSges(3,3) * t130;
t5 = m(3) * t126 + t84 * mrSges(3,2) - t152 * t85 + (t144 * t105 - t109 * t75) * t137 + t124;
t132 = t109 * t142;
t121 = -m(4) * (-pkin(2) * t95 + t122 - t132) - t104 * t11 - t108 * t9;
t6 = m(3) * (t132 - t145) + (t75 - t76) * t96 + t152 * t95 - t150 * t84 - t143 * t131 + t121;
t8 = m(3) * t125 + t150 * t85 + (-mrSges(3,2) + mrSges(4,3)) * t95 - t144 * t96 + t143 * t130 + t160;
t135 = t102 * t5 + t6 * t138 + t8 * t139;
t2 = m(2) * t127 - t111 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t105 * t6 + t109 * t8;
t1 = m(2) * t129 + qJDD(1) * mrSges(2,1) - t111 * mrSges(2,2) - t101 * t5 + (t105 * t8 + t109 * t6) * t102;
t3 = [-m(1) * g(1) - t1 * t106 + t110 * t2, t2, t8, t85 * mrSges(4,2) - t77 * t131 + t124, t11, -t46 * mrSges(6,2) - t67 * t53 - t123, t14; -m(1) * g(2) + t1 * t110 + t106 * t2, t1, t6, -t85 * mrSges(4,1) - t95 * mrSges(4,3) - t81 * t130 - t96 * t77 - t160, t9, t46 * mrSges(6,1) - t73 * mrSges(6,3) + t67 * t50 - t90 * t54 - t115, t13; (-m(1) - m(2)) * g(3) + t135, -m(2) * g(3) + t135, t5, t84 * mrSges(4,1) + t95 * mrSges(4,2) + t81 * t131 + t96 * t76 - t121, t114, t47 * mrSges(6,1) + t73 * mrSges(6,2) + t68 * t50 + t90 * t53 - t120, t119;];
f_new  = t3;
