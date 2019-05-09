% Calculate vector of cutting forces with Newton-Euler
% S6RRPRPR13
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
% Datum: 2019-05-06 16:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRPR13_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR13_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR13_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR13_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR13_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR13_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR13_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 16:32:03
% EndTime: 2019-05-06 16:32:14
% DurationCPUTime: 4.89s
% Computational Cost: add. (62153->216), mult. (140996->285), div. (0->0), fcn. (103779->12), ass. (0->112)
t159 = -2 * qJD(3);
t110 = cos(pkin(6));
t102 = qJD(1) * t110 + qJD(2);
t113 = sin(qJ(2));
t108 = sin(pkin(6));
t144 = qJD(1) * t108;
t136 = t113 * t144;
t158 = (pkin(2) * t102 + t159) * t136;
t100 = t102 ^ 2;
t101 = qJDD(1) * t110 + qJDD(2);
t117 = cos(qJ(2));
t146 = t108 * t113;
t119 = qJD(1) ^ 2;
t114 = sin(qJ(1));
t118 = cos(qJ(1));
t134 = t114 * g(1) - g(2) * t118;
t85 = pkin(8) * t108 * t119 + qJDD(1) * pkin(1) + t134;
t149 = t110 * t85;
t131 = -g(1) * t118 - g(2) * t114;
t140 = qJDD(1) * t108;
t86 = -pkin(1) * t119 + pkin(8) * t140 + t131;
t129 = -g(3) * t146 + t113 * t149 + t117 * t86;
t143 = qJD(1) * t117;
t135 = t108 * t143;
t87 = (-pkin(2) * t117 - qJ(3) * t113) * t144;
t157 = t100 * pkin(2) - t101 * qJ(3) + t102 * t159 - t87 * t135 - t129;
t156 = -pkin(2) - pkin(9);
t155 = t110 * g(3);
t154 = mrSges(3,1) - mrSges(4,2);
t153 = mrSges(3,3) + mrSges(4,1);
t112 = sin(qJ(4));
t116 = cos(qJ(4));
t147 = t108 ^ 2 * t119;
t137 = t117 ^ 2 * t147;
t90 = pkin(3) * t136 - pkin(9) * t102;
t91 = (qJD(2) * t143 + qJDD(1) * t113) * t108;
t92 = -qJD(2) * t136 + t117 * t140;
t37 = -pkin(3) * t137 - t155 - t91 * qJ(3) + t156 * t92 + (-t85 + (-qJ(3) * t102 * t117 - t113 * t90) * qJD(1)) * t108 + t158;
t145 = t108 * t117;
t148 = g(3) * t145 + t113 * t86;
t127 = -t100 * qJ(3) + t87 * t136 + qJDD(3) + t148;
t39 = t91 * pkin(3) + t156 * t101 + (-pkin(3) * t102 * t144 - pkin(9) * t113 * t147 - t149) * t117 + t127;
t152 = t112 * t39 + t116 * t37;
t84 = mrSges(4,1) * t136 + mrSges(4,2) * t102;
t151 = mrSges(3,1) * t102 - mrSges(3,3) * t136 - t84;
t88 = (mrSges(4,2) * t117 - mrSges(4,3) * t113) * t144;
t150 = t88 + (-mrSges(3,1) * t117 + mrSges(3,2) * t113) * t144;
t107 = sin(pkin(11));
t109 = cos(pkin(11));
t111 = sin(qJ(6));
t115 = cos(qJ(6));
t74 = t102 * t112 + t116 * t135;
t75 = t102 * t116 - t112 * t135;
t59 = pkin(4) * t74 - qJ(5) * t75;
t80 = qJDD(4) + t91;
t96 = qJD(4) + t136;
t94 = t96 ^ 2;
t24 = -pkin(4) * t94 + qJ(5) * t80 - t59 * t74 + t152;
t121 = t92 * pkin(3) - pkin(9) * t137 + t102 * t90 - t157;
t57 = qJD(4) * t75 + t101 * t112 + t116 * t92;
t58 = -qJD(4) * t74 + t101 * t116 - t112 * t92;
t27 = (t74 * t96 - t58) * qJ(5) + (t75 * t96 + t57) * pkin(4) + t121;
t65 = t107 * t96 + t109 * t75;
t132 = -0.2e1 * qJD(5) * t65 - t107 * t24 + t109 * t27;
t48 = t107 * t80 + t109 * t58;
t64 = -t107 * t75 + t109 * t96;
t18 = (t64 * t74 - t48) * pkin(10) + (t64 * t65 + t57) * pkin(5) + t132;
t139 = 0.2e1 * qJD(5) * t64 + t107 * t27 + t109 * t24;
t47 = -t107 * t58 + t109 * t80;
t53 = pkin(5) * t74 - t65 * pkin(10);
t63 = t64 ^ 2;
t19 = -t63 * pkin(5) + t47 * pkin(10) - t53 * t74 + t139;
t45 = -t111 * t65 + t115 * t64;
t30 = t45 * qJD(6) + t111 * t47 + t115 * t48;
t46 = t111 * t64 + t115 * t65;
t32 = -mrSges(7,1) * t45 + mrSges(7,2) * t46;
t73 = qJD(6) + t74;
t40 = -mrSges(7,2) * t73 + t45 * mrSges(7,3);
t55 = qJDD(6) + t57;
t16 = m(7) * (-t111 * t19 + t115 * t18) - t30 * mrSges(7,3) + t55 * mrSges(7,1) - t46 * t32 + t73 * t40;
t29 = -t46 * qJD(6) - t111 * t48 + t115 * t47;
t41 = mrSges(7,1) * t73 - t46 * mrSges(7,3);
t17 = m(7) * (t111 * t18 + t115 * t19) + t29 * mrSges(7,3) - t55 * mrSges(7,2) + t45 * t32 - t73 * t41;
t49 = -mrSges(6,1) * t64 + mrSges(6,2) * t65;
t51 = -mrSges(6,2) * t74 + t64 * mrSges(6,3);
t13 = m(6) * t132 + t57 * mrSges(6,1) - t48 * mrSges(6,3) + t111 * t17 + t115 * t16 - t65 * t49 + t74 * t51;
t52 = mrSges(6,1) * t74 - t65 * mrSges(6,3);
t14 = m(6) * t139 - t57 * mrSges(6,2) + t47 * mrSges(6,3) - t111 * t16 + t115 * t17 + t64 * t49 - t74 * t52;
t60 = mrSges(5,1) * t74 + mrSges(5,2) * t75;
t67 = mrSges(5,1) * t96 - mrSges(5,3) * t75;
t10 = m(5) * t152 - t80 * mrSges(5,2) - t57 * mrSges(5,3) - t107 * t13 + t109 * t14 - t74 * t60 - t96 * t67;
t130 = -t108 * t85 - t155;
t133 = -t112 * t37 + t116 * t39;
t23 = -pkin(4) * t80 - qJ(5) * t94 + t75 * t59 + qJDD(5) - t133;
t125 = t29 * mrSges(7,1) + t45 * t40 - m(7) * (-t47 * pkin(5) - t63 * pkin(10) + t65 * t53 + t23) - t30 * mrSges(7,2) - t46 * t41;
t120 = m(6) * t23 - t47 * mrSges(6,1) + t48 * mrSges(6,2) - t64 * t51 + t65 * t52 - t125;
t66 = -mrSges(5,2) * t96 - mrSges(5,3) * t74;
t15 = m(5) * t133 + t80 * mrSges(5,1) - t58 * mrSges(5,3) - t75 * t60 + t96 * t66 - t120;
t83 = -mrSges(4,1) * t135 - mrSges(4,3) * t102;
t128 = t116 * t10 - t112 * t15 + m(4) * (-t92 * pkin(2) + (-t102 * t135 - t91) * qJ(3) + t130 + t158) + t83 * t135 - t91 * mrSges(4,3);
t82 = -mrSges(3,2) * t102 + mrSges(3,3) * t135;
t5 = m(3) * t130 + t91 * mrSges(3,2) - t154 * t92 + (t151 * t113 - t117 * t82) * t144 + t128;
t138 = t117 * t149;
t126 = -m(4) * (-t101 * pkin(2) + t127 - t138) - t112 * t10 - t116 * t15;
t6 = m(3) * (t138 - t148) - t153 * t91 + (t82 - t83) * t102 + t154 * t101 - t150 * t136 + t126;
t124 = m(5) * t121 + t57 * mrSges(5,1) + t58 * mrSges(5,2) + t107 * t14 + t109 * t13 + t74 * t66 + t75 * t67;
t122 = -m(4) * t157 + t124;
t8 = m(3) * t129 + t153 * t92 - t151 * t102 + (-mrSges(3,2) + mrSges(4,3)) * t101 + t150 * t135 + t122;
t141 = t110 * t5 + t6 * t145 + t8 * t146;
t2 = m(2) * t131 - t119 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t113 * t6 + t117 * t8;
t1 = m(2) * t134 + qJDD(1) * mrSges(2,1) - t119 * mrSges(2,2) - t108 * t5 + (t113 * t8 + t117 * t6) * t110;
t3 = [-m(1) * g(1) - t1 * t114 + t118 * t2, t2, t8, t92 * mrSges(4,2) - t84 * t136 + t128, t10, t14, t17; -m(1) * g(2) + t1 * t118 + t114 * t2, t1, t6, -t92 * mrSges(4,1) - t101 * mrSges(4,3) - t102 * t84 - t88 * t135 - t122, t15, t13, t16; (-m(1) - m(2)) * g(3) + t141, -m(2) * g(3) + t141, t5, t91 * mrSges(4,1) + t101 * mrSges(4,2) + t102 * t83 + t88 * t136 - t126, t124, t120, -t125;];
f_new  = t3;
