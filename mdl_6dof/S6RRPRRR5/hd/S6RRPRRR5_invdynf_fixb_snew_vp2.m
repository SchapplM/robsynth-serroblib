% Calculate vector of cutting forces with Newton-Euler
% S6RRPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-06 21:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRRR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 21:16:43
% EndTime: 2019-05-06 21:17:03
% DurationCPUTime: 9.32s
% Computational Cost: add. (155411->213), mult. (406458->295), div. (0->0), fcn. (325574->14), ass. (0->117)
t154 = -2 * qJD(3);
t109 = sin(pkin(12));
t111 = cos(pkin(12));
t116 = sin(qJ(2));
t121 = cos(qJ(2));
t110 = sin(pkin(6));
t144 = qJD(1) * t110;
t92 = (t109 * t116 - t111 * t121) * t144;
t142 = qJD(1) * qJD(2);
t101 = (qJDD(1) * t116 + t121 * t142) * t110;
t112 = cos(pkin(6));
t104 = qJDD(1) * t112 + qJDD(2);
t105 = qJD(1) * t112 + qJD(2);
t123 = qJD(1) ^ 2;
t117 = sin(qJ(1));
t122 = cos(qJ(1));
t137 = t117 * g(1) - g(2) * t122;
t152 = pkin(8) * t110;
t98 = qJDD(1) * pkin(1) + t123 * t152 + t137;
t149 = t112 * t98;
t132 = -g(1) * t122 - g(2) * t117;
t99 = -pkin(1) * t123 + qJDD(1) * t152 + t132;
t134 = -t116 * t99 + t121 * t149;
t147 = t110 ^ 2 * t123;
t51 = pkin(2) * t104 - qJ(3) * t101 + (pkin(2) * t116 * t147 + (qJ(3) * qJD(1) * t105 - g(3)) * t110) * t121 + t134;
t102 = (qJDD(1) * t121 - t116 * t142) * t110;
t146 = t110 * t116;
t130 = -g(3) * t146 + t116 * t149 + t121 * t99;
t140 = t121 ^ 2 * t147;
t139 = t116 * t144;
t95 = pkin(2) * t105 - qJ(3) * t139;
t55 = -pkin(2) * t140 + qJ(3) * t102 - t105 * t95 + t130;
t93 = (t109 * t121 + t111 * t116) * t144;
t153 = -t109 * t55 + t111 * t51 + t93 * t154;
t114 = sin(qJ(5));
t119 = cos(qJ(5));
t115 = sin(qJ(4));
t120 = cos(qJ(4));
t103 = t105 ^ 2;
t141 = t109 * t51 + t111 * t55 + t92 * t154;
t71 = pkin(3) * t92 - pkin(9) * t93;
t35 = -pkin(3) * t103 + pkin(9) * t104 - t71 * t92 + t141;
t131 = -g(3) * t112 - t110 * t98;
t126 = -pkin(2) * t102 - qJ(3) * t140 + t95 * t139 + qJDD(3) + t131;
t75 = -t101 * t109 + t102 * t111;
t76 = t101 * t111 + t102 * t109;
t39 = (t105 * t92 - t76) * pkin(9) + (t105 * t93 - t75) * pkin(3) + t126;
t150 = t115 * t39 + t120 * t35;
t80 = t105 * t120 - t115 * t93;
t81 = t105 * t115 + t120 * t93;
t62 = -pkin(4) * t80 - pkin(10) * t81;
t74 = qJDD(4) - t75;
t91 = qJD(4) + t92;
t90 = t91 ^ 2;
t25 = -pkin(4) * t90 + pkin(10) * t74 + t62 * t80 + t150;
t34 = -pkin(3) * t104 - pkin(9) * t103 + t93 * t71 - t153;
t58 = -t81 * qJD(4) + t104 * t120 - t115 * t76;
t59 = qJD(4) * t80 + t104 * t115 + t120 * t76;
t28 = (-t80 * t91 - t59) * pkin(10) + (t81 * t91 - t58) * pkin(4) + t34;
t151 = t114 * t28 + t119 * t25;
t145 = t110 * t121;
t113 = sin(qJ(6));
t118 = cos(qJ(6));
t135 = -t114 * t25 + t119 * t28;
t65 = -t114 * t81 + t119 * t91;
t41 = qJD(5) * t65 + t114 * t74 + t119 * t59;
t57 = qJDD(5) - t58;
t66 = t114 * t91 + t119 * t81;
t79 = qJD(5) - t80;
t19 = (t65 * t79 - t41) * pkin(11) + (t65 * t66 + t57) * pkin(5) + t135;
t40 = -qJD(5) * t66 - t114 * t59 + t119 * t74;
t54 = pkin(5) * t79 - pkin(11) * t66;
t64 = t65 ^ 2;
t20 = -pkin(5) * t64 + pkin(11) * t40 - t54 * t79 + t151;
t44 = -t113 * t66 + t118 * t65;
t31 = qJD(6) * t44 + t113 * t40 + t118 * t41;
t45 = t113 * t65 + t118 * t66;
t37 = -mrSges(7,1) * t44 + mrSges(7,2) * t45;
t77 = qJD(6) + t79;
t42 = -mrSges(7,2) * t77 + mrSges(7,3) * t44;
t56 = qJDD(6) + t57;
t17 = m(7) * (-t113 * t20 + t118 * t19) - t31 * mrSges(7,3) + t56 * mrSges(7,1) - t45 * t37 + t77 * t42;
t30 = -qJD(6) * t45 - t113 * t41 + t118 * t40;
t43 = mrSges(7,1) * t77 - mrSges(7,3) * t45;
t18 = m(7) * (t113 * t19 + t118 * t20) + t30 * mrSges(7,3) - t56 * mrSges(7,2) + t44 * t37 - t77 * t43;
t46 = -mrSges(6,1) * t65 + mrSges(6,2) * t66;
t52 = -mrSges(6,2) * t79 + mrSges(6,3) * t65;
t13 = m(6) * t135 + t57 * mrSges(6,1) - t41 * mrSges(6,3) + t113 * t18 + t118 * t17 - t66 * t46 + t79 * t52;
t53 = mrSges(6,1) * t79 - mrSges(6,3) * t66;
t14 = m(6) * t151 - t57 * mrSges(6,2) + t40 * mrSges(6,3) - t113 * t17 + t118 * t18 + t65 * t46 - t79 * t53;
t67 = -mrSges(5,2) * t91 + mrSges(5,3) * t80;
t68 = mrSges(5,1) * t91 - mrSges(5,3) * t81;
t125 = m(5) * t34 - t58 * mrSges(5,1) + t59 * mrSges(5,2) + t114 * t14 + t119 * t13 - t80 * t67 + t81 * t68;
t70 = mrSges(4,1) * t92 + mrSges(4,2) * t93;
t82 = -mrSges(4,2) * t105 - mrSges(4,3) * t92;
t10 = m(4) * t153 + t104 * mrSges(4,1) - t76 * mrSges(4,3) + t105 * t82 - t93 * t70 - t125;
t100 = (-mrSges(3,1) * t121 + mrSges(3,2) * t116) * t144;
t61 = -mrSges(5,1) * t80 + mrSges(5,2) * t81;
t12 = m(5) * t150 - t74 * mrSges(5,2) + t58 * mrSges(5,3) - t114 * t13 + t119 * t14 + t80 * t61 - t91 * t68;
t133 = -t115 * t35 + t120 * t39;
t24 = -pkin(4) * t74 - pkin(10) * t90 + t81 * t62 - t133;
t128 = t30 * mrSges(7,1) + t44 * t42 - m(7) * (-pkin(5) * t40 - pkin(11) * t64 + t54 * t66 + t24) - t31 * mrSges(7,2) - t45 * t43;
t124 = m(6) * t24 - t40 * mrSges(6,1) + t41 * mrSges(6,2) - t65 * t52 + t66 * t53 - t128;
t16 = m(5) * t133 + t74 * mrSges(5,1) - t59 * mrSges(5,3) - t81 * t61 + t91 * t67 - t124;
t83 = mrSges(4,1) * t105 - mrSges(4,3) * t93;
t7 = m(4) * t141 - t104 * mrSges(4,2) + t75 * mrSges(4,3) - t105 * t83 - t115 * t16 + t120 * t12 - t92 * t70;
t138 = t121 * t144;
t97 = -mrSges(3,2) * t105 + mrSges(3,3) * t138;
t5 = m(3) * (-g(3) * t145 + t134) - t101 * mrSges(3,3) + t104 * mrSges(3,1) - t100 * t139 + t105 * t97 + t109 * t7 + t111 * t10;
t96 = mrSges(3,1) * t105 - mrSges(3,3) * t139;
t6 = m(3) * t130 - t104 * mrSges(3,2) + t102 * mrSges(3,3) - t109 * t10 + t100 * t138 - t105 * t96 + t111 * t7;
t127 = m(4) * t126 - t75 * mrSges(4,1) + t76 * mrSges(4,2) + t115 * t12 + t120 * t16 + t92 * t82 + t93 * t83;
t9 = m(3) * t131 + t101 * mrSges(3,2) - t102 * mrSges(3,1) + (t116 * t96 - t121 * t97) * t144 + t127;
t143 = t112 * t9 + t5 * t145 + t6 * t146;
t2 = m(2) * t132 - t123 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t116 * t5 + t121 * t6;
t1 = m(2) * t137 + qJDD(1) * mrSges(2,1) - t123 * mrSges(2,2) - t110 * t9 + (t116 * t6 + t121 * t5) * t112;
t3 = [-m(1) * g(1) - t1 * t117 + t122 * t2, t2, t6, t7, t12, t14, t18; -m(1) * g(2) + t1 * t122 + t117 * t2, t1, t5, t10, t16, t13, t17; (-m(1) - m(2)) * g(3) + t143, -m(2) * g(3) + t143, t9, t127, t125, t124, -t128;];
f_new  = t3;
