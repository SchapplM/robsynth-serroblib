% Calculate vector of cutting forces with Newton-Euler
% S6RRPRRR10
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
% Datum: 2019-05-06 23:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRRR10_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR10_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR10_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR10_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR10_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR10_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR10_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 23:36:54
% EndTime: 2019-05-06 23:37:19
% DurationCPUTime: 9.39s
% Computational Cost: add. (172639->215), mult. (390390->295), div. (0->0), fcn. (318973->14), ass. (0->115)
t114 = cos(pkin(6));
t150 = t114 * g(3);
t116 = sin(qJ(5));
t121 = cos(qJ(5));
t112 = sin(pkin(6));
t123 = cos(qJ(2));
t142 = qJD(1) * t123;
t137 = t112 * t142;
t104 = qJD(4) - t137;
t103 = t104 ^ 2;
t117 = sin(qJ(4));
t122 = cos(qJ(4));
t118 = sin(qJ(2));
t143 = qJD(1) * t112;
t138 = t118 * t143;
t140 = qJDD(1) * t112;
t102 = -qJD(2) * t138 + t123 * t140;
t111 = sin(pkin(12));
t113 = cos(pkin(12));
t108 = qJD(1) * t114 + qJD(2);
t106 = t108 ^ 2;
t107 = qJDD(1) * t114 + qJDD(2);
t125 = qJD(1) ^ 2;
t119 = sin(qJ(1));
t124 = cos(qJ(1));
t136 = t119 * g(1) - g(2) * t124;
t97 = pkin(8) * t112 * t125 + qJDD(1) * pkin(1) + t136;
t146 = t114 * t97;
t132 = -g(1) * t124 - g(2) * t119;
t98 = -pkin(1) * t125 + pkin(8) * t140 + t132;
t147 = t118 * t146 + t123 * t98;
t99 = (-pkin(2) * t123 - qJ(3) * t118) * t143;
t64 = -t106 * pkin(2) + t107 * qJ(3) + (-g(3) * t118 + t99 * t142) * t112 + t147;
t101 = (qJD(2) * t142 + qJDD(1) * t118) * t112;
t65 = -t102 * pkin(2) - t150 - t101 * qJ(3) + (-t97 + (pkin(2) * t118 - qJ(3) * t123) * t108 * qJD(1)) * t112;
t91 = t108 * t111 + t113 * t138;
t133 = -0.2e1 * qJD(3) * t91 - t111 * t64 + t113 * t65;
t82 = t101 * t113 + t107 * t111;
t90 = t108 * t113 - t111 * t138;
t33 = (-t90 * t137 - t82) * pkin(9) + (t90 * t91 - t102) * pkin(3) + t133;
t139 = 0.2e1 * qJD(3) * t90 + t111 * t65 + t113 * t64;
t81 = -t101 * t111 + t107 * t113;
t83 = -pkin(3) * t137 - pkin(9) * t91;
t89 = t90 ^ 2;
t37 = -pkin(3) * t89 + pkin(9) * t81 + t83 * t137 + t139;
t148 = t117 * t33 + t122 * t37;
t75 = -t117 * t91 + t122 * t90;
t76 = t117 * t90 + t122 * t91;
t59 = -pkin(4) * t75 - pkin(10) * t76;
t94 = qJDD(4) - t102;
t25 = -pkin(4) * t103 + pkin(10) * t94 + t59 * t75 + t148;
t144 = t112 * t123;
t131 = -g(3) * t144 - t118 * t98 + t123 * t146;
t63 = -t107 * pkin(2) - t106 * qJ(3) + t99 * t138 + qJDD(3) - t131;
t128 = -t81 * pkin(3) - t89 * pkin(9) + t91 * t83 + t63;
t52 = -t76 * qJD(4) - t117 * t82 + t122 * t81;
t53 = qJD(4) * t75 + t117 * t81 + t122 * t82;
t31 = (-t104 * t75 - t53) * pkin(10) + (t104 * t76 - t52) * pkin(4) + t128;
t149 = t116 * t31 + t121 * t25;
t145 = t112 * t118;
t100 = (-mrSges(3,1) * t123 + mrSges(3,2) * t118) * t143;
t115 = sin(qJ(6));
t120 = cos(qJ(6));
t135 = -t116 * t25 + t121 * t31;
t67 = t104 * t121 - t116 * t76;
t40 = t67 * qJD(5) + t116 * t94 + t121 * t53;
t51 = qJDD(5) - t52;
t68 = t104 * t116 + t121 * t76;
t74 = qJD(5) - t75;
t19 = (t67 * t74 - t40) * pkin(11) + (t67 * t68 + t51) * pkin(5) + t135;
t39 = -t68 * qJD(5) - t116 * t53 + t121 * t94;
t56 = pkin(5) * t74 - t68 * pkin(11);
t66 = t67 ^ 2;
t20 = -t66 * pkin(5) + t39 * pkin(11) - t56 * t74 + t149;
t46 = -t115 * t68 + t120 * t67;
t28 = t46 * qJD(6) + t115 * t39 + t120 * t40;
t47 = t115 * t67 + t120 * t68;
t38 = -mrSges(7,1) * t46 + mrSges(7,2) * t47;
t72 = qJD(6) + t74;
t43 = -mrSges(7,2) * t72 + t46 * mrSges(7,3);
t49 = qJDD(6) + t51;
t17 = m(7) * (-t115 * t20 + t120 * t19) - t28 * mrSges(7,3) + t49 * mrSges(7,1) - t47 * t38 + t72 * t43;
t27 = -t47 * qJD(6) - t115 * t40 + t120 * t39;
t44 = mrSges(7,1) * t72 - t47 * mrSges(7,3);
t18 = m(7) * (t115 * t19 + t120 * t20) + t27 * mrSges(7,3) - t49 * mrSges(7,2) + t46 * t38 - t72 * t44;
t48 = -mrSges(6,1) * t67 + mrSges(6,2) * t68;
t54 = -mrSges(6,2) * t74 + t67 * mrSges(6,3);
t14 = m(6) * t135 + t51 * mrSges(6,1) - t40 * mrSges(6,3) + t115 * t18 + t120 * t17 - t68 * t48 + t74 * t54;
t55 = mrSges(6,1) * t74 - t68 * mrSges(6,3);
t15 = m(6) * t149 - t51 * mrSges(6,2) + t39 * mrSges(6,3) - t115 * t17 + t120 * t18 + t67 * t48 - t74 * t55;
t69 = -mrSges(5,2) * t104 + mrSges(5,3) * t75;
t70 = mrSges(5,1) * t104 - mrSges(5,3) * t76;
t129 = -m(5) * t128 + t52 * mrSges(5,1) - t53 * mrSges(5,2) - t116 * t15 - t121 * t14 + t75 * t69 - t76 * t70;
t79 = mrSges(4,2) * t137 + mrSges(4,3) * t90;
t80 = -mrSges(4,1) * t137 - mrSges(4,3) * t91;
t126 = m(4) * t63 - t81 * mrSges(4,1) + t82 * mrSges(4,2) - t90 * t79 + t91 * t80 - t129;
t96 = -mrSges(3,2) * t108 + mrSges(3,3) * t137;
t10 = m(3) * t131 + t107 * mrSges(3,1) - t101 * mrSges(3,3) - t100 * t138 + t108 * t96 - t126;
t58 = -mrSges(5,1) * t75 + mrSges(5,2) * t76;
t11 = m(5) * t148 - t94 * mrSges(5,2) + t52 * mrSges(5,3) - t104 * t70 - t116 * t14 + t121 * t15 + t75 * t58;
t134 = -t117 * t37 + t122 * t33;
t24 = -pkin(4) * t94 - pkin(10) * t103 + t76 * t59 - t134;
t130 = t27 * mrSges(7,1) + t46 * t43 - m(7) * (-t39 * pkin(5) - t66 * pkin(11) + t68 * t56 + t24) - t28 * mrSges(7,2) - t47 * t44;
t127 = m(6) * t24 - t39 * mrSges(6,1) + t40 * mrSges(6,2) - t67 * t54 + t68 * t55 - t130;
t16 = m(5) * t134 + t94 * mrSges(5,1) - t53 * mrSges(5,3) + t104 * t69 - t76 * t58 - t127;
t77 = -mrSges(4,1) * t90 + mrSges(4,2) * t91;
t7 = m(4) * t133 - t102 * mrSges(4,1) - t82 * mrSges(4,3) + t117 * t11 + t122 * t16 - t79 * t137 - t91 * t77;
t8 = m(4) * t139 + t102 * mrSges(4,2) + t81 * mrSges(4,3) + t122 * t11 - t117 * t16 + t80 * t137 + t90 * t77;
t95 = mrSges(3,1) * t108 - mrSges(3,3) * t138;
t4 = m(3) * (-g(3) * t145 + t147) + t102 * mrSges(3,3) - t107 * mrSges(3,2) + t100 * t137 - t108 * t95 + t113 * t8 - t111 * t7;
t6 = m(3) * (-t112 * t97 - t150) + t101 * mrSges(3,2) - t102 * mrSges(3,1) + t111 * t8 + t113 * t7 + (t118 * t95 - t123 * t96) * t143;
t141 = t10 * t144 + t114 * t6 + t4 * t145;
t2 = m(2) * t132 - t125 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t118 * t10 + t123 * t4;
t1 = m(2) * t136 + qJDD(1) * mrSges(2,1) - t125 * mrSges(2,2) - t112 * t6 + (t10 * t123 + t118 * t4) * t114;
t3 = [-m(1) * g(1) - t1 * t119 + t124 * t2, t2, t4, t8, t11, t15, t18; -m(1) * g(2) + t1 * t124 + t119 * t2, t1, t10, t7, t16, t14, t17; (-m(1) - m(2)) * g(3) + t141, -m(2) * g(3) + t141, t6, t126, -t129, t127, -t130;];
f_new  = t3;
