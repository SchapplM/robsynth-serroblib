% Calculate vector of cutting forces with Newton-Euler
% S6RRPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-05-06 10:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPPRR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 10:19:41
% EndTime: 2019-05-06 10:19:51
% DurationCPUTime: 4.20s
% Computational Cost: add. (48356->213), mult. (129877->281), div. (0->0), fcn. (99165->12), ass. (0->112)
t155 = -2 * qJD(4);
t106 = cos(pkin(6));
t100 = t106 * qJD(1) + qJD(2);
t104 = sin(pkin(11));
t144 = cos(pkin(11));
t105 = sin(pkin(6));
t109 = sin(qJ(2));
t113 = cos(qJ(2));
t115 = qJD(1) ^ 2;
t110 = sin(qJ(1));
t114 = cos(qJ(1));
t132 = t110 * g(1) - t114 * g(2);
t153 = pkin(8) * t105;
t90 = qJDD(1) * pkin(1) + t115 * t153 + t132;
t145 = t106 * t90;
t130 = -t114 * g(1) - t110 * g(2);
t91 = -t115 * pkin(1) + qJDD(1) * t153 + t130;
t131 = -t109 * t91 + t113 * t145;
t142 = t105 ^ 2 * t115;
t136 = qJD(1) * qJD(2);
t93 = (qJDD(1) * t109 + t113 * t136) * t105;
t99 = t106 * qJDD(1) + qJDD(2);
t36 = t99 * pkin(2) - t93 * qJ(3) + (pkin(2) * t109 * t142 + (qJ(3) * qJD(1) * t100 - g(3)) * t105) * t113 + t131;
t141 = t105 * t109;
t126 = -g(3) * t141 + t109 * t145 + t113 * t91;
t135 = t113 ^ 2 * t142;
t139 = qJD(1) * t105;
t134 = t109 * t139;
t87 = t100 * pkin(2) - qJ(3) * t134;
t94 = (qJDD(1) * t113 - t109 * t136) * t105;
t39 = -pkin(2) * t135 + t94 * qJ(3) - t100 * t87 + t126;
t149 = t104 * t36 + t144 * t39;
t133 = t113 * t139;
t83 = t104 * t134 - t144 * t133;
t84 = (t104 * t113 + t144 * t109) * t139;
t56 = t83 * pkin(3) - t84 * qJ(4);
t98 = t100 ^ 2;
t154 = t98 * pkin(3) - t99 * qJ(4) + t100 * t155 + t83 * t56 - t149;
t152 = mrSges(4,1) - mrSges(5,2);
t151 = -mrSges(4,3) - mrSges(5,1);
t108 = sin(qJ(5));
t112 = cos(qJ(5));
t146 = t100 * t83;
t129 = -t104 * t39 + t144 * t36;
t28 = -t99 * pkin(3) - t98 * qJ(4) + qJDD(4) - t129 + ((2 * qJD(3)) + t56) * t84;
t63 = t104 * t94 + t144 * t93;
t22 = (t83 * t84 - t99) * pkin(9) + (t63 + t146) * pkin(4) + t28;
t128 = -t106 * g(3) - t105 * t90;
t118 = -t94 * pkin(2) - qJ(3) * t135 + t87 * t134 + qJDD(3) + t128;
t116 = (-t63 + t146) * qJ(4) + t118 + (pkin(3) * t100 + t155) * t84;
t62 = t104 * t93 - t144 * t94;
t72 = t84 * pkin(4) - t100 * pkin(9);
t82 = t83 ^ 2;
t26 = -t82 * pkin(4) - t84 * t72 + (pkin(3) + pkin(9)) * t62 + t116;
t150 = t108 * t22 + t112 * t26;
t58 = -t83 * mrSges(5,2) - t84 * mrSges(5,3);
t148 = -t83 * mrSges(4,1) - t84 * mrSges(4,2) - t58;
t70 = t83 * mrSges(5,1) - t100 * mrSges(5,3);
t147 = -t100 * mrSges(4,2) - t83 * mrSges(4,3) - t70;
t143 = qJD(3) * t83;
t140 = t105 * t113;
t107 = sin(qJ(6));
t111 = cos(qJ(6));
t76 = -0.2e1 * t143;
t119 = -t62 * pkin(4) - t82 * pkin(9) + t100 * t72 - t154 + t76;
t66 = -t108 * t100 + t112 * t83;
t67 = t112 * t100 + t108 * t83;
t47 = -t66 * pkin(5) - t67 * pkin(10);
t61 = qJDD(5) + t63;
t81 = qJD(5) + t84;
t80 = t81 ^ 2;
t19 = -t80 * pkin(5) + t61 * pkin(10) + t66 * t47 + t150;
t43 = -t67 * qJD(5) - t108 * t99 + t112 * t62;
t44 = t66 * qJD(5) + t108 * t62 + t112 * t99;
t20 = (-t66 * t81 - t44) * pkin(10) + (t67 * t81 - t43) * pkin(5) + t119;
t50 = -t107 * t67 + t111 * t81;
t31 = t50 * qJD(6) + t107 * t61 + t111 * t44;
t51 = t107 * t81 + t111 * t67;
t32 = -t50 * mrSges(7,1) + t51 * mrSges(7,2);
t65 = qJD(6) - t66;
t37 = -t65 * mrSges(7,2) + t50 * mrSges(7,3);
t42 = qJDD(6) - t43;
t16 = m(7) * (-t107 * t19 + t111 * t20) - t31 * mrSges(7,3) + t42 * mrSges(7,1) - t51 * t32 + t65 * t37;
t30 = -t51 * qJD(6) - t107 * t44 + t111 * t61;
t38 = t65 * mrSges(7,1) - t51 * mrSges(7,3);
t17 = m(7) * (t107 * t20 + t111 * t19) + t30 * mrSges(7,3) - t42 * mrSges(7,2) + t50 * t32 - t65 * t38;
t52 = -t81 * mrSges(6,2) + t66 * mrSges(6,3);
t53 = t81 * mrSges(6,1) - t67 * mrSges(6,3);
t122 = m(6) * t119 - t43 * mrSges(6,1) + t44 * mrSges(6,2) + t107 * t17 + t111 * t16 - t66 * t52 + t67 * t53;
t121 = -m(5) * (0.2e1 * t143 + t154) + t122;
t69 = t100 * mrSges(4,1) - t84 * mrSges(4,3);
t71 = t84 * mrSges(5,1) + t100 * mrSges(5,2);
t10 = m(4) * (t76 + t149) + (-mrSges(4,2) + mrSges(5,3)) * t99 + t148 * t83 + t151 * t62 + (-t69 + t71) * t100 + t121;
t46 = -t66 * mrSges(6,1) + t67 * mrSges(6,2);
t12 = m(6) * t150 - t61 * mrSges(6,2) + t43 * mrSges(6,3) - t107 * t16 + t111 * t17 + t66 * t46 - t81 * t53;
t127 = -t108 * t26 + t112 * t22;
t120 = m(7) * (-t61 * pkin(5) - t80 * pkin(10) + t67 * t47 - t127) - t30 * mrSges(7,1) + t31 * mrSges(7,2) - t50 * t37 + t51 * t38;
t13 = m(6) * t127 + t61 * mrSges(6,1) - t44 * mrSges(6,3) - t67 * t46 + t81 * t52 - t120;
t123 = -m(5) * t28 - t108 * t12 - t112 * t13;
t7 = m(4) * t129 + t152 * t99 + (-0.2e1 * m(4) * qJD(3) + t148) * t84 + t151 * t63 + t147 * t100 + t123;
t89 = -t100 * mrSges(3,2) + mrSges(3,3) * t133;
t92 = (-mrSges(3,1) * t113 + mrSges(3,2) * t109) * t139;
t5 = m(3) * (-g(3) * t140 + t131) - t93 * mrSges(3,3) + t99 * mrSges(3,1) - t92 * t134 + t100 * t89 + t104 * t10 + t144 * t7;
t88 = t100 * mrSges(3,1) - mrSges(3,3) * t134;
t6 = m(3) * t126 - t99 * mrSges(3,2) + t94 * mrSges(3,3) + t144 * t10 - t100 * t88 - t104 * t7 + t92 * t133;
t124 = -t108 * t13 + t112 * t12 + m(5) * (t62 * pkin(3) + t116) - t84 * t71 - t63 * mrSges(5,3);
t117 = m(4) * t118 + t63 * mrSges(4,2) + t147 * t83 + t152 * t62 + t84 * t69 + t124;
t9 = (t109 * t88 - t113 * t89) * t139 + m(3) * t128 + t93 * mrSges(3,2) - t94 * mrSges(3,1) + t117;
t137 = t106 * t9 + t5 * t140 + t6 * t141;
t2 = m(2) * t130 - t115 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t109 * t5 + t113 * t6;
t1 = m(2) * t132 + qJDD(1) * mrSges(2,1) - t115 * mrSges(2,2) - t105 * t9 + (t109 * t6 + t113 * t5) * t106;
t3 = [-m(1) * g(1) - t110 * t1 + t114 * t2, t2, t6, t10, -t62 * mrSges(5,2) - t83 * t70 + t124, t12, t17; -m(1) * g(2) + t114 * t1 + t110 * t2, t1, t5, t7, t62 * mrSges(5,1) - t99 * mrSges(5,3) - t100 * t71 + t83 * t58 - t121, t13, t16; (-m(1) - m(2)) * g(3) + t137, -m(2) * g(3) + t137, t9, t117, t63 * mrSges(5,1) + t99 * mrSges(5,2) + t100 * t70 + t84 * t58 - t123, t122, t120;];
f_new  = t3;
