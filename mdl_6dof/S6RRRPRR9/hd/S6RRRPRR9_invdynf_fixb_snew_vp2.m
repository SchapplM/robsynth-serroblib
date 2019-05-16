% Calculate vector of cutting forces with Newton-Euler
% S6RRRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-05-07 13:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPRR9_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR9_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR9_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR9_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_invdynf_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR9_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR9_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR9_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 12:53:07
% EndTime: 2019-05-07 12:54:08
% DurationCPUTime: 21.71s
% Computational Cost: add. (382566->227), mult. (1003467->327), div. (0->0), fcn. (851766->16), ass. (0->129)
t105 = sin(pkin(13));
t108 = cos(pkin(13));
t113 = sin(qJ(3));
t118 = cos(qJ(3));
t106 = sin(pkin(7));
t147 = t106 * t118;
t109 = cos(pkin(7));
t107 = sin(pkin(6));
t114 = sin(qJ(2));
t119 = cos(qJ(2));
t139 = qJD(1) * qJD(2);
t100 = (qJDD(1) * t114 + t119 * t139) * t107;
t110 = cos(pkin(6));
t102 = t110 * qJDD(1) + qJDD(2);
t103 = t110 * qJD(1) + qJD(2);
t121 = qJD(1) ^ 2;
t115 = sin(qJ(1));
t120 = cos(qJ(1));
t134 = t115 * g(1) - t120 * g(2);
t157 = pkin(9) * t107;
t97 = qJDD(1) * pkin(1) + t121 * t157 + t134;
t150 = t110 * t97;
t129 = -t120 * g(1) - t115 * g(2);
t98 = -t121 * pkin(1) + qJDD(1) * t157 + t129;
t131 = -t114 * t98 + t119 * t150;
t142 = qJD(1) * t114;
t155 = pkin(10) * t109;
t141 = qJD(1) * t119;
t135 = t107 * t141;
t90 = (t103 * t106 + t109 * t135) * pkin(10);
t143 = qJD(1) * t107;
t156 = pkin(10) * t106;
t94 = (-pkin(2) * t119 - t114 * t156) * t143;
t59 = -t100 * t155 + t102 * pkin(2) + t103 * t90 + (-g(3) * t119 - t142 * t94) * t107 + t131;
t151 = t109 * t59;
t101 = (qJDD(1) * t119 - t114 * t139) * t107;
t126 = t101 * t109 + t102 * t106;
t152 = t114 * t150 + t119 * t98;
t136 = t107 * t142;
t93 = t103 * pkin(2) - t136 * t155;
t60 = -t103 * t93 + (-g(3) * t114 + t141 * t94) * t107 + t126 * pkin(10) + t152;
t154 = t110 * g(3);
t64 = -t100 * t156 - t101 * pkin(2) - t154 + (-t97 + (t114 * t93 - t119 * t90) * qJD(1)) * t107;
t130 = -t113 * t60 + t118 * t151 + t64 * t147;
t144 = t109 * t119;
t84 = t103 * t147 + (-t113 * t114 + t118 * t144) * t143;
t72 = t84 * qJD(3) + t118 * t100 + t113 * t126;
t148 = t106 * t113;
t85 = t103 * t148 + (t113 * t144 + t114 * t118) * t143;
t86 = -t106 * t101 + t109 * t102 + qJDD(3);
t91 = t109 * t103 - t106 * t135 + qJD(3);
t31 = (t84 * t91 - t72) * qJ(4) + (t84 * t85 + t86) * pkin(3) + t130;
t137 = t113 * t151 + t118 * t60 + t64 * t148;
t71 = -t85 * qJD(3) - t113 * t100 + t118 * t126;
t81 = t91 * pkin(3) - t85 * qJ(4);
t83 = t84 ^ 2;
t34 = -t83 * pkin(3) + t71 * qJ(4) - t91 * t81 + t137;
t78 = t105 * t84 + t108 * t85;
t158 = -0.2e1 * qJD(4) * t78 - t105 * t34 + t108 * t31;
t112 = sin(qJ(5));
t117 = cos(qJ(5));
t77 = -t105 * t85 + t108 * t84;
t138 = 0.2e1 * qJD(4) * t77 + t105 * t31 + t108 * t34;
t55 = -t77 * pkin(4) - t78 * pkin(11);
t89 = t91 ^ 2;
t25 = -t89 * pkin(4) + t86 * pkin(11) + t77 * t55 + t138;
t133 = -t106 * t59 + t109 * t64;
t124 = -t71 * pkin(3) - t83 * qJ(4) + t85 * t81 + qJDD(4) + t133;
t49 = -t105 * t72 + t108 * t71;
t50 = t105 * t71 + t108 * t72;
t27 = (-t77 * t91 - t50) * pkin(11) + (t78 * t91 - t49) * pkin(4) + t124;
t153 = t112 * t27 + t117 * t25;
t146 = t107 * t114;
t145 = t107 * t119;
t111 = sin(qJ(6));
t116 = cos(qJ(6));
t66 = -t112 * t78 + t117 * t91;
t67 = t112 * t91 + t117 * t78;
t44 = -t66 * pkin(5) - t67 * pkin(12);
t46 = qJDD(5) - t49;
t76 = qJD(5) - t77;
t75 = t76 ^ 2;
t21 = -t75 * pkin(5) + t46 * pkin(12) + t66 * t44 + t153;
t24 = -t86 * pkin(4) - t89 * pkin(11) + t78 * t55 - t158;
t38 = -t67 * qJD(5) - t112 * t50 + t117 * t86;
t39 = t66 * qJD(5) + t112 * t86 + t117 * t50;
t22 = (-t66 * t76 - t39) * pkin(12) + (t67 * t76 - t38) * pkin(5) + t24;
t47 = -t111 * t67 + t116 * t76;
t29 = t47 * qJD(6) + t111 * t46 + t116 * t39;
t48 = t111 * t76 + t116 * t67;
t36 = -t47 * mrSges(7,1) + t48 * mrSges(7,2);
t37 = qJDD(6) - t38;
t65 = qJD(6) - t66;
t40 = -t65 * mrSges(7,2) + t47 * mrSges(7,3);
t18 = m(7) * (-t111 * t21 + t116 * t22) - t29 * mrSges(7,3) + t37 * mrSges(7,1) - t48 * t36 + t65 * t40;
t28 = -t48 * qJD(6) - t111 * t39 + t116 * t46;
t41 = t65 * mrSges(7,1) - t48 * mrSges(7,3);
t19 = m(7) * (t111 * t22 + t116 * t21) + t28 * mrSges(7,3) - t37 * mrSges(7,2) + t47 * t36 - t65 * t41;
t43 = -t66 * mrSges(6,1) + t67 * mrSges(6,2);
t52 = t76 * mrSges(6,1) - t67 * mrSges(6,3);
t15 = m(6) * t153 - t46 * mrSges(6,2) + t38 * mrSges(6,3) - t111 * t18 + t116 * t19 + t66 * t43 - t76 * t52;
t127 = -t112 * t25 + t117 * t27;
t123 = m(7) * (-t46 * pkin(5) - t75 * pkin(12) + t67 * t44 - t127) - t28 * mrSges(7,1) + t29 * mrSges(7,2) - t47 * t40 + t48 * t41;
t51 = -t76 * mrSges(6,2) + t66 * mrSges(6,3);
t17 = m(6) * t127 + t46 * mrSges(6,1) - t39 * mrSges(6,3) - t67 * t43 + t76 * t51 - t123;
t68 = -t91 * mrSges(5,2) + t77 * mrSges(5,3);
t69 = t91 * mrSges(5,1) - t78 * mrSges(5,3);
t125 = m(5) * t124 - t49 * mrSges(5,1) + t50 * mrSges(5,2) + t112 * t15 + t117 * t17 - t77 * t68 + t78 * t69;
t80 = -t91 * mrSges(4,2) + t84 * mrSges(4,3);
t82 = t91 * mrSges(4,1) - t85 * mrSges(4,3);
t12 = m(4) * t133 - t71 * mrSges(4,1) + t72 * mrSges(4,2) - t84 * t80 + t85 * t82 + t125;
t54 = -t77 * mrSges(5,1) + t78 * mrSges(5,2);
t11 = m(5) * t138 - t86 * mrSges(5,2) + t49 * mrSges(5,3) - t112 * t17 + t117 * t15 + t77 * t54 - t91 * t69;
t122 = m(6) * t24 - t38 * mrSges(6,1) + t39 * mrSges(6,2) + t111 * t19 + t116 * t18 - t66 * t51 + t67 * t52;
t13 = m(5) * t158 + t86 * mrSges(5,1) - t50 * mrSges(5,3) - t78 * t54 + t91 * t68 - t122;
t79 = -t84 * mrSges(4,1) + t85 * mrSges(4,2);
t10 = m(4) * t137 - t86 * mrSges(4,2) + t71 * mrSges(4,3) - t105 * t13 + t108 * t11 + t84 * t79 - t91 * t82;
t9 = m(4) * t130 + t86 * mrSges(4,1) - t72 * mrSges(4,3) + t105 * t11 + t108 * t13 - t85 * t79 + t91 * t80;
t128 = t113 * t10 + t118 * t9;
t96 = -t103 * mrSges(3,2) + mrSges(3,3) * t135;
t99 = (-mrSges(3,1) * t119 + mrSges(3,2) * t114) * t143;
t4 = m(3) * (-g(3) * t145 + t131) - t100 * mrSges(3,3) + t102 * mrSges(3,1) - t99 * t136 + t103 * t96 - t106 * t12 + t128 * t109;
t95 = t103 * mrSges(3,1) - mrSges(3,3) * t136;
t6 = m(3) * (-t107 * t97 - t154) + t100 * mrSges(3,2) - t101 * mrSges(3,1) + t109 * t12 + t128 * t106 + (t114 * t95 - t119 * t96) * t143;
t8 = m(3) * (-g(3) * t146 + t152) + t101 * mrSges(3,3) - t102 * mrSges(3,2) + t99 * t135 - t103 * t95 + t118 * t10 - t113 * t9;
t140 = t110 * t6 + t4 * t145 + t8 * t146;
t2 = m(2) * t129 - t121 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t114 * t4 + t119 * t8;
t1 = m(2) * t134 + qJDD(1) * mrSges(2,1) - t121 * mrSges(2,2) - t107 * t6 + (t114 * t8 + t119 * t4) * t110;
t3 = [-m(1) * g(1) - t115 * t1 + t120 * t2, t2, t8, t10, t11, t15, t19; -m(1) * g(2) + t120 * t1 + t115 * t2, t1, t4, t9, t13, t17, t18; (-m(1) - m(2)) * g(3) + t140, -m(2) * g(3) + t140, t6, t12, t125, t122, t123;];
f_new  = t3;
