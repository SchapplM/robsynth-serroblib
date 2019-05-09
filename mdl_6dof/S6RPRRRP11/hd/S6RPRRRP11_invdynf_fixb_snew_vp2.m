% Calculate vector of cutting forces with Newton-Euler
% S6RPRRRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 02:13
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRRP11_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP11_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP11_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP11_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP11_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP11_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP11_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 02:01:26
% EndTime: 2019-05-06 02:01:46
% DurationCPUTime: 9.79s
% Computational Cost: add. (153018->215), mult. (475463->305), div. (0->0), fcn. (404596->14), ass. (0->117)
t100 = sin(pkin(6));
t106 = sin(qJ(3));
t110 = cos(qJ(3));
t101 = cos(pkin(12));
t102 = cos(pkin(7));
t145 = t101 * t102;
t103 = cos(pkin(6));
t99 = sin(pkin(7));
t150 = t103 * t99;
t98 = sin(pkin(12));
t114 = (-t106 * t98 + t110 * t145) * t100 + t110 * t150;
t79 = t114 * qJD(1);
t144 = t102 * t106;
t149 = t106 * t99;
t116 = t103 * t149 + (t101 * t144 + t110 * t98) * t100;
t80 = t116 * qJD(1);
t70 = -t80 * qJD(3) + t114 * qJDD(1);
t112 = qJD(1) ^ 2;
t107 = sin(qJ(1));
t111 = cos(qJ(1));
t127 = -t111 * g(1) - t107 * g(2);
t142 = qJD(1) * t100;
t147 = qJ(2) * t100;
t159 = -t112 * pkin(1) + 0.2e1 * qJD(2) * t142 + qJDD(1) * t147 + t127;
t123 = -pkin(9) * t98 * t99 - pkin(2) * t101;
t148 = pkin(9) * qJDD(1);
t118 = qJD(1) * t123 * t142 + t102 * t148;
t134 = t107 * g(1) - t111 * g(2);
t90 = qJDD(1) * pkin(1) + t112 * t147 + t134;
t151 = t103 * t90;
t128 = t101 * t151 - t159 * t98;
t87 = (t100 * t145 + t150) * qJD(1) * pkin(9);
t54 = (pkin(2) * qJDD(1) + qJD(1) * t87) * t103 + (-g(3) * t101 - t118 * t98) * t100 + t128;
t136 = t159 * t101 + t98 * t151;
t152 = t100 * t98;
t92 = (-pkin(9) * t102 * t152 + pkin(2) * t103) * qJD(1);
t55 = (-qJD(1) * t92 + t99 * t148) * t103 + (-g(3) * t98 + t118 * t101) * t100 + t136;
t129 = -t103 * g(3) + qJDD(2);
t64 = (-t90 + t123 * qJDD(1) + (-t101 * t87 + t92 * t98) * qJD(1)) * t100 + t129;
t158 = -t106 * t55 + (t102 * t54 + t64 * t99) * t110;
t105 = sin(qJ(4));
t109 = cos(qJ(4));
t137 = t110 * t55 + t54 * t144 + t64 * t149;
t69 = -t79 * pkin(3) - t80 * pkin(10);
t146 = t100 * t101;
t120 = t102 * t103 - t99 * t146;
t88 = t120 * qJD(1) + qJD(3);
t84 = t88 ^ 2;
t85 = t120 * qJDD(1) + qJDD(3);
t29 = -t84 * pkin(3) + t85 * pkin(10) + t79 * t69 + t137;
t133 = t102 * t64 - t99 * t54;
t71 = t79 * qJD(3) + t116 * qJDD(1);
t31 = (-t79 * t88 - t71) * pkin(10) + (t80 * t88 - t70) * pkin(3) + t133;
t130 = -t105 * t29 + t109 * t31;
t73 = -t105 * t80 + t109 * t88;
t74 = t105 * t88 + t109 * t80;
t57 = -t73 * pkin(4) - t74 * pkin(11);
t67 = qJDD(4) - t70;
t78 = qJD(4) - t79;
t77 = t78 ^ 2;
t21 = -t67 * pkin(4) - t77 * pkin(11) + t74 * t57 - t130;
t104 = sin(qJ(5));
t108 = cos(qJ(5));
t50 = t73 * qJD(4) + t105 * t85 + t109 * t71;
t63 = t104 * t78 + t108 * t74;
t34 = -t63 * qJD(5) - t104 * t50 + t108 * t67;
t62 = -t104 * t74 + t108 * t78;
t35 = t62 * qJD(5) + t104 * t67 + t108 * t50;
t72 = qJD(5) - t73;
t44 = t72 * pkin(5) - t63 * qJ(6);
t45 = t72 * mrSges(7,1) - t63 * mrSges(7,3);
t61 = t62 ^ 2;
t138 = m(7) * (-t34 * pkin(5) - t61 * qJ(6) + t63 * t44 + qJDD(6) + t21) + t35 * mrSges(7,2) + t63 * t45;
t42 = -t72 * mrSges(7,2) + t62 * mrSges(7,3);
t43 = -t72 * mrSges(6,2) + t62 * mrSges(6,3);
t46 = t72 * mrSges(6,1) - t63 * mrSges(6,3);
t157 = m(6) * t21 + t35 * mrSges(6,2) - (t43 + t42) * t62 - (mrSges(6,1) + mrSges(7,1)) * t34 + t63 * t46 + t138;
t154 = t105 * t31 + t109 * t29;
t22 = -t77 * pkin(4) + t67 * pkin(11) + t73 * t57 + t154;
t28 = -t85 * pkin(3) - t84 * pkin(10) + t80 * t69 - t158;
t49 = -t74 * qJD(4) - t105 * t71 + t109 * t85;
t25 = (-t73 * t78 - t50) * pkin(11) + (t74 * t78 - t49) * pkin(4) + t28;
t155 = t104 * t25 + t108 * t22;
t131 = -t104 * t22 + t108 * t25;
t48 = qJDD(5) - t49;
t140 = m(7) * (-0.2e1 * qJD(6) * t63 + (t62 * t72 - t35) * qJ(6) + (t62 * t63 + t48) * pkin(5) + t131) + t72 * t42 + t48 * mrSges(7,1);
t39 = -t62 * mrSges(7,1) + t63 * mrSges(7,2);
t40 = -t62 * mrSges(6,1) + t63 * mrSges(6,2);
t13 = m(6) * t131 + t48 * mrSges(6,1) + t72 * t43 + (-t40 - t39) * t63 + (-mrSges(6,3) - mrSges(7,3)) * t35 + t140;
t139 = m(7) * (-t61 * pkin(5) + t34 * qJ(6) + 0.2e1 * qJD(6) * t62 - t72 * t44 + t155) + t34 * mrSges(7,3) + t62 * t39;
t14 = m(6) * t155 + t34 * mrSges(6,3) + t62 * t40 + (-t46 - t45) * t72 + (-mrSges(6,2) - mrSges(7,2)) * t48 + t139;
t56 = -t73 * mrSges(5,1) + t74 * mrSges(5,2);
t66 = t78 * mrSges(5,1) - t74 * mrSges(5,3);
t12 = m(5) * t154 - t67 * mrSges(5,2) + t49 * mrSges(5,3) - t104 * t13 + t108 * t14 + t73 * t56 - t78 * t66;
t65 = -t78 * mrSges(5,2) + t73 * mrSges(5,3);
t15 = m(5) * t130 + t67 * mrSges(5,1) - t50 * mrSges(5,3) - t74 * t56 + t78 * t65 - t157;
t75 = -t88 * mrSges(4,2) + t79 * mrSges(4,3);
t76 = t88 * mrSges(4,1) - t80 * mrSges(4,3);
t10 = m(4) * t133 - t70 * mrSges(4,1) + t71 * mrSges(4,2) + t105 * t12 + t109 * t15 - t79 * t75 + t80 * t76;
t122 = t103 * mrSges(3,1) - mrSges(3,3) * t152;
t113 = m(5) * t28 - t49 * mrSges(5,1) + t50 * mrSges(5,2) + t104 * t14 + t108 * t13 - t73 * t65 + t74 * t66;
t68 = -t79 * mrSges(4,1) + t80 * mrSges(4,2);
t11 = m(4) * t158 + t85 * mrSges(4,1) - t71 * mrSges(4,3) - t80 * t68 + t88 * t75 - t113;
t9 = m(4) * t137 - t85 * mrSges(4,2) + t70 * mrSges(4,3) - t105 * t15 + t109 * t12 + t79 * t68 - t88 * t76;
t124 = t106 * t9 + t110 * t11;
t126 = -mrSges(3,1) * t101 + mrSges(3,2) * t98;
t89 = t126 * t142;
t121 = -t103 * mrSges(3,2) + mrSges(3,3) * t146;
t94 = t121 * qJD(1);
t4 = m(3) * (-g(3) * t146 + t128) - t99 * t10 + t124 * t102 + t122 * qJDD(1) + (t103 * t94 - t89 * t152) * qJD(1);
t93 = t122 * qJD(1);
t6 = m(3) * t129 + t102 * t10 + t124 * t99 + (-m(3) * t90 + t126 * qJDD(1) + (-t101 * t94 + t93 * t98) * qJD(1)) * t100;
t8 = m(3) * (-g(3) * t152 + t136) + t110 * t9 - t106 * t11 + t121 * qJDD(1) + (-t103 * t93 + t89 * t146) * qJD(1);
t141 = t103 * t6 + t4 * t146 + t8 * t152;
t2 = m(2) * t127 - t112 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t101 * t8 - t98 * t4;
t1 = m(2) * t134 + qJDD(1) * mrSges(2,1) - t112 * mrSges(2,2) - t100 * t6 + (t101 * t4 + t98 * t8) * t103;
t3 = [-m(1) * g(1) - t107 * t1 + t111 * t2, t2, t8, t9, t12, t14, -t48 * mrSges(7,2) - t72 * t45 + t139; -m(1) * g(2) + t111 * t1 + t107 * t2, t1, t4, t11, t15, t13, -t35 * mrSges(7,3) - t63 * t39 + t140; (-m(1) - m(2)) * g(3) + t141, -m(2) * g(3) + t141, t6, t10, t113, t157, -t34 * mrSges(7,1) - t62 * t42 + t138;];
f_new  = t3;
