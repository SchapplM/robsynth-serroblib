% Calculate vector of cutting forces with Newton-Euler
% S6RPRRPR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-05-06 00:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRPR11_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR11_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR11_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR11_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_invdynf_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR11_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR11_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR11_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 00:05:25
% EndTime: 2019-05-06 00:05:58
% DurationCPUTime: 19.91s
% Computational Cost: add. (318553->219), mult. (997305->322), div. (0->0), fcn. (856010->16), ass. (0->126)
t105 = sin(pkin(12));
t107 = sin(pkin(6));
t114 = sin(qJ(3));
t117 = cos(qJ(3));
t109 = cos(pkin(12));
t110 = cos(pkin(7));
t151 = t109 * t110;
t106 = sin(pkin(7));
t111 = cos(pkin(6));
t154 = t106 * t111;
t122 = (-t105 * t114 + t117 * t151) * t107 + t117 * t154;
t83 = t122 * qJD(1);
t150 = t110 * t114;
t153 = t106 * t114;
t124 = t111 * t153 + (t105 * t117 + t109 * t150) * t107;
t84 = t124 * qJD(1);
t74 = -t84 * qJD(3) + t122 * qJDD(1);
t131 = -pkin(9) * t105 * t106 - pkin(2) * t109;
t148 = qJD(1) * t107;
t157 = pkin(9) * qJDD(1);
t128 = qJD(1) * t131 * t148 + t110 * t157;
t142 = qJD(2) * t148;
t152 = t107 * t109;
t119 = qJD(1) ^ 2;
t115 = sin(qJ(1));
t118 = cos(qJ(1));
t141 = t115 * g(1) - t118 * g(2);
t156 = qJ(2) * t107;
t96 = qJDD(1) * pkin(1) + t119 * t156 + t141;
t158 = t111 * t96;
t132 = -g(3) * t152 - 0.2e1 * t105 * t142 + t109 * t158;
t93 = (t107 * t151 + t154) * qJD(1) * pkin(9);
t136 = -t118 * g(1) - t115 * g(2);
t97 = -t119 * pkin(1) + qJDD(1) * t156 + t136;
t56 = (pkin(2) * qJDD(1) + qJD(1) * t93) * t111 + (-t107 * t128 - t97) * t105 + t132;
t144 = t105 * t158 + (0.2e1 * t142 + t97) * t109;
t155 = t105 * t107;
t98 = (-pkin(9) * t110 * t155 + pkin(2) * t111) * qJD(1);
t57 = (-qJD(1) * t98 + t106 * t157) * t111 + (-g(3) * t105 + t128 * t109) * t107 + t144;
t139 = -t111 * g(3) + qJDD(2);
t68 = (-t96 + t131 * qJDD(1) + (t105 * t98 - t109 * t93) * qJD(1)) * t107 + t139;
t161 = -t114 * t57 + (t106 * t68 + t110 * t56) * t117;
t160 = cos(qJ(4));
t113 = sin(qJ(4));
t145 = t117 * t57 + t56 * t150 + t68 * t153;
t73 = -t83 * pkin(3) - t84 * pkin(10);
t126 = -t106 * t152 + t110 * t111;
t94 = t126 * qJD(1) + qJD(3);
t90 = t94 ^ 2;
t91 = t126 * qJDD(1) + qJDD(3);
t34 = -t90 * pkin(3) + t91 * pkin(10) + t83 * t73 + t145;
t140 = -t106 * t56 + t110 * t68;
t75 = t83 * qJD(3) + t124 * qJDD(1);
t37 = (-t83 * t94 - t75) * pkin(10) + (t84 * t94 - t74) * pkin(3) + t140;
t159 = t113 * t37 + t160 * t34;
t104 = sin(pkin(13));
t108 = cos(pkin(13));
t112 = sin(qJ(6));
t116 = cos(qJ(6));
t77 = t113 * t84 - t160 * t94;
t78 = t113 * t94 + t160 * t84;
t58 = t77 * pkin(4) - t78 * qJ(5);
t71 = qJDD(4) - t74;
t82 = qJD(4) - t83;
t81 = t82 ^ 2;
t24 = -t81 * pkin(4) + t71 * qJ(5) - t77 * t58 + t159;
t33 = -t91 * pkin(3) - t90 * pkin(10) + t84 * t73 - t161;
t51 = t78 * qJD(4) + t113 * t75 - t160 * t91;
t52 = -t77 * qJD(4) + t113 * t91 + t160 * t75;
t27 = (t77 * t82 - t52) * qJ(5) + (t78 * t82 + t51) * pkin(4) + t33;
t67 = t104 * t82 + t108 * t78;
t138 = -0.2e1 * qJD(5) * t67 - t104 * t24 + t108 * t27;
t42 = t104 * t71 + t108 * t52;
t66 = -t104 * t78 + t108 * t82;
t18 = (t66 * t77 - t42) * pkin(11) + (t66 * t67 + t51) * pkin(5) + t138;
t146 = 0.2e1 * qJD(5) * t66 + t104 * t27 + t108 * t24;
t41 = -t104 * t52 + t108 * t71;
t49 = t77 * pkin(5) - t67 * pkin(11);
t65 = t66 ^ 2;
t19 = -t65 * pkin(5) + t41 * pkin(11) - t77 * t49 + t146;
t43 = -t112 * t67 + t116 * t66;
t30 = t43 * qJD(6) + t112 * t41 + t116 * t42;
t44 = t112 * t66 + t116 * t67;
t38 = -t43 * mrSges(7,1) + t44 * mrSges(7,2);
t76 = qJD(6) + t77;
t39 = -t76 * mrSges(7,2) + t43 * mrSges(7,3);
t50 = qJDD(6) + t51;
t16 = m(7) * (-t112 * t19 + t116 * t18) - t30 * mrSges(7,3) + t50 * mrSges(7,1) - t44 * t38 + t76 * t39;
t29 = -t44 * qJD(6) - t112 * t42 + t116 * t41;
t40 = t76 * mrSges(7,1) - t44 * mrSges(7,3);
t17 = m(7) * (t112 * t18 + t116 * t19) + t29 * mrSges(7,3) - t50 * mrSges(7,2) + t43 * t38 - t76 * t40;
t45 = -t66 * mrSges(6,1) + t67 * mrSges(6,2);
t47 = -t77 * mrSges(6,2) + t66 * mrSges(6,3);
t13 = m(6) * t138 + t51 * mrSges(6,1) - t42 * mrSges(6,3) + t112 * t17 + t116 * t16 - t67 * t45 + t77 * t47;
t48 = t77 * mrSges(6,1) - t67 * mrSges(6,3);
t14 = m(6) * t146 - t51 * mrSges(6,2) + t41 * mrSges(6,3) - t112 * t16 + t116 * t17 + t66 * t45 - t77 * t48;
t59 = t77 * mrSges(5,1) + t78 * mrSges(5,2);
t70 = t82 * mrSges(5,1) - t78 * mrSges(5,3);
t12 = m(5) * t159 - t71 * mrSges(5,2) - t51 * mrSges(5,3) - t104 * t13 + t108 * t14 - t77 * t59 - t82 * t70;
t137 = -t113 * t34 + t160 * t37;
t23 = -t71 * pkin(4) - t81 * qJ(5) + t78 * t58 + qJDD(5) - t137;
t125 = t29 * mrSges(7,1) + t43 * t39 - m(7) * (-t41 * pkin(5) - t65 * pkin(11) + t67 * t49 + t23) - t30 * mrSges(7,2) - t44 * t40;
t120 = m(6) * t23 - t41 * mrSges(6,1) + t42 * mrSges(6,2) - t66 * t47 + t67 * t48 - t125;
t69 = -t82 * mrSges(5,2) - t77 * mrSges(5,3);
t15 = m(5) * t137 + t71 * mrSges(5,1) - t52 * mrSges(5,3) - t78 * t59 + t82 * t69 - t120;
t79 = -t94 * mrSges(4,2) + t83 * mrSges(4,3);
t80 = t94 * mrSges(4,1) - t84 * mrSges(4,3);
t10 = m(4) * t140 - t74 * mrSges(4,1) + t75 * mrSges(4,2) + t113 * t12 + t160 * t15 - t83 * t79 + t84 * t80;
t129 = -t111 * mrSges(3,2) + mrSges(3,3) * t152;
t100 = t129 * qJD(1);
t130 = t111 * mrSges(3,1) - mrSges(3,3) * t155;
t121 = m(5) * t33 + t51 * mrSges(5,1) + t52 * mrSges(5,2) + t104 * t14 + t108 * t13 + t77 * t69 + t78 * t70;
t72 = -t83 * mrSges(4,1) + t84 * mrSges(4,2);
t11 = m(4) * t161 + t91 * mrSges(4,1) - t75 * mrSges(4,3) - t84 * t72 + t94 * t79 - t121;
t9 = m(4) * t145 - t91 * mrSges(4,2) + t74 * mrSges(4,3) - t113 * t15 + t160 * t12 + t83 * t72 - t94 * t80;
t135 = t117 * t11 + t114 * t9;
t134 = -mrSges(3,1) * t109 + mrSges(3,2) * t105;
t95 = t134 * t148;
t4 = m(3) * (-t105 * t97 + t132) - t106 * t10 + t135 * t110 + t130 * qJDD(1) + (t111 * t100 - t95 * t155) * qJD(1);
t99 = t130 * qJD(1);
t6 = m(3) * t139 + t110 * t10 + t135 * t106 + (-m(3) * t96 + t134 * qJDD(1) + (-t100 * t109 + t105 * t99) * qJD(1)) * t107;
t8 = m(3) * (-g(3) * t155 + t144) + t117 * t9 - t114 * t11 + t129 * qJDD(1) + (-t111 * t99 + t95 * t152) * qJD(1);
t147 = t111 * t6 + t4 * t152 + t8 * t155;
t2 = m(2) * t136 - t119 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t105 * t4 + t109 * t8;
t1 = m(2) * t141 + qJDD(1) * mrSges(2,1) - t119 * mrSges(2,2) - t107 * t6 + (t105 * t8 + t109 * t4) * t111;
t3 = [-m(1) * g(1) - t115 * t1 + t118 * t2, t2, t8, t9, t12, t14, t17; -m(1) * g(2) + t118 * t1 + t115 * t2, t1, t4, t11, t15, t13, t16; (-m(1) - m(2)) * g(3) + t147, -m(2) * g(3) + t147, t6, t10, t121, t120, -t125;];
f_new  = t3;
