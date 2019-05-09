% Calculate vector of cutting forces with Newton-Euler
% S6RRRRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-08 02:41
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRPR14_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR14_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR14_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR14_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_invdynf_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR14_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR14_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR14_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 02:07:09
% EndTime: 2019-05-08 02:08:11
% DurationCPUTime: 23.41s
% Computational Cost: add. (432459->228), mult. (1072537->326), div. (0->0), fcn. (908944->16), ass. (0->129)
t107 = sin(pkin(7));
t114 = sin(qJ(3));
t118 = cos(qJ(3));
t110 = cos(pkin(7));
t108 = sin(pkin(6));
t115 = sin(qJ(2));
t119 = cos(qJ(2));
t141 = qJD(1) * qJD(2);
t101 = (qJDD(1) * t115 + t119 * t141) * t108;
t111 = cos(pkin(6));
t103 = qJDD(1) * t111 + qJDD(2);
t104 = qJD(1) * t111 + qJD(2);
t121 = qJD(1) ^ 2;
t116 = sin(qJ(1));
t120 = cos(qJ(1));
t135 = t116 * g(1) - g(2) * t120;
t158 = pkin(9) * t108;
t98 = qJDD(1) * pkin(1) + t121 * t158 + t135;
t151 = t111 * t98;
t130 = -g(1) * t120 - g(2) * t116;
t99 = -pkin(1) * t121 + qJDD(1) * t158 + t130;
t133 = -t115 * t99 + t119 * t151;
t144 = qJD(1) * t115;
t156 = pkin(10) * t110;
t143 = qJD(1) * t119;
t136 = t108 * t143;
t150 = t104 * t107;
t91 = (t110 * t136 + t150) * pkin(10);
t145 = qJD(1) * t108;
t157 = pkin(10) * t107;
t95 = (-pkin(2) * t119 - t115 * t157) * t145;
t56 = -t101 * t156 + t103 * pkin(2) + t104 * t91 + (-g(3) * t119 - t95 * t144) * t108 + t133;
t152 = t110 * t56;
t102 = (qJDD(1) * t119 - t115 * t141) * t108;
t126 = t102 * t110 + t103 * t107;
t153 = t115 * t151 + t119 * t99;
t137 = t108 * t144;
t94 = pkin(2) * t104 - t137 * t156;
t57 = -t104 * t94 + (-g(3) * t115 + t95 * t143) * t108 + t126 * pkin(10) + t153;
t155 = t111 * g(3);
t62 = -t101 * t157 - t102 * pkin(2) - t155 + (-t98 + (t115 * t94 - t119 * t91) * qJD(1)) * t108;
t160 = -t114 * t57 + (t107 * t62 + t152) * t118;
t146 = t110 * t119;
t149 = t107 * t114;
t85 = t104 * t149 + (t114 * t146 + t115 * t118) * t145;
t72 = -t85 * qJD(3) - t114 * t101 + t126 * t118;
t84 = (-t114 * t115 + t118 * t146) * t145 + t118 * t150;
t159 = cos(qJ(4));
t113 = sin(qJ(4));
t139 = t114 * t152 + t118 * t57 + t62 * t149;
t75 = -pkin(3) * t84 - pkin(11) * t85;
t86 = -t102 * t107 + t103 * t110 + qJDD(3);
t92 = t104 * t110 - t107 * t136 + qJD(3);
t90 = t92 ^ 2;
t34 = -pkin(3) * t90 + pkin(11) * t86 + t75 * t84 + t139;
t134 = -t107 * t56 + t110 * t62;
t73 = t84 * qJD(3) + t118 * t101 + t126 * t114;
t36 = (-t84 * t92 - t73) * pkin(11) + (t85 * t92 - t72) * pkin(3) + t134;
t154 = t113 * t36 + t159 * t34;
t148 = t108 * t115;
t147 = t108 * t119;
t106 = sin(pkin(13));
t109 = cos(pkin(13));
t112 = sin(qJ(6));
t117 = cos(qJ(6));
t77 = t113 * t85 - t159 * t92;
t78 = t113 * t92 + t159 * t85;
t58 = pkin(4) * t77 - qJ(5) * t78;
t71 = qJDD(4) - t72;
t82 = qJD(4) - t84;
t81 = t82 ^ 2;
t24 = -pkin(4) * t81 + qJ(5) * t71 - t58 * t77 + t154;
t33 = -t86 * pkin(3) - t90 * pkin(11) + t85 * t75 - t160;
t47 = qJD(4) * t78 + t113 * t73 - t159 * t86;
t48 = -t77 * qJD(4) + t113 * t86 + t159 * t73;
t27 = (t77 * t82 - t48) * qJ(5) + (t78 * t82 + t47) * pkin(4) + t33;
t68 = t106 * t82 + t109 * t78;
t132 = -0.2e1 * qJD(5) * t68 - t106 * t24 + t109 * t27;
t40 = t106 * t71 + t109 * t48;
t67 = -t106 * t78 + t109 * t82;
t18 = (t67 * t77 - t40) * pkin(12) + (t67 * t68 + t47) * pkin(5) + t132;
t140 = 0.2e1 * qJD(5) * t67 + t106 * t27 + t109 * t24;
t39 = -t106 * t48 + t109 * t71;
t52 = pkin(5) * t77 - pkin(12) * t68;
t66 = t67 ^ 2;
t19 = -pkin(5) * t66 + pkin(12) * t39 - t52 * t77 + t140;
t43 = -t112 * t68 + t117 * t67;
t30 = qJD(6) * t43 + t112 * t39 + t117 * t40;
t44 = t112 * t67 + t117 * t68;
t38 = -mrSges(7,1) * t43 + mrSges(7,2) * t44;
t76 = qJD(6) + t77;
t41 = -mrSges(7,2) * t76 + mrSges(7,3) * t43;
t46 = qJDD(6) + t47;
t16 = m(7) * (-t112 * t19 + t117 * t18) - t30 * mrSges(7,3) + t46 * mrSges(7,1) - t44 * t38 + t76 * t41;
t29 = -qJD(6) * t44 - t112 * t40 + t117 * t39;
t42 = mrSges(7,1) * t76 - mrSges(7,3) * t44;
t17 = m(7) * (t112 * t18 + t117 * t19) + t29 * mrSges(7,3) - t46 * mrSges(7,2) + t43 * t38 - t76 * t42;
t45 = -mrSges(6,1) * t67 + mrSges(6,2) * t68;
t50 = -mrSges(6,2) * t77 + mrSges(6,3) * t67;
t13 = m(6) * t132 + t47 * mrSges(6,1) - t40 * mrSges(6,3) + t112 * t17 + t117 * t16 - t68 * t45 + t77 * t50;
t51 = mrSges(6,1) * t77 - mrSges(6,3) * t68;
t14 = m(6) * t140 - t47 * mrSges(6,2) + t39 * mrSges(6,3) - t112 * t16 + t117 * t17 + t67 * t45 - t77 * t51;
t59 = mrSges(5,1) * t77 + mrSges(5,2) * t78;
t70 = mrSges(5,1) * t82 - mrSges(5,3) * t78;
t12 = m(5) * t154 - t71 * mrSges(5,2) - t47 * mrSges(5,3) - t106 * t13 + t109 * t14 - t77 * t59 - t82 * t70;
t131 = -t113 * t34 + t159 * t36;
t23 = -t71 * pkin(4) - t81 * qJ(5) + t78 * t58 + qJDD(5) - t131;
t124 = t29 * mrSges(7,1) + t43 * t41 - m(7) * (-t39 * pkin(5) - t66 * pkin(12) + t68 * t52 + t23) - t30 * mrSges(7,2) - t44 * t42;
t122 = m(6) * t23 - t39 * mrSges(6,1) + t40 * mrSges(6,2) - t67 * t50 + t68 * t51 - t124;
t69 = -mrSges(5,2) * t82 - mrSges(5,3) * t77;
t15 = m(5) * t131 + t71 * mrSges(5,1) - t48 * mrSges(5,3) - t78 * t59 + t82 * t69 - t122;
t79 = -mrSges(4,2) * t92 + mrSges(4,3) * t84;
t80 = mrSges(4,1) * t92 - mrSges(4,3) * t85;
t10 = m(4) * t134 - t72 * mrSges(4,1) + t73 * mrSges(4,2) + t113 * t12 + t159 * t15 - t84 * t79 + t85 * t80;
t100 = (-mrSges(3,1) * t119 + mrSges(3,2) * t115) * t145;
t123 = m(5) * t33 + t47 * mrSges(5,1) + t48 * mrSges(5,2) + t106 * t14 + t109 * t13 + t77 * t69 + t78 * t70;
t74 = -mrSges(4,1) * t84 + mrSges(4,2) * t85;
t11 = m(4) * t160 + t86 * mrSges(4,1) - t73 * mrSges(4,3) - t85 * t74 + t92 * t79 - t123;
t9 = m(4) * t139 - t86 * mrSges(4,2) + t72 * mrSges(4,3) - t113 * t15 + t159 * t12 + t84 * t74 - t92 * t80;
t129 = t11 * t118 + t114 * t9;
t97 = -mrSges(3,2) * t104 + mrSges(3,3) * t136;
t4 = m(3) * (-g(3) * t147 + t133) - t101 * mrSges(3,3) + t103 * mrSges(3,1) - t100 * t137 + t104 * t97 - t107 * t10 + t129 * t110;
t96 = mrSges(3,1) * t104 - mrSges(3,3) * t137;
t6 = m(3) * (-t108 * t98 - t155) + t101 * mrSges(3,2) - t102 * mrSges(3,1) + t110 * t10 + t129 * t107 + (t115 * t96 - t119 * t97) * t145;
t8 = m(3) * (-g(3) * t148 + t153) + t102 * mrSges(3,3) - t103 * mrSges(3,2) + t100 * t136 - t104 * t96 + t118 * t9 - t114 * t11;
t142 = t111 * t6 + t4 * t147 + t8 * t148;
t2 = m(2) * t130 - t121 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t115 * t4 + t119 * t8;
t1 = m(2) * t135 + qJDD(1) * mrSges(2,1) - t121 * mrSges(2,2) - t108 * t6 + (t115 * t8 + t119 * t4) * t111;
t3 = [-m(1) * g(1) - t1 * t116 + t120 * t2, t2, t8, t9, t12, t14, t17; -m(1) * g(2) + t1 * t120 + t116 * t2, t1, t4, t11, t15, t13, t16; (-m(1) - m(2)) * g(3) + t142, -m(2) * g(3) + t142, t6, t10, t123, t122, -t124;];
f_new  = t3;
