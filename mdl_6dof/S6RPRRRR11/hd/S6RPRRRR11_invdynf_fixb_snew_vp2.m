% Calculate vector of cutting forces with Newton-Euler
% S6RPRRRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-05-06 05:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRRR11_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR11_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR11_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR11_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_invdynf_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR11_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR11_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR11_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 05:30:28
% EndTime: 2019-05-06 05:31:05
% DurationCPUTime: 20.25s
% Computational Cost: add. (329664->220), mult. (1022471->320), div. (0->0), fcn. (881008->16), ass. (0->128)
t105 = sin(pkin(13));
t107 = sin(pkin(6));
t114 = sin(qJ(3));
t119 = cos(qJ(3));
t108 = cos(pkin(13));
t109 = cos(pkin(7));
t152 = t108 * t109;
t106 = sin(pkin(7));
t110 = cos(pkin(6));
t155 = t106 * t110;
t124 = (-t105 * t114 + t119 * t152) * t107 + t119 * t155;
t86 = t124 * qJD(1);
t151 = t109 * t114;
t154 = t106 * t114;
t126 = t110 * t154 + (t105 * t119 + t108 * t151) * t107;
t87 = t126 * qJD(1);
t74 = -t87 * qJD(3) + t124 * qJDD(1);
t133 = -pkin(9) * t105 * t106 - pkin(2) * t108;
t149 = qJD(1) * t107;
t158 = pkin(9) * qJDD(1);
t130 = qJD(1) * t133 * t149 + t109 * t158;
t144 = qJD(2) * t149;
t153 = t107 * t108;
t121 = qJD(1) ^ 2;
t115 = sin(qJ(1));
t120 = cos(qJ(1));
t143 = t115 * g(1) - t120 * g(2);
t157 = qJ(2) * t107;
t97 = qJDD(1) * pkin(1) + t121 * t157 + t143;
t159 = t110 * t97;
t134 = -g(3) * t153 - 0.2e1 * t105 * t144 + t108 * t159;
t94 = (t107 * t152 + t155) * qJD(1) * pkin(9);
t138 = -t120 * g(1) - t115 * g(2);
t98 = -t121 * pkin(1) + qJDD(1) * t157 + t138;
t57 = (pkin(2) * qJDD(1) + qJD(1) * t94) * t110 + (-t130 * t107 - t98) * t105 + t134;
t146 = t105 * t159 + (0.2e1 * t144 + t98) * t108;
t156 = t105 * t107;
t99 = (-pkin(9) * t109 * t156 + pkin(2) * t110) * qJD(1);
t58 = (-qJD(1) * t99 + t106 * t158) * t110 + (-g(3) * t105 + t130 * t108) * t107 + t146;
t139 = -t110 * g(3) + qJDD(2);
t67 = (-t97 + t133 * qJDD(1) + (t105 * t99 - t108 * t94) * qJD(1)) * t107 + t139;
t162 = -t114 * t58 + (t106 * t67 + t109 * t57) * t119;
t112 = sin(qJ(5));
t117 = cos(qJ(5));
t113 = sin(qJ(4));
t118 = cos(qJ(4));
t147 = t119 * t58 + t57 * t151 + t67 * t154;
t73 = -t86 * pkin(3) - t87 * pkin(10);
t128 = -t106 * t153 + t109 * t110;
t95 = t128 * qJD(1) + qJD(3);
t91 = t95 ^ 2;
t92 = t128 * qJDD(1) + qJDD(3);
t34 = -t91 * pkin(3) + t92 * pkin(10) + t86 * t73 + t147;
t142 = -t106 * t57 + t109 * t67;
t75 = t86 * qJD(3) + t126 * qJDD(1);
t36 = (-t86 * t95 - t75) * pkin(10) + (t87 * t95 - t74) * pkin(3) + t142;
t160 = t113 * t36 + t118 * t34;
t79 = -t113 * t87 + t118 * t95;
t80 = t113 * t95 + t118 * t87;
t60 = -t79 * pkin(4) - t80 * pkin(11);
t71 = qJDD(4) - t74;
t85 = qJD(4) - t86;
t84 = t85 ^ 2;
t24 = -t84 * pkin(4) + t71 * pkin(11) + t79 * t60 + t160;
t33 = -t92 * pkin(3) - t91 * pkin(10) + t87 * t73 - t162;
t52 = -t80 * qJD(4) - t113 * t75 + t118 * t92;
t53 = t79 * qJD(4) + t113 * t92 + t118 * t75;
t27 = (-t79 * t85 - t53) * pkin(11) + (t80 * t85 - t52) * pkin(4) + t33;
t161 = t112 * t27 + t117 * t24;
t111 = sin(qJ(6));
t116 = cos(qJ(6));
t141 = -t112 * t24 + t117 * t27;
t65 = -t112 * t80 + t117 * t85;
t40 = t65 * qJD(5) + t112 * t71 + t117 * t53;
t51 = qJDD(5) - t52;
t66 = t112 * t85 + t117 * t80;
t78 = qJD(5) - t79;
t18 = (t65 * t78 - t40) * pkin(12) + (t65 * t66 + t51) * pkin(5) + t141;
t39 = -t66 * qJD(5) - t112 * t53 + t117 * t71;
t49 = t78 * pkin(5) - t66 * pkin(12);
t64 = t65 ^ 2;
t19 = -t64 * pkin(5) + t39 * pkin(12) - t78 * t49 + t161;
t43 = -t111 * t66 + t116 * t65;
t30 = t43 * qJD(6) + t111 * t39 + t116 * t40;
t44 = t111 * t65 + t116 * t66;
t38 = -t43 * mrSges(7,1) + t44 * mrSges(7,2);
t76 = qJD(6) + t78;
t41 = -t76 * mrSges(7,2) + t43 * mrSges(7,3);
t50 = qJDD(6) + t51;
t16 = m(7) * (-t111 * t19 + t116 * t18) - t30 * mrSges(7,3) + t50 * mrSges(7,1) - t44 * t38 + t76 * t41;
t29 = -t44 * qJD(6) - t111 * t40 + t116 * t39;
t42 = t76 * mrSges(7,1) - t44 * mrSges(7,3);
t17 = m(7) * (t111 * t18 + t116 * t19) + t29 * mrSges(7,3) - t50 * mrSges(7,2) + t43 * t38 - t76 * t42;
t45 = -t65 * mrSges(6,1) + t66 * mrSges(6,2);
t47 = -t78 * mrSges(6,2) + t65 * mrSges(6,3);
t13 = m(6) * t141 + t51 * mrSges(6,1) - t40 * mrSges(6,3) + t111 * t17 + t116 * t16 - t66 * t45 + t78 * t47;
t48 = t78 * mrSges(6,1) - t66 * mrSges(6,3);
t14 = m(6) * t161 - t51 * mrSges(6,2) + t39 * mrSges(6,3) - t111 * t16 + t116 * t17 + t65 * t45 - t78 * t48;
t59 = -t79 * mrSges(5,1) + t80 * mrSges(5,2);
t69 = t85 * mrSges(5,1) - t80 * mrSges(5,3);
t12 = m(5) * t160 - t71 * mrSges(5,2) + t52 * mrSges(5,3) - t112 * t13 + t117 * t14 + t79 * t59 - t85 * t69;
t140 = -t113 * t34 + t118 * t36;
t23 = -t71 * pkin(4) - t84 * pkin(11) + t80 * t60 - t140;
t127 = t29 * mrSges(7,1) + t43 * t41 - m(7) * (-t39 * pkin(5) - t64 * pkin(12) + t66 * t49 + t23) - t30 * mrSges(7,2) - t44 * t42;
t122 = m(6) * t23 - t39 * mrSges(6,1) + t40 * mrSges(6,2) - t65 * t47 + t66 * t48 - t127;
t68 = -t85 * mrSges(5,2) + t79 * mrSges(5,3);
t15 = m(5) * t140 + t71 * mrSges(5,1) - t53 * mrSges(5,3) - t80 * t59 + t85 * t68 - t122;
t81 = -t95 * mrSges(4,2) + t86 * mrSges(4,3);
t82 = t95 * mrSges(4,1) - t87 * mrSges(4,3);
t10 = m(4) * t142 - t74 * mrSges(4,1) + t75 * mrSges(4,2) + t113 * t12 + t118 * t15 - t86 * t81 + t87 * t82;
t131 = -t110 * mrSges(3,2) + mrSges(3,3) * t153;
t101 = t131 * qJD(1);
t132 = t110 * mrSges(3,1) - mrSges(3,3) * t156;
t123 = m(5) * t33 - t52 * mrSges(5,1) + t53 * mrSges(5,2) + t112 * t14 + t117 * t13 - t79 * t68 + t80 * t69;
t72 = -t86 * mrSges(4,1) + t87 * mrSges(4,2);
t11 = m(4) * t162 + t92 * mrSges(4,1) - t75 * mrSges(4,3) - t87 * t72 + t95 * t81 - t123;
t9 = m(4) * t147 - t92 * mrSges(4,2) + t74 * mrSges(4,3) - t113 * t15 + t118 * t12 + t86 * t72 - t95 * t82;
t137 = t119 * t11 + t114 * t9;
t136 = -mrSges(3,1) * t108 + mrSges(3,2) * t105;
t96 = t136 * t149;
t4 = m(3) * (-t105 * t98 + t134) - t106 * t10 + t137 * t109 + t132 * qJDD(1) + (t110 * t101 - t96 * t156) * qJD(1);
t100 = t132 * qJD(1);
t6 = m(3) * t139 + t109 * t10 + t137 * t106 + (-m(3) * t97 + t136 * qJDD(1) + (t100 * t105 - t101 * t108) * qJD(1)) * t107;
t8 = m(3) * (-g(3) * t156 + t146) + t119 * t9 - t114 * t11 + t131 * qJDD(1) + (-t110 * t100 + t96 * t153) * qJD(1);
t148 = t110 * t6 + t4 * t153 + t8 * t156;
t2 = m(2) * t138 - t121 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t105 * t4 + t108 * t8;
t1 = m(2) * t143 + qJDD(1) * mrSges(2,1) - t121 * mrSges(2,2) - t107 * t6 + (t105 * t8 + t108 * t4) * t110;
t3 = [-m(1) * g(1) - t115 * t1 + t120 * t2, t2, t8, t9, t12, t14, t17; -m(1) * g(2) + t120 * t1 + t115 * t2, t1, t4, t11, t15, t13, t16; (-m(1) - m(2)) * g(3) + t148, -m(2) * g(3) + t148, t6, t10, t123, t122, -t127;];
f_new  = t3;
