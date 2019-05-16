% Calculate vector of cutting forces with Newton-Euler
% S6RRRRPR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-07 23:44
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRPR11_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR11_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR11_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR11_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR11_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR11_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR11_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 23:22:21
% EndTime: 2019-05-07 23:22:55
% DurationCPUTime: 11.98s
% Computational Cost: add. (230213->215), mult. (495510->296), div. (0->0), fcn. (400431->14), ass. (0->116)
t109 = sin(pkin(6));
t115 = sin(qJ(2));
t120 = cos(qJ(2));
t138 = qJD(1) * qJD(2);
t99 = (-qJDD(1) * t120 + t115 * t138) * t109;
t149 = pkin(8) * t109;
t111 = cos(pkin(6));
t148 = t111 * g(3);
t113 = sin(qJ(4));
t118 = cos(qJ(4));
t140 = qJD(1) * t120;
t135 = t109 * t140;
t102 = qJD(3) - t135;
t100 = t102 ^ 2;
t114 = sin(qJ(3));
t119 = cos(qJ(3));
t105 = t111 * qJD(1) + qJD(2);
t103 = t105 ^ 2;
t104 = t111 * qJDD(1) + qJDD(2);
t122 = qJD(1) ^ 2;
t116 = sin(qJ(1));
t121 = cos(qJ(1));
t134 = t116 * g(1) - t121 * g(2);
t94 = qJDD(1) * pkin(1) + t122 * t149 + t134;
t144 = t111 * t94;
t130 = -t121 * g(1) - t116 * g(2);
t95 = -t122 * pkin(1) + qJDD(1) * t149 + t130;
t145 = t115 * t144 + t120 * t95;
t141 = qJD(1) * t109;
t97 = (-pkin(2) * t120 - pkin(9) * t115) * t141;
t58 = -t103 * pkin(2) + t104 * pkin(9) + (-g(3) * t115 + t97 * t140) * t109 + t145;
t98 = (qJDD(1) * t115 + t120 * t138) * t109;
t59 = t99 * pkin(2) - t98 * pkin(9) - t148 + (-t94 + (pkin(2) * t115 - pkin(9) * t120) * t105 * qJD(1)) * t109;
t146 = t114 * t59 + t119 * t58;
t136 = t115 * t141;
t86 = t119 * t105 - t114 * t136;
t87 = t114 * t105 + t119 * t136;
t75 = -t86 * pkin(3) - t87 * pkin(10);
t91 = qJDD(3) + t99;
t35 = -t100 * pkin(3) + t91 * pkin(10) + t86 * t75 + t146;
t142 = t109 * t120;
t129 = -g(3) * t142 - t115 * t95 + t120 * t144;
t57 = -t104 * pkin(2) - t103 * pkin(9) + t97 * t136 - t129;
t72 = -t87 * qJD(3) + t119 * t104 - t114 * t98;
t73 = t86 * qJD(3) + t114 * t104 + t119 * t98;
t38 = (-t102 * t86 - t73) * pkin(10) + (t102 * t87 - t72) * pkin(3) + t57;
t147 = t113 * t38 + t118 * t35;
t143 = t109 * t115;
t132 = -t114 * t58 + t119 * t59;
t34 = -t91 * pkin(3) - t100 * pkin(10) + t87 * t75 - t132;
t78 = t113 * t102 + t118 * t87;
t48 = -t78 * qJD(4) - t113 * t73 + t118 * t91;
t85 = qJD(4) - t86;
t68 = t85 * pkin(4) - t78 * qJ(5);
t77 = t118 * t102 - t113 * t87;
t76 = t77 ^ 2;
t125 = -t48 * pkin(4) - t76 * qJ(5) + t78 * t68 + qJDD(5) + t34;
t112 = sin(qJ(6));
t117 = cos(qJ(6));
t108 = sin(pkin(12));
t110 = cos(pkin(12));
t49 = t77 * qJD(4) + t113 * t91 + t118 * t73;
t40 = -t108 * t49 + t110 * t48;
t41 = t108 * t48 + t110 * t49;
t63 = -t108 * t78 + t110 * t77;
t64 = t108 * t77 + t110 * t78;
t46 = t112 * t63 + t117 * t64;
t27 = -t46 * qJD(6) - t112 * t41 + t117 * t40;
t45 = -t112 * t64 + t117 * t63;
t28 = t45 * qJD(6) + t112 * t40 + t117 * t41;
t83 = qJD(6) + t85;
t42 = -t83 * mrSges(7,2) + t45 * mrSges(7,3);
t43 = t83 * mrSges(7,1) - t46 * mrSges(7,3);
t52 = t85 * pkin(5) - t64 * pkin(11);
t62 = t63 ^ 2;
t127 = t27 * mrSges(7,1) + t45 * t42 - m(7) * (-t40 * pkin(5) - t62 * pkin(11) + t64 * t52 + t125) - t28 * mrSges(7,2) - t46 * t43;
t50 = -t85 * mrSges(6,2) + t63 * mrSges(6,3);
t51 = t85 * mrSges(6,1) - t64 * mrSges(6,3);
t126 = -m(6) * t125 + t40 * mrSges(6,1) - t41 * mrSges(6,2) + t63 * t50 - t64 * t51 + t127;
t67 = -t85 * mrSges(5,2) + t77 * mrSges(5,3);
t69 = t85 * mrSges(5,1) - t78 * mrSges(5,3);
t123 = m(5) * t34 - t48 * mrSges(5,1) + t49 * mrSges(5,2) - t77 * t67 + t78 * t69 - t126;
t74 = -t86 * mrSges(4,1) + t87 * mrSges(4,2);
t79 = -t102 * mrSges(4,2) + t86 * mrSges(4,3);
t16 = m(4) * t132 + t91 * mrSges(4,1) - t73 * mrSges(4,3) + t102 * t79 - t87 * t74 - t123;
t133 = -t113 * t35 + t118 * t38;
t71 = qJDD(4) - t72;
t23 = (t77 * t85 - t49) * qJ(5) + (t77 * t78 + t71) * pkin(4) + t133;
t25 = -t76 * pkin(4) + t48 * qJ(5) - t85 * t68 + t147;
t131 = -0.2e1 * qJD(5) * t64 - t108 * t25 + t110 * t23;
t17 = (t63 * t85 - t41) * pkin(11) + (t63 * t64 + t71) * pkin(5) + t131;
t137 = 0.2e1 * qJD(5) * t63 + t108 * t23 + t110 * t25;
t18 = -t62 * pkin(5) + t40 * pkin(11) - t85 * t52 + t137;
t32 = -t45 * mrSges(7,1) + t46 * mrSges(7,2);
t70 = qJDD(6) + t71;
t14 = m(7) * (-t112 * t18 + t117 * t17) - t28 * mrSges(7,3) + t70 * mrSges(7,1) - t46 * t32 + t83 * t42;
t15 = m(7) * (t112 * t17 + t117 * t18) + t27 * mrSges(7,3) - t70 * mrSges(7,2) + t45 * t32 - t83 * t43;
t47 = -t63 * mrSges(6,1) + t64 * mrSges(6,2);
t12 = m(6) * t131 + t71 * mrSges(6,1) - t41 * mrSges(6,3) + t112 * t15 + t117 * t14 - t64 * t47 + t85 * t50;
t13 = m(6) * t137 - t71 * mrSges(6,2) + t40 * mrSges(6,3) - t112 * t14 + t117 * t15 + t63 * t47 - t85 * t51;
t65 = -t77 * mrSges(5,1) + t78 * mrSges(5,2);
t10 = m(5) * t133 + t71 * mrSges(5,1) - t49 * mrSges(5,3) + t108 * t13 + t110 * t12 - t78 * t65 + t85 * t67;
t11 = m(5) * t147 - t71 * mrSges(5,2) + t48 * mrSges(5,3) - t108 * t12 + t110 * t13 + t77 * t65 - t85 * t69;
t80 = t102 * mrSges(4,1) - t87 * mrSges(4,3);
t9 = m(4) * t146 - t91 * mrSges(4,2) + t72 * mrSges(4,3) - t113 * t10 - t102 * t80 + t118 * t11 + t86 * t74;
t92 = t105 * mrSges(3,1) - mrSges(3,3) * t136;
t96 = (-mrSges(3,1) * t120 + mrSges(3,2) * t115) * t141;
t4 = m(3) * (-g(3) * t143 + t145) - t99 * mrSges(3,3) - t104 * mrSges(3,2) + t96 * t135 - t105 * t92 + t119 * t9 - t114 * t16;
t93 = -t105 * mrSges(3,2) + mrSges(3,3) * t135;
t6 = m(3) * (-t109 * t94 - t148) + t98 * mrSges(3,2) + t99 * mrSges(3,1) + t114 * t9 + t119 * t16 + (t115 * t92 - t120 * t93) * t141;
t124 = m(4) * t57 - t72 * mrSges(4,1) + t73 * mrSges(4,2) + t118 * t10 + t113 * t11 - t86 * t79 + t87 * t80;
t8 = m(3) * t129 + t104 * mrSges(3,1) - t98 * mrSges(3,3) + t105 * t93 - t96 * t136 - t124;
t139 = t111 * t6 + t8 * t142 + t4 * t143;
t2 = m(2) * t130 - t122 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t115 * t8 + t120 * t4;
t1 = m(2) * t134 + qJDD(1) * mrSges(2,1) - t122 * mrSges(2,2) - t109 * t6 + (t115 * t4 + t120 * t8) * t111;
t3 = [-m(1) * g(1) - t116 * t1 + t121 * t2, t2, t4, t9, t11, t13, t15; -m(1) * g(2) + t121 * t1 + t116 * t2, t1, t8, t16, t10, t12, t14; (-m(1) - m(2)) * g(3) + t139, -m(2) * g(3) + t139, t6, t124, t123, -t126, -t127;];
f_new  = t3;
