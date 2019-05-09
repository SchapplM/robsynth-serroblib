% Calculate vector of cutting forces with Newton-Euler
% S6RRRRRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 12:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRRR7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR7_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR7_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR7_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR7_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR7_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 12:16:25
% EndTime: 2019-05-08 12:17:01
% DurationCPUTime: 12.50s
% Computational Cost: add. (237903->216), mult. (507031->294), div. (0->0), fcn. (411931->14), ass. (0->118)
t108 = sin(pkin(6));
t114 = sin(qJ(2));
t120 = cos(qJ(2));
t137 = qJD(1) * qJD(2);
t99 = (-qJDD(1) * t120 + t114 * t137) * t108;
t149 = pkin(8) * t108;
t109 = cos(pkin(6));
t148 = t109 * g(3);
t111 = sin(qJ(5));
t117 = cos(qJ(5));
t112 = sin(qJ(4));
t118 = cos(qJ(4));
t139 = qJD(1) * t120;
t135 = t108 * t139;
t102 = qJD(3) - t135;
t100 = t102 ^ 2;
t113 = sin(qJ(3));
t119 = cos(qJ(3));
t105 = t109 * qJD(1) + qJD(2);
t103 = t105 ^ 2;
t104 = t109 * qJDD(1) + qJDD(2);
t122 = qJD(1) ^ 2;
t115 = sin(qJ(1));
t121 = cos(qJ(1));
t134 = t115 * g(1) - t121 * g(2);
t94 = qJDD(1) * pkin(1) + t122 * t149 + t134;
t143 = t109 * t94;
t130 = -t121 * g(1) - t115 * g(2);
t95 = -t122 * pkin(1) + qJDD(1) * t149 + t130;
t144 = t114 * t143 + t120 * t95;
t140 = qJD(1) * t108;
t97 = (-pkin(2) * t120 - pkin(9) * t114) * t140;
t58 = -t103 * pkin(2) + t104 * pkin(9) + (-g(3) * t114 + t139 * t97) * t108 + t144;
t98 = (qJDD(1) * t114 + t120 * t137) * t108;
t59 = t99 * pkin(2) - t98 * pkin(9) - t148 + (-t94 + (pkin(2) * t114 - pkin(9) * t120) * t105 * qJD(1)) * t108;
t145 = t113 * t59 + t119 * t58;
t136 = t114 * t140;
t86 = t119 * t105 - t113 * t136;
t87 = t113 * t105 + t119 * t136;
t74 = -t86 * pkin(3) - t87 * pkin(10);
t91 = qJDD(3) + t99;
t38 = -t100 * pkin(3) + t91 * pkin(10) + t86 * t74 + t145;
t141 = t108 * t120;
t129 = -g(3) * t141 - t114 * t95 + t120 * t143;
t57 = -t104 * pkin(2) - t103 * pkin(9) + t97 * t136 - t129;
t71 = -t87 * qJD(3) + t119 * t104 - t113 * t98;
t72 = t86 * qJD(3) + t113 * t104 + t119 * t98;
t41 = (-t102 * t86 - t72) * pkin(10) + (t102 * t87 - t71) * pkin(3) + t57;
t132 = -t112 * t38 + t118 * t41;
t76 = t118 * t102 - t112 * t87;
t49 = t76 * qJD(4) + t112 * t91 + t118 * t72;
t70 = qJDD(4) - t71;
t77 = t112 * t102 + t118 * t87;
t85 = qJD(4) - t86;
t23 = (t76 * t85 - t49) * pkin(11) + (t76 * t77 + t70) * pkin(4) + t132;
t146 = t112 * t41 + t118 * t38;
t48 = -t77 * qJD(4) - t112 * t72 + t118 * t91;
t67 = t85 * pkin(4) - t77 * pkin(11);
t75 = t76 ^ 2;
t28 = -t75 * pkin(4) + t48 * pkin(11) - t85 * t67 + t146;
t147 = t111 * t23 + t117 * t28;
t142 = t108 * t114;
t131 = -t113 * t58 + t119 * t59;
t37 = -t91 * pkin(3) - t100 * pkin(10) + t87 * t74 - t131;
t126 = -t48 * pkin(4) - t75 * pkin(11) + t77 * t67 + t37;
t110 = sin(qJ(6));
t116 = cos(qJ(6));
t62 = t111 * t76 + t117 * t77;
t34 = -t62 * qJD(5) - t111 * t49 + t117 * t48;
t61 = -t111 * t77 + t117 * t76;
t35 = t61 * qJD(5) + t111 * t48 + t117 * t49;
t46 = t110 * t61 + t116 * t62;
t25 = -t46 * qJD(6) - t110 * t35 + t116 * t34;
t45 = -t110 * t62 + t116 * t61;
t26 = t45 * qJD(6) + t110 * t34 + t116 * t35;
t83 = qJD(5) + t85;
t80 = qJD(6) + t83;
t42 = -t80 * mrSges(7,2) + t45 * mrSges(7,3);
t43 = t80 * mrSges(7,1) - t46 * mrSges(7,3);
t52 = t83 * pkin(5) - t62 * pkin(12);
t60 = t61 ^ 2;
t127 = t25 * mrSges(7,1) + t45 * t42 - m(7) * (-t34 * pkin(5) - t60 * pkin(12) + t62 * t52 + t126) - t26 * mrSges(7,2) - t46 * t43;
t50 = -t83 * mrSges(6,2) + t61 * mrSges(6,3);
t51 = t83 * mrSges(6,1) - t62 * mrSges(6,3);
t125 = -m(6) * t126 + t34 * mrSges(6,1) - t35 * mrSges(6,2) + t61 * t50 - t62 * t51 + t127;
t65 = -t85 * mrSges(5,2) + t76 * mrSges(5,3);
t66 = t85 * mrSges(5,1) - t77 * mrSges(5,3);
t123 = m(5) * t37 - t48 * mrSges(5,1) + t49 * mrSges(5,2) - t76 * t65 + t77 * t66 - t125;
t73 = -t86 * mrSges(4,1) + t87 * mrSges(4,2);
t78 = -t102 * mrSges(4,2) + t86 * mrSges(4,3);
t14 = m(4) * t131 + t91 * mrSges(4,1) - t72 * mrSges(4,3) + t102 * t78 - t87 * t73 - t123;
t133 = -t111 * t28 + t117 * t23;
t69 = qJDD(5) + t70;
t17 = (t61 * t83 - t35) * pkin(12) + (t61 * t62 + t69) * pkin(5) + t133;
t18 = -t60 * pkin(5) + t34 * pkin(12) - t83 * t52 + t147;
t32 = -t45 * mrSges(7,1) + t46 * mrSges(7,2);
t68 = qJDD(6) + t69;
t15 = m(7) * (-t110 * t18 + t116 * t17) - t26 * mrSges(7,3) + t68 * mrSges(7,1) - t46 * t32 + t80 * t42;
t16 = m(7) * (t110 * t17 + t116 * t18) + t25 * mrSges(7,3) - t68 * mrSges(7,2) + t45 * t32 - t80 * t43;
t47 = -t61 * mrSges(6,1) + t62 * mrSges(6,2);
t12 = m(6) * t133 + t69 * mrSges(6,1) - t35 * mrSges(6,3) + t110 * t16 + t116 * t15 - t62 * t47 + t83 * t50;
t13 = m(6) * t147 - t69 * mrSges(6,2) + t34 * mrSges(6,3) - t110 * t15 + t116 * t16 + t61 * t47 - t83 * t51;
t63 = -t76 * mrSges(5,1) + t77 * mrSges(5,2);
t10 = m(5) * t132 + t70 * mrSges(5,1) - t49 * mrSges(5,3) + t111 * t13 + t117 * t12 - t77 * t63 + t85 * t65;
t11 = m(5) * t146 - t70 * mrSges(5,2) + t48 * mrSges(5,3) - t111 * t12 + t117 * t13 + t76 * t63 - t85 * t66;
t79 = t102 * mrSges(4,1) - t87 * mrSges(4,3);
t9 = m(4) * t145 - t91 * mrSges(4,2) + t71 * mrSges(4,3) - t112 * t10 - t102 * t79 + t118 * t11 + t86 * t73;
t92 = t105 * mrSges(3,1) - mrSges(3,3) * t136;
t96 = (-mrSges(3,1) * t120 + mrSges(3,2) * t114) * t140;
t4 = m(3) * (-g(3) * t142 + t144) - t99 * mrSges(3,3) - t104 * mrSges(3,2) + t96 * t135 - t105 * t92 + t119 * t9 - t113 * t14;
t93 = -t105 * mrSges(3,2) + mrSges(3,3) * t135;
t6 = m(3) * (-t108 * t94 - t148) + t98 * mrSges(3,2) + t99 * mrSges(3,1) + t113 * t9 + t119 * t14 + (t114 * t92 - t120 * t93) * t140;
t124 = m(4) * t57 - t71 * mrSges(4,1) + t72 * mrSges(4,2) + t118 * t10 + t112 * t11 - t86 * t78 + t87 * t79;
t8 = m(3) * t129 + t104 * mrSges(3,1) - t98 * mrSges(3,3) + t105 * t93 - t136 * t96 - t124;
t138 = t109 * t6 + t8 * t141 + t4 * t142;
t2 = m(2) * t130 - t122 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t114 * t8 + t120 * t4;
t1 = m(2) * t134 + qJDD(1) * mrSges(2,1) - t122 * mrSges(2,2) - t108 * t6 + (t114 * t4 + t120 * t8) * t109;
t3 = [-m(1) * g(1) - t115 * t1 + t121 * t2, t2, t4, t9, t11, t13, t16; -m(1) * g(2) + t121 * t1 + t115 * t2, t1, t8, t14, t10, t12, t15; (-m(1) - m(2)) * g(3) + t138, -m(2) * g(3) + t138, t6, t124, t123, -t125, -t127;];
f_new  = t3;
