% Calculate vector of cutting forces with Newton-Euler
% S6RRRPRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-05-07 15:13
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPRR12_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR12_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR12_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR12_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR12_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR12_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR12_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 14:55:54
% EndTime: 2019-05-07 14:56:24
% DurationCPUTime: 12.07s
% Computational Cost: add. (221886->215), mult. (484903->296), div. (0->0), fcn. (390970->14), ass. (0->116)
t108 = sin(pkin(6));
t114 = sin(qJ(2));
t119 = cos(qJ(2));
t137 = qJD(1) * qJD(2);
t97 = (-qJDD(1) * t119 + t114 * t137) * t108;
t148 = pkin(8) * t108;
t110 = cos(pkin(6));
t147 = t110 * g(3);
t112 = sin(qJ(5));
t117 = cos(qJ(5));
t107 = sin(pkin(12));
t109 = cos(pkin(12));
t139 = qJD(1) * t119;
t134 = t108 * t139;
t101 = qJD(3) - t134;
t100 = t101 ^ 2;
t113 = sin(qJ(3));
t118 = cos(qJ(3));
t104 = t110 * qJD(1) + qJD(2);
t102 = t104 ^ 2;
t103 = t110 * qJDD(1) + qJDD(2);
t121 = qJD(1) ^ 2;
t115 = sin(qJ(1));
t120 = cos(qJ(1));
t133 = t115 * g(1) - t120 * g(2);
t92 = qJDD(1) * pkin(1) + t121 * t148 + t133;
t143 = t110 * t92;
t129 = -t120 * g(1) - t115 * g(2);
t93 = -t121 * pkin(1) + qJDD(1) * t148 + t129;
t144 = t114 * t143 + t119 * t93;
t140 = qJD(1) * t108;
t95 = (-pkin(2) * t119 - pkin(9) * t114) * t140;
t56 = -t102 * pkin(2) + t103 * pkin(9) + (-g(3) * t114 + t95 * t139) * t108 + t144;
t96 = (qJDD(1) * t114 + t119 * t137) * t108;
t57 = t97 * pkin(2) - t96 * pkin(9) - t147 + (-t92 + (pkin(2) * t114 - pkin(9) * t119) * t104 * qJD(1)) * t108;
t145 = t113 * t57 + t118 * t56;
t135 = t114 * t140;
t85 = -t118 * t104 + t113 * t135;
t86 = t113 * t104 + t118 * t135;
t72 = t85 * pkin(3) - t86 * qJ(4);
t89 = qJDD(3) + t97;
t35 = -t100 * pkin(3) + t89 * qJ(4) - t85 * t72 + t145;
t141 = t108 * t119;
t128 = -g(3) * t141 - t114 * t93 + t119 * t143;
t55 = -t103 * pkin(2) - t102 * pkin(9) + t95 * t135 - t128;
t70 = t86 * qJD(3) - t118 * t103 + t113 * t96;
t71 = -t85 * qJD(3) + t113 * t103 + t118 * t96;
t38 = (t101 * t85 - t71) * qJ(4) + (t101 * t86 + t70) * pkin(3) + t55;
t78 = t107 * t101 + t109 * t86;
t130 = -0.2e1 * qJD(4) * t78 - t107 * t35 + t109 * t38;
t62 = t107 * t89 + t109 * t71;
t77 = t109 * t101 - t107 * t86;
t23 = (t77 * t85 - t62) * pkin(10) + (t77 * t78 + t70) * pkin(4) + t130;
t136 = 0.2e1 * qJD(4) * t77 + t107 * t38 + t109 * t35;
t61 = -t107 * t71 + t109 * t89;
t67 = t85 * pkin(4) - t78 * pkin(10);
t76 = t77 ^ 2;
t25 = -t76 * pkin(4) + t61 * pkin(10) - t85 * t67 + t136;
t146 = t112 * t23 + t117 * t25;
t142 = t108 * t114;
t131 = -t113 * t56 + t118 * t57;
t34 = -t89 * pkin(3) - t100 * qJ(4) + t86 * t72 + qJDD(4) - t131;
t124 = -t61 * pkin(4) - t76 * pkin(10) + t78 * t67 + t34;
t111 = sin(qJ(6));
t116 = cos(qJ(6));
t60 = t112 * t77 + t117 * t78;
t40 = -t60 * qJD(5) - t112 * t62 + t117 * t61;
t59 = -t112 * t78 + t117 * t77;
t41 = t59 * qJD(5) + t112 * t61 + t117 * t62;
t46 = t111 * t59 + t116 * t60;
t27 = -t46 * qJD(6) - t111 * t41 + t116 * t40;
t45 = -t111 * t60 + t116 * t59;
t28 = t45 * qJD(6) + t111 * t40 + t116 * t41;
t84 = qJD(5) + t85;
t83 = qJD(6) + t84;
t42 = -t83 * mrSges(7,2) + t45 * mrSges(7,3);
t43 = t83 * mrSges(7,1) - t46 * mrSges(7,3);
t50 = t84 * pkin(5) - t60 * pkin(11);
t58 = t59 ^ 2;
t126 = t27 * mrSges(7,1) + t45 * t42 - m(7) * (-t40 * pkin(5) - t58 * pkin(11) + t60 * t50 + t124) - t28 * mrSges(7,2) - t46 * t43;
t48 = -t84 * mrSges(6,2) + t59 * mrSges(6,3);
t49 = t84 * mrSges(6,1) - t60 * mrSges(6,3);
t125 = -m(6) * t124 + t40 * mrSges(6,1) - t41 * mrSges(6,2) + t59 * t48 - t60 * t49 + t126;
t65 = -t85 * mrSges(5,2) + t77 * mrSges(5,3);
t66 = t85 * mrSges(5,1) - t78 * mrSges(5,3);
t122 = m(5) * t34 - t61 * mrSges(5,1) + t62 * mrSges(5,2) - t77 * t65 + t78 * t66 - t125;
t73 = t85 * mrSges(4,1) + t86 * mrSges(4,2);
t79 = -t101 * mrSges(4,2) - t85 * mrSges(4,3);
t16 = m(4) * t131 + t89 * mrSges(4,1) - t71 * mrSges(4,3) + t101 * t79 - t86 * t73 - t122;
t132 = -t112 * t25 + t117 * t23;
t69 = qJDD(5) + t70;
t17 = (t59 * t84 - t41) * pkin(11) + (t59 * t60 + t69) * pkin(5) + t132;
t18 = -t58 * pkin(5) + t40 * pkin(11) - t84 * t50 + t146;
t32 = -t45 * mrSges(7,1) + t46 * mrSges(7,2);
t68 = qJDD(6) + t69;
t14 = m(7) * (-t111 * t18 + t116 * t17) - t28 * mrSges(7,3) + t68 * mrSges(7,1) - t46 * t32 + t83 * t42;
t15 = m(7) * (t111 * t17 + t116 * t18) + t27 * mrSges(7,3) - t68 * mrSges(7,2) + t45 * t32 - t83 * t43;
t47 = -t59 * mrSges(6,1) + t60 * mrSges(6,2);
t12 = m(6) * t132 + t69 * mrSges(6,1) - t41 * mrSges(6,3) + t111 * t15 + t116 * t14 - t60 * t47 + t84 * t48;
t13 = m(6) * t146 - t69 * mrSges(6,2) + t40 * mrSges(6,3) - t111 * t14 + t116 * t15 + t59 * t47 - t84 * t49;
t63 = -t77 * mrSges(5,1) + t78 * mrSges(5,2);
t10 = m(5) * t130 + t70 * mrSges(5,1) - t62 * mrSges(5,3) + t112 * t13 + t117 * t12 - t78 * t63 + t85 * t65;
t11 = m(5) * t136 - t70 * mrSges(5,2) + t61 * mrSges(5,3) - t112 * t12 + t117 * t13 + t77 * t63 - t85 * t66;
t80 = t101 * mrSges(4,1) - t86 * mrSges(4,3);
t9 = m(4) * t145 - t89 * mrSges(4,2) - t70 * mrSges(4,3) - t107 * t10 - t101 * t80 + t109 * t11 - t85 * t73;
t90 = t104 * mrSges(3,1) - mrSges(3,3) * t135;
t94 = (-mrSges(3,1) * t119 + mrSges(3,2) * t114) * t140;
t4 = m(3) * (-g(3) * t142 + t144) - t97 * mrSges(3,3) - t103 * mrSges(3,2) + t94 * t134 - t104 * t90 + t118 * t9 - t113 * t16;
t91 = -t104 * mrSges(3,2) + mrSges(3,3) * t134;
t6 = m(3) * (-t108 * t92 - t147) + t96 * mrSges(3,2) + t97 * mrSges(3,1) + t113 * t9 + t118 * t16 + (t114 * t90 - t119 * t91) * t140;
t123 = m(4) * t55 + t70 * mrSges(4,1) + t71 * mrSges(4,2) + t109 * t10 + t107 * t11 + t85 * t79 + t86 * t80;
t8 = m(3) * t128 + t103 * mrSges(3,1) - t96 * mrSges(3,3) + t104 * t91 - t94 * t135 - t123;
t138 = t110 * t6 + t8 * t141 + t4 * t142;
t2 = m(2) * t129 - t121 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t114 * t8 + t119 * t4;
t1 = m(2) * t133 + qJDD(1) * mrSges(2,1) - t121 * mrSges(2,2) - t108 * t6 + (t114 * t4 + t119 * t8) * t110;
t3 = [-m(1) * g(1) - t115 * t1 + t120 * t2, t2, t4, t9, t11, t13, t15; -m(1) * g(2) + t120 * t1 + t115 * t2, t1, t8, t16, t10, t12, t14; (-m(1) - m(2)) * g(3) + t138, -m(2) * g(3) + t138, t6, t123, t122, -t125, -t126;];
f_new  = t3;
