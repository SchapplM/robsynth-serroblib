% Calculate vector of cutting forces with Newton-Euler
% S6RRRRPR5
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-05-07 20:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRPR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 20:25:39
% EndTime: 2019-05-07 20:25:48
% DurationCPUTime: 2.99s
% Computational Cost: add. (39568->207), mult. (79133->259), div. (0->0), fcn. (56112->10), ass. (0->101)
t103 = qJD(2) + qJD(3);
t101 = t103 ^ 2;
t102 = qJDD(2) + qJDD(3);
t107 = sin(qJ(3));
t111 = cos(qJ(3));
t108 = sin(qJ(2));
t112 = cos(qJ(2));
t132 = qJD(1) * qJD(2);
t114 = qJD(1) ^ 2;
t109 = sin(qJ(1));
t113 = cos(qJ(1));
t126 = -t113 * g(1) - t109 * g(2);
t88 = -t114 * pkin(1) + qJDD(1) * pkin(7) + t126;
t135 = t108 * t88;
t144 = pkin(2) * t114;
t92 = t108 * qJDD(1) + t112 * t132;
t55 = qJDD(2) * pkin(2) - t92 * pkin(8) - t135 + (pkin(8) * t132 + t108 * t144 - g(3)) * t112;
t104 = t112 ^ 2;
t131 = -t108 * g(3) + t112 * t88;
t93 = t112 * qJDD(1) - t108 * t132;
t134 = qJD(1) * t108;
t96 = qJD(2) * pkin(2) - pkin(8) * t134;
t56 = t93 * pkin(8) - qJD(2) * t96 - t104 * t144 + t131;
t139 = -t107 * t56 + t111 * t55;
t133 = qJD(1) * t112;
t85 = -t107 * t134 + t111 * t133;
t86 = (t107 * t112 + t108 * t111) * qJD(1);
t74 = -t85 * pkin(3) - t86 * pkin(9);
t119 = t102 * pkin(3) + t101 * pkin(9) - t86 * t74 + t139;
t106 = sin(qJ(4));
t145 = cos(qJ(4));
t76 = -t103 * t145 + t106 * t86;
t84 = qJD(4) - t85;
t143 = t76 * t84;
t65 = t85 * qJD(3) + t107 * t93 + t111 * t92;
t39 = -t76 * qJD(4) + t106 * t102 + t145 * t65;
t150 = (-t39 + t143) * qJ(5) - t119;
t105 = sin(qJ(6));
t110 = cos(qJ(6));
t129 = t109 * g(1) - t113 * g(2);
t124 = -qJDD(1) * pkin(1) - t129;
t116 = -t93 * pkin(2) + t96 * t134 + (-pkin(8) * t104 - pkin(7)) * t114 + t124;
t64 = -t86 * qJD(3) - t107 * t92 + t111 * t93;
t29 = (-t103 * t85 - t65) * pkin(9) + (t103 * t86 - t64) * pkin(3) + t116;
t138 = t107 * t55 + t111 * t56;
t33 = -t101 * pkin(3) + t102 * pkin(9) + t85 * t74 + t138;
t127 = -t106 * t33 + t145 * t29;
t77 = t106 * t103 + t145 * t86;
t52 = t76 * pkin(4) - t77 * qJ(5);
t62 = qJDD(4) - t64;
t83 = t84 ^ 2;
t21 = -t62 * pkin(4) - t83 * qJ(5) + t77 * t52 + qJDD(5) - t127;
t16 = (-t39 - t143) * pkin(10) + (t76 * t77 - t62) * pkin(5) + t21;
t140 = t106 * t29 + t145 * t33;
t147 = 2 * qJD(5);
t122 = -t83 * pkin(4) + t62 * qJ(5) + t147 * t84 - t76 * t52 + t140;
t38 = t77 * qJD(4) - t102 * t145 + t106 * t65;
t71 = -t84 * pkin(5) - t77 * pkin(10);
t75 = t76 ^ 2;
t17 = -t75 * pkin(5) + t38 * pkin(10) + t84 * t71 + t122;
t46 = -t105 * t77 + t110 * t76;
t27 = t46 * qJD(6) + t105 * t38 + t110 * t39;
t47 = t105 * t76 + t110 * t77;
t36 = -t46 * mrSges(7,1) + t47 * mrSges(7,2);
t82 = qJD(6) - t84;
t42 = -t82 * mrSges(7,2) + t46 * mrSges(7,3);
t59 = qJDD(6) - t62;
t14 = m(7) * (-t105 * t17 + t110 * t16) - t27 * mrSges(7,3) + t59 * mrSges(7,1) - t47 * t36 + t82 * t42;
t26 = -t47 * qJD(6) - t105 * t39 + t110 * t38;
t43 = t82 * mrSges(7,1) - t47 * mrSges(7,3);
t15 = m(7) * (t105 * t16 + t110 * t17) + t26 * mrSges(7,3) - t59 * mrSges(7,2) + t46 * t36 - t82 * t43;
t120 = -m(6) * t21 - t105 * t15 - t110 * t14;
t53 = t76 * mrSges(6,1) - t77 * mrSges(6,3);
t137 = -t76 * mrSges(5,1) - t77 * mrSges(5,2) - t53;
t141 = -mrSges(5,3) - mrSges(6,2);
t67 = -t76 * mrSges(6,2) + t84 * mrSges(6,3);
t68 = -t84 * mrSges(5,2) - t76 * mrSges(5,3);
t11 = m(5) * t127 + (t68 + t67) * t84 + t137 * t77 + (mrSges(5,1) + mrSges(6,1)) * t62 + t141 * t39 + t120;
t78 = -t103 * mrSges(4,2) + t85 * mrSges(4,3);
t79 = t103 * mrSges(4,1) - t86 * mrSges(4,3);
t70 = -t84 * mrSges(6,1) + t77 * mrSges(6,2);
t123 = m(6) * t122 + t62 * mrSges(6,3) - t105 * t14 + t110 * t15 + t84 * t70;
t69 = t84 * mrSges(5,1) - t77 * mrSges(5,3);
t9 = m(5) * t140 - t62 * mrSges(5,2) + t137 * t76 + t141 * t38 - t84 * t69 + t123;
t118 = -m(4) * t116 + t64 * mrSges(4,1) - t65 * mrSges(4,2) - t106 * t9 - t145 * t11 + t85 * t78 - t86 * t79;
t94 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t134;
t95 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t133;
t149 = (t108 * t94 - t112 * t95) * qJD(1) + m(3) * (-t114 * pkin(7) + t124) - t93 * mrSges(3,1) + t92 * mrSges(3,2) - t118;
t128 = m(7) * (-t75 * pkin(10) + (-pkin(4) - pkin(5)) * t38 + (-pkin(4) * t84 + t147 + t71) * t77 - t150) + t27 * mrSges(7,2) - t26 * mrSges(7,1) + t47 * t43 - t46 * t42;
t121 = m(6) * (-0.2e1 * qJD(5) * t77 + (t77 * t84 + t38) * pkin(4) + t150) + t38 * mrSges(6,1) + t76 * t67 - t128;
t148 = -m(5) * t119 + t38 * mrSges(5,1) + (t69 - t70) * t77 + (mrSges(5,2) - mrSges(6,3)) * t39 + t76 * t68 + t121;
t73 = -t85 * mrSges(4,1) + t86 * mrSges(4,2);
t12 = m(4) * t139 + t102 * mrSges(4,1) - t65 * mrSges(4,3) + t103 * t78 - t86 * t73 - t148;
t7 = m(4) * t138 - t102 * mrSges(4,2) + t64 * mrSges(4,3) - t103 * t79 - t106 * t11 + t145 * t9 + t85 * t73;
t91 = (-mrSges(3,1) * t112 + mrSges(3,2) * t108) * qJD(1);
t4 = m(3) * (-t112 * g(3) - t135) - t92 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t91 * t134 + qJD(2) * t95 + t107 * t7 + t111 * t12;
t5 = m(3) * t131 - qJDD(2) * mrSges(3,2) + t93 * mrSges(3,3) - qJD(2) * t94 - t107 * t12 + t111 * t7 + t133 * t91;
t146 = t108 * t5 + t112 * t4;
t6 = m(2) * t129 + qJDD(1) * mrSges(2,1) - t114 * mrSges(2,2) - t149;
t1 = m(2) * t126 - t114 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t108 * t4 + t112 * t5;
t2 = [-m(1) * g(1) + t113 * t1 - t109 * t6, t1, t5, t7, t9, -t38 * mrSges(6,2) - t76 * t53 + t123, t15; -m(1) * g(2) + t109 * t1 + t113 * t6, t6, t4, t12, t11, -t39 * mrSges(6,3) - t77 * t70 + t121, t14; (-m(1) - m(2)) * g(3) + t146, -m(2) * g(3) + t146, t149, -t118, t148, -t62 * mrSges(6,1) + t39 * mrSges(6,2) + t77 * t53 - t84 * t67 - t120, t128;];
f_new  = t2;
