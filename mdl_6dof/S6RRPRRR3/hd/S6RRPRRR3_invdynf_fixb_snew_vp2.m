% Calculate vector of cutting forces with Newton-Euler
% S6RRPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-06 20:23
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRRR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 20:16:17
% EndTime: 2019-05-06 20:16:34
% DurationCPUTime: 6.00s
% Computational Cost: add. (105232->208), mult. (241546->275), div. (0->0), fcn. (177682->12), ass. (0->106)
t104 = sin(pkin(11));
t105 = cos(pkin(11));
t109 = sin(qJ(2));
t114 = cos(qJ(2));
t134 = qJD(1) * qJD(2);
t117 = qJD(1) ^ 2;
t110 = sin(qJ(1));
t115 = cos(qJ(1));
t127 = -t115 * g(1) - t110 * g(2);
t93 = -t117 * pkin(1) + qJDD(1) * pkin(7) + t127;
t138 = t109 * t93;
t141 = pkin(2) * t117;
t96 = t109 * qJDD(1) + t114 * t134;
t60 = qJDD(2) * pkin(2) - t96 * qJ(3) - t138 + (qJ(3) * t134 + t109 * t141 - g(3)) * t114;
t103 = t114 ^ 2;
t132 = -t109 * g(3) + t114 * t93;
t97 = t114 * qJDD(1) - t109 * t134;
t136 = qJD(1) * t109;
t98 = qJD(2) * pkin(2) - qJ(3) * t136;
t61 = t97 * qJ(3) - qJD(2) * t98 - t103 * t141 + t132;
t90 = (t104 * t114 + t105 * t109) * qJD(1);
t144 = -0.2e1 * qJD(3) * t90 - t104 * t61 + t105 * t60;
t135 = qJD(1) * t114;
t100 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t135;
t107 = sin(qJ(5));
t112 = cos(qJ(5));
t106 = sin(qJ(6));
t111 = cos(qJ(6));
t108 = sin(qJ(4));
t113 = cos(qJ(4));
t116 = qJD(2) ^ 2;
t89 = -t104 * t136 + t105 * t135;
t133 = 0.2e1 * qJD(3) * t89 + t104 * t60 + t105 * t61;
t72 = -t89 * pkin(3) - t90 * pkin(8);
t38 = -t116 * pkin(3) + qJDD(2) * pkin(8) + t89 * t72 + t133;
t131 = t110 * g(1) - t115 * g(2);
t125 = -qJDD(1) * pkin(1) - t131;
t119 = -t97 * pkin(2) + qJDD(3) + t98 * t136 + (-qJ(3) * t103 - pkin(7)) * t117 + t125;
t76 = -t104 * t96 + t105 * t97;
t77 = t104 * t97 + t105 * t96;
t41 = (-qJD(2) * t89 - t77) * pkin(8) + (qJD(2) * t90 - t76) * pkin(3) + t119;
t128 = -t108 * t38 + t113 * t41;
t79 = t113 * qJD(2) - t108 * t90;
t53 = t79 * qJD(4) + t108 * qJDD(2) + t113 * t77;
t75 = qJDD(4) - t76;
t80 = t108 * qJD(2) + t113 * t90;
t88 = qJD(4) - t89;
t23 = (t79 * t88 - t53) * pkin(9) + (t79 * t80 + t75) * pkin(4) + t128;
t139 = t108 * t41 + t113 * t38;
t52 = -t80 * qJD(4) + t113 * qJDD(2) - t108 * t77;
t68 = t88 * pkin(4) - t80 * pkin(9);
t78 = t79 ^ 2;
t28 = -t78 * pkin(4) + t52 * pkin(9) - t88 * t68 + t139;
t129 = -t107 * t28 + t112 * t23;
t58 = -t107 * t80 + t112 * t79;
t35 = t58 * qJD(5) + t107 * t52 + t112 * t53;
t59 = t107 * t79 + t112 * t80;
t73 = qJDD(5) + t75;
t84 = qJD(5) + t88;
t17 = (t58 * t84 - t35) * pkin(10) + (t58 * t59 + t73) * pkin(5) + t129;
t140 = t107 * t23 + t112 * t28;
t34 = -t59 * qJD(5) - t107 * t53 + t112 * t52;
t50 = t84 * pkin(5) - t59 * pkin(10);
t57 = t58 ^ 2;
t18 = -t57 * pkin(5) + t34 * pkin(10) - t84 * t50 + t140;
t45 = -t106 * t59 + t111 * t58;
t26 = t45 * qJD(6) + t106 * t34 + t111 * t35;
t46 = t106 * t58 + t111 * t59;
t32 = -t45 * mrSges(7,1) + t46 * mrSges(7,2);
t83 = qJD(6) + t84;
t42 = -t83 * mrSges(7,2) + t45 * mrSges(7,3);
t71 = qJDD(6) + t73;
t15 = m(7) * (-t106 * t18 + t111 * t17) - t26 * mrSges(7,3) + t71 * mrSges(7,1) - t46 * t32 + t83 * t42;
t25 = -t46 * qJD(6) - t106 * t35 + t111 * t34;
t43 = t83 * mrSges(7,1) - t46 * mrSges(7,3);
t16 = m(7) * (t106 * t17 + t111 * t18) + t25 * mrSges(7,3) - t71 * mrSges(7,2) + t45 * t32 - t83 * t43;
t47 = -t58 * mrSges(6,1) + t59 * mrSges(6,2);
t48 = -t84 * mrSges(6,2) + t58 * mrSges(6,3);
t12 = m(6) * t129 + t73 * mrSges(6,1) - t35 * mrSges(6,3) + t106 * t16 + t111 * t15 - t59 * t47 + t84 * t48;
t49 = t84 * mrSges(6,1) - t59 * mrSges(6,3);
t13 = m(6) * t140 - t73 * mrSges(6,2) + t34 * mrSges(6,3) - t106 * t15 + t111 * t16 + t58 * t47 - t84 * t49;
t62 = -t79 * mrSges(5,1) + t80 * mrSges(5,2);
t66 = -t88 * mrSges(5,2) + t79 * mrSges(5,3);
t10 = m(5) * t128 + t75 * mrSges(5,1) - t53 * mrSges(5,3) + t107 * t13 + t112 * t12 - t80 * t62 + t88 * t66;
t67 = t88 * mrSges(5,1) - t80 * mrSges(5,3);
t11 = m(5) * t139 - t75 * mrSges(5,2) + t52 * mrSges(5,3) - t107 * t12 + t112 * t13 + t79 * t62 - t88 * t67;
t81 = -qJD(2) * mrSges(4,2) + t89 * mrSges(4,3);
t82 = qJD(2) * mrSges(4,1) - t90 * mrSges(4,3);
t123 = -m(4) * t119 + t76 * mrSges(4,1) - t77 * mrSges(4,2) - t113 * t10 - t108 * t11 + t89 * t81 - t90 * t82;
t99 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t136;
t143 = -(t114 * t100 - t109 * t99) * qJD(1) + m(3) * (-t117 * pkin(7) + t125) - t97 * mrSges(3,1) + t96 * mrSges(3,2) - t123;
t37 = -qJDD(2) * pkin(3) - t116 * pkin(8) + t90 * t72 - t144;
t121 = -t52 * pkin(4) - t78 * pkin(9) + t80 * t68 + t37;
t124 = t25 * mrSges(7,1) + t45 * t42 - m(7) * (-t34 * pkin(5) - t57 * pkin(10) + t59 * t50 + t121) - t26 * mrSges(7,2) - t46 * t43;
t122 = -m(6) * t121 + t34 * mrSges(6,1) - t35 * mrSges(6,2) + t58 * t48 - t59 * t49 + t124;
t118 = m(5) * t37 - t52 * mrSges(5,1) + t53 * mrSges(5,2) - t79 * t66 + t80 * t67 - t122;
t70 = -t89 * mrSges(4,1) + t90 * mrSges(4,2);
t14 = m(4) * t144 + qJDD(2) * mrSges(4,1) - t77 * mrSges(4,3) + qJD(2) * t81 - t90 * t70 - t118;
t7 = m(4) * t133 - qJDD(2) * mrSges(4,2) + t76 * mrSges(4,3) - qJD(2) * t82 - t108 * t10 + t113 * t11 + t89 * t70;
t95 = (-mrSges(3,1) * t114 + mrSges(3,2) * t109) * qJD(1);
t4 = m(3) * (-t114 * g(3) - t138) - t96 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t95 * t136 + qJD(2) * t100 + t104 * t7 + t105 * t14;
t5 = m(3) * t132 - qJDD(2) * mrSges(3,2) + t97 * mrSges(3,3) - qJD(2) * t99 - t104 * t14 + t105 * t7 + t95 * t135;
t142 = t109 * t5 + t114 * t4;
t6 = m(2) * t131 + qJDD(1) * mrSges(2,1) - t117 * mrSges(2,2) - t143;
t1 = m(2) * t127 - t117 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t109 * t4 + t114 * t5;
t2 = [-m(1) * g(1) + t115 * t1 - t110 * t6, t1, t5, t7, t11, t13, t16; -m(1) * g(2) + t110 * t1 + t115 * t6, t6, t4, t14, t10, t12, t15; (-m(1) - m(2)) * g(3) + t142, -m(2) * g(3) + t142, t143, -t123, t118, -t122, -t124;];
f_new  = t2;
