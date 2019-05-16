% Calculate vector of cutting forces with Newton-Euler
% S6RPRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-05-06 03:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRRR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 03:12:32
% EndTime: 2019-05-06 03:12:44
% DurationCPUTime: 6.58s
% Computational Cost: add. (106077->196), mult. (264737->253), div. (0->0), fcn. (213598->12), ass. (0->106)
t109 = qJD(1) ^ 2;
t98 = cos(pkin(11));
t95 = t98 ^ 2;
t97 = sin(pkin(11));
t133 = t97 ^ 2 + t95;
t139 = t133 * mrSges(3,3);
t102 = sin(qJ(3));
t107 = cos(qJ(3));
t123 = -t98 * mrSges(3,1) + t97 * mrSges(3,2);
t118 = qJDD(1) * mrSges(3,3) + t109 * t123;
t130 = qJD(1) * qJD(2);
t129 = -t98 * g(3) - 0.2e1 * t97 * t130;
t101 = sin(qJ(4));
t106 = cos(qJ(4));
t132 = pkin(7) * qJDD(1);
t137 = pkin(2) * t109;
t103 = sin(qJ(1));
t108 = cos(qJ(1));
t122 = -t108 * g(1) - t103 * g(2);
t85 = -t109 * pkin(1) + qJDD(1) * qJ(2) + t122;
t64 = (t98 * t137 - t132 - t85) * t97 + t129;
t125 = -t97 * g(3) + (0.2e1 * t130 + t85) * t98;
t65 = t98 * t132 - t95 * t137 + t125;
t126 = -t102 * t65 + t107 * t64;
t119 = -t102 * t97 + t107 * t98;
t83 = t119 * qJD(1);
t120 = t102 * t98 + t107 * t97;
t84 = t120 * qJD(1);
t72 = -t83 * mrSges(4,1) + t84 * mrSges(4,2);
t131 = t83 * qJD(3);
t76 = t120 * qJDD(1) + t131;
t77 = -qJD(3) * mrSges(4,2) + t83 * mrSges(4,3);
t100 = sin(qJ(5));
t105 = cos(qJ(5));
t104 = cos(qJ(6));
t37 = (-t76 + t131) * pkin(8) + (t83 * t84 + qJDD(3)) * pkin(3) + t126;
t134 = t102 * t64 + t107 * t65;
t75 = -t84 * qJD(3) + t119 * qJDD(1);
t79 = qJD(3) * pkin(3) - t84 * pkin(8);
t82 = t83 ^ 2;
t42 = -t82 * pkin(3) + t75 * pkin(8) - qJD(3) * t79 + t134;
t127 = -t101 * t42 + t106 * t37;
t67 = -t101 * t84 + t106 * t83;
t47 = t67 * qJD(4) + t101 * t75 + t106 * t76;
t68 = t101 * t83 + t106 * t84;
t93 = qJDD(3) + qJDD(4);
t96 = qJD(3) + qJD(4);
t21 = (t67 * t96 - t47) * pkin(9) + (t67 * t68 + t93) * pkin(4) + t127;
t135 = t101 * t37 + t106 * t42;
t46 = -t68 * qJD(4) - t101 * t76 + t106 * t75;
t62 = t96 * pkin(4) - t68 * pkin(9);
t66 = t67 ^ 2;
t23 = -t66 * pkin(4) + t46 * pkin(9) - t96 * t62 + t135;
t136 = t100 * t21 + t105 * t23;
t55 = -t100 * t68 + t105 * t67;
t56 = t100 * t67 + t105 * t68;
t41 = -t55 * pkin(5) - t56 * pkin(10);
t91 = qJD(5) + t96;
t89 = t91 ^ 2;
t90 = qJDD(5) + t93;
t18 = -t89 * pkin(5) + t90 * pkin(10) + t55 * t41 + t136;
t128 = t103 * g(1) - t108 * g(2);
t124 = qJDD(2) - t128;
t114 = (-pkin(2) * t98 - pkin(1)) * qJDD(1) + (-t133 * pkin(7) - qJ(2)) * t109 + t124;
t113 = -t75 * pkin(3) - t82 * pkin(8) + t84 * t79 + t114;
t110 = -t46 * pkin(4) - t66 * pkin(9) + t68 * t62 + t113;
t30 = -t56 * qJD(5) - t100 * t47 + t105 * t46;
t31 = t55 * qJD(5) + t100 * t46 + t105 * t47;
t19 = t110 + (-t55 * t91 - t31) * pkin(10) + (t56 * t91 - t30) * pkin(5);
t99 = sin(qJ(6));
t48 = t104 * t91 - t99 * t56;
t25 = t48 * qJD(6) + t104 * t31 + t99 * t90;
t29 = qJDD(6) - t30;
t49 = t104 * t56 + t99 * t91;
t32 = -t48 * mrSges(7,1) + t49 * mrSges(7,2);
t52 = qJD(6) - t55;
t35 = -t52 * mrSges(7,2) + t48 * mrSges(7,3);
t15 = m(7) * (t104 * t19 - t99 * t18) - t25 * mrSges(7,3) + t29 * mrSges(7,1) - t49 * t32 + t52 * t35;
t24 = -t49 * qJD(6) + t104 * t90 - t99 * t31;
t36 = t52 * mrSges(7,1) - t49 * mrSges(7,3);
t16 = m(7) * (t104 * t18 + t99 * t19) + t24 * mrSges(7,3) - t29 * mrSges(7,2) + t48 * t32 - t52 * t36;
t39 = -t55 * mrSges(6,1) + t56 * mrSges(6,2);
t51 = t91 * mrSges(6,1) - t56 * mrSges(6,3);
t11 = m(6) * t136 - t90 * mrSges(6,2) + t30 * mrSges(6,3) + t104 * t16 - t99 * t15 + t55 * t39 - t91 * t51;
t121 = -t100 * t23 + t105 * t21;
t116 = m(7) * (-t90 * pkin(5) - t89 * pkin(10) + t56 * t41 - t121) - t24 * mrSges(7,1) + t25 * mrSges(7,2) - t48 * t35 + t49 * t36;
t50 = -t91 * mrSges(6,2) + t55 * mrSges(6,3);
t12 = m(6) * t121 + t90 * mrSges(6,1) - t31 * mrSges(6,3) - t56 * t39 + t91 * t50 - t116;
t57 = -t67 * mrSges(5,1) + t68 * mrSges(5,2);
t60 = -t96 * mrSges(5,2) + t67 * mrSges(5,3);
t8 = m(5) * t127 + t93 * mrSges(5,1) - t47 * mrSges(5,3) + t100 * t11 + t105 * t12 - t68 * t57 + t96 * t60;
t61 = t96 * mrSges(5,1) - t68 * mrSges(5,3);
t9 = m(5) * t135 - t93 * mrSges(5,2) + t46 * mrSges(5,3) - t100 * t12 + t105 * t11 + t67 * t57 - t96 * t61;
t6 = m(4) * t126 + qJDD(3) * mrSges(4,1) - t76 * mrSges(4,3) + qJD(3) * t77 + t101 * t9 + t106 * t8 - t84 * t72;
t78 = qJD(3) * mrSges(4,1) - t84 * mrSges(4,3);
t7 = m(4) * t134 - qJDD(3) * mrSges(4,2) + t75 * mrSges(4,3) - qJD(3) * t78 - t101 * t8 + t106 * t9 + t83 * t72;
t4 = m(3) * t129 + t102 * t7 + t107 * t6 + (-m(3) * t85 - t118) * t97;
t5 = m(3) * t125 - t102 * t6 + t107 * t7 + t118 * t98;
t138 = t98 * t4 + t97 * t5;
t117 = -m(6) * t110 + t30 * mrSges(6,1) - t31 * mrSges(6,2) - t104 * t15 - t99 * t16 + t55 * t50 - t56 * t51;
t115 = -m(5) * t113 + t46 * mrSges(5,1) - t47 * mrSges(5,2) + t67 * t60 - t68 * t61 + t117;
t112 = -m(4) * t114 + t75 * mrSges(4,1) - t76 * mrSges(4,2) + t83 * t77 - t84 * t78 + t115;
t111 = m(3) * (-qJDD(1) * pkin(1) - t109 * qJ(2) + t124) - t112;
t10 = -t111 + (-mrSges(2,2) + t139) * t109 + (mrSges(2,1) - t123) * qJDD(1) + m(2) * t128;
t1 = m(2) * t122 - t109 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t97 * t4 + t98 * t5;
t2 = [-m(1) * g(1) + t108 * t1 - t103 * t10, t1, t5, t7, t9, t11, t16; -m(1) * g(2) + t103 * t1 + t108 * t10, t10, t4, t6, t8, t12, t15; (-m(1) - m(2)) * g(3) + t138, -m(2) * g(3) + t138, t123 * qJDD(1) - t109 * t139 + t111, -t112, -t115, -t117, t116;];
f_new  = t2;
