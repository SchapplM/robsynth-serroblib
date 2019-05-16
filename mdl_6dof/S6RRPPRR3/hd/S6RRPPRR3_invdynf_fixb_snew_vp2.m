% Calculate vector of cutting forces with Newton-Euler
% S6RRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-05-06 10:07
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPPRR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:59:20
% EndTime: 2019-05-06 09:59:40
% DurationCPUTime: 8.12s
% Computational Cost: add. (129175->212), mult. (349078->296), div. (0->0), fcn. (278144->14), ass. (0->115)
t154 = -2 * qJD(3);
t110 = sin(pkin(11));
t113 = cos(pkin(11));
t111 = sin(pkin(6));
t117 = sin(qJ(2));
t121 = cos(qJ(2));
t143 = qJD(1) * qJD(2);
t100 = (qJDD(1) * t117 + t121 * t143) * t111;
t114 = cos(pkin(6));
t104 = t114 * qJDD(1) + qJDD(2);
t105 = t114 * qJD(1) + qJD(2);
t123 = qJD(1) ^ 2;
t118 = sin(qJ(1));
t122 = cos(qJ(1));
t137 = t118 * g(1) - t122 * g(2);
t152 = pkin(8) * t111;
t97 = qJDD(1) * pkin(1) + t123 * t152 + t137;
t150 = t114 * t97;
t133 = -t122 * g(1) - t118 * g(2);
t98 = -t123 * pkin(1) + qJDD(1) * t152 + t133;
t135 = -t117 * t98 + t121 * t150;
t148 = t111 ^ 2 * t123;
t53 = t104 * pkin(2) - t100 * qJ(3) + (pkin(2) * t117 * t148 + (qJ(3) * qJD(1) * t105 - g(3)) * t111) * t121 + t135;
t101 = (qJDD(1) * t121 - t117 * t143) * t111;
t147 = t111 * t117;
t130 = -g(3) * t147 + t117 * t150 + t121 * t98;
t140 = t121 ^ 2 * t148;
t145 = qJD(1) * t111;
t139 = t117 * t145;
t94 = t105 * pkin(2) - qJ(3) * t139;
t56 = -pkin(2) * t140 + t101 * qJ(3) - t105 * t94 + t130;
t91 = (t110 * t121 + t113 * t117) * t145;
t153 = -t110 * t56 + t113 * t53 + t91 * t154;
t116 = sin(qJ(5));
t120 = cos(qJ(5));
t109 = sin(pkin(12));
t112 = cos(pkin(12));
t103 = t105 ^ 2;
t138 = t121 * t145;
t90 = t110 * t139 - t113 * t138;
t141 = t110 * t53 + t113 * t56 + t90 * t154;
t70 = t90 * pkin(3) - t91 * qJ(4);
t34 = -t103 * pkin(3) + t104 * qJ(4) - t90 * t70 + t141;
t132 = -t114 * g(3) - t111 * t97;
t126 = -t101 * pkin(2) - qJ(3) * t140 + t94 * t139 + qJDD(3) + t132;
t74 = t110 * t100 - t113 * t101;
t75 = t113 * t100 + t110 * t101;
t37 = (t105 * t90 - t75) * qJ(4) + (t105 * t91 + t74) * pkin(3) + t126;
t80 = t109 * t105 + t112 * t91;
t134 = -0.2e1 * qJD(4) * t80 - t109 * t34 + t112 * t37;
t68 = t109 * t104 + t112 * t75;
t79 = t112 * t105 - t109 * t91;
t25 = (t79 * t90 - t68) * pkin(9) + (t79 * t80 + t74) * pkin(4) + t134;
t142 = 0.2e1 * qJD(4) * t79 + t109 * t37 + t112 * t34;
t66 = t90 * pkin(4) - t80 * pkin(9);
t67 = t112 * t104 - t109 * t75;
t78 = t79 ^ 2;
t27 = -t78 * pkin(4) + t67 * pkin(9) - t90 * t66 + t142;
t151 = t116 * t25 + t120 * t27;
t146 = t111 * t121;
t115 = sin(qJ(6));
t119 = cos(qJ(6));
t33 = -t104 * pkin(3) - t103 * qJ(4) + t91 * t70 + qJDD(4) - t153;
t125 = -t67 * pkin(4) - t78 * pkin(9) + t80 * t66 + t33;
t60 = -t116 * t80 + t120 * t79;
t61 = t116 * t79 + t120 * t80;
t47 = -t60 * pkin(5) - t61 * pkin(10);
t73 = qJDD(5) + t74;
t89 = qJD(5) + t90;
t88 = t89 ^ 2;
t22 = -t88 * pkin(5) + t73 * pkin(10) + t60 * t47 + t151;
t41 = -t61 * qJD(5) - t116 * t68 + t120 * t67;
t42 = t60 * qJD(5) + t116 * t67 + t120 * t68;
t23 = (-t60 * t89 - t42) * pkin(10) + (t61 * t89 - t41) * pkin(5) + t125;
t49 = -t115 * t61 + t119 * t89;
t31 = t49 * qJD(6) + t115 * t73 + t119 * t42;
t50 = t115 * t89 + t119 * t61;
t38 = -t49 * mrSges(7,1) + t50 * mrSges(7,2);
t40 = qJDD(6) - t41;
t59 = qJD(6) - t60;
t43 = -t59 * mrSges(7,2) + t49 * mrSges(7,3);
t19 = m(7) * (-t115 * t22 + t119 * t23) - t31 * mrSges(7,3) + t40 * mrSges(7,1) - t50 * t38 + t59 * t43;
t30 = -t50 * qJD(6) - t115 * t42 + t119 * t73;
t44 = t59 * mrSges(7,1) - t50 * mrSges(7,3);
t20 = m(7) * (t115 * t23 + t119 * t22) + t30 * mrSges(7,3) - t40 * mrSges(7,2) + t49 * t38 - t59 * t44;
t54 = -t89 * mrSges(6,2) + t60 * mrSges(6,3);
t55 = t89 * mrSges(6,1) - t61 * mrSges(6,3);
t128 = -m(6) * t125 + t41 * mrSges(6,1) - t42 * mrSges(6,2) - t115 * t20 - t119 * t19 + t60 * t54 - t61 * t55;
t64 = -t90 * mrSges(5,2) + t79 * mrSges(5,3);
t65 = t90 * mrSges(5,1) - t80 * mrSges(5,3);
t124 = m(5) * t33 - t67 * mrSges(5,1) + t68 * mrSges(5,2) - t79 * t64 + t80 * t65 - t128;
t71 = t90 * mrSges(4,1) + t91 * mrSges(4,2);
t81 = -t105 * mrSges(4,2) - t90 * mrSges(4,3);
t14 = m(4) * t153 + t104 * mrSges(4,1) - t75 * mrSges(4,3) + t105 * t81 - t91 * t71 - t124;
t46 = -t60 * mrSges(6,1) + t61 * mrSges(6,2);
t15 = m(6) * t151 - t73 * mrSges(6,2) + t41 * mrSges(6,3) - t115 * t19 + t119 * t20 + t60 * t46 - t89 * t55;
t131 = -t116 * t27 + t120 * t25;
t127 = m(7) * (-t73 * pkin(5) - t88 * pkin(10) + t61 * t47 - t131) - t30 * mrSges(7,1) + t31 * mrSges(7,2) - t49 * t43 + t50 * t44;
t16 = m(6) * t131 + t73 * mrSges(6,1) - t42 * mrSges(6,3) - t61 * t46 + t89 * t54 - t127;
t62 = -t79 * mrSges(5,1) + t80 * mrSges(5,2);
t12 = m(5) * t134 + t74 * mrSges(5,1) - t68 * mrSges(5,3) + t116 * t15 + t120 * t16 - t80 * t62 + t90 * t64;
t13 = m(5) * t142 - t74 * mrSges(5,2) + t67 * mrSges(5,3) - t116 * t16 + t120 * t15 + t79 * t62 - t90 * t65;
t82 = t105 * mrSges(4,1) - t91 * mrSges(4,3);
t7 = m(4) * t141 - t104 * mrSges(4,2) - t74 * mrSges(4,3) - t105 * t82 - t109 * t12 + t112 * t13 - t90 * t71;
t96 = -t105 * mrSges(3,2) + mrSges(3,3) * t138;
t99 = (-mrSges(3,1) * t121 + mrSges(3,2) * t117) * t145;
t5 = m(3) * (-g(3) * t146 + t135) - t100 * mrSges(3,3) + t104 * mrSges(3,1) - t99 * t139 + t105 * t96 + t110 * t7 + t113 * t14;
t95 = t105 * mrSges(3,1) - mrSges(3,3) * t139;
t6 = m(3) * t130 - t104 * mrSges(3,2) + t101 * mrSges(3,3) - t105 * t95 - t110 * t14 + t113 * t7 + t99 * t138;
t129 = m(4) * t126 + t74 * mrSges(4,1) + t75 * mrSges(4,2) + t109 * t13 + t112 * t12 + t90 * t81 + t91 * t82;
t9 = m(3) * t132 + t100 * mrSges(3,2) - t101 * mrSges(3,1) + (t117 * t95 - t121 * t96) * t145 + t129;
t144 = t114 * t9 + t5 * t146 + t6 * t147;
t2 = m(2) * t133 - t123 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t117 * t5 + t121 * t6;
t1 = m(2) * t137 + qJDD(1) * mrSges(2,1) - t123 * mrSges(2,2) - t111 * t9 + (t117 * t6 + t121 * t5) * t114;
t3 = [-m(1) * g(1) - t118 * t1 + t122 * t2, t2, t6, t7, t13, t15, t20; -m(1) * g(2) + t122 * t1 + t118 * t2, t1, t5, t14, t12, t16, t19; (-m(1) - m(2)) * g(3) + t144, -m(2) * g(3) + t144, t9, t129, t124, -t128, t127;];
f_new  = t3;
