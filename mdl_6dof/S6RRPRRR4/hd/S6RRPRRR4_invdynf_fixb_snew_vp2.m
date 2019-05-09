% Calculate vector of cutting forces with Newton-Euler
% S6RRPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-06 20:52
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRRR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 20:39:18
% EndTime: 2019-05-06 20:39:36
% DurationCPUTime: 8.44s
% Computational Cost: add. (141009->213), mult. (368134->295), div. (0->0), fcn. (295919->14), ass. (0->116)
t109 = sin(pkin(12));
t111 = cos(pkin(12));
t110 = sin(pkin(6));
t116 = sin(qJ(2));
t121 = cos(qJ(2));
t142 = qJD(1) * qJD(2);
t100 = (qJDD(1) * t116 + t121 * t142) * t110;
t112 = cos(pkin(6));
t104 = qJDD(1) * t112 + qJDD(2);
t105 = qJD(1) * t112 + qJD(2);
t123 = qJD(1) ^ 2;
t117 = sin(qJ(1));
t122 = cos(qJ(1));
t137 = t117 * g(1) - g(2) * t122;
t152 = pkin(8) * t110;
t97 = qJDD(1) * pkin(1) + t123 * t152 + t137;
t149 = t112 * t97;
t133 = -g(1) * t122 - g(2) * t117;
t98 = -pkin(1) * t123 + qJDD(1) * t152 + t133;
t134 = -t116 * t98 + t121 * t149;
t147 = t110 ^ 2 * t123;
t53 = t104 * pkin(2) - t100 * qJ(3) + (pkin(2) * t116 * t147 + (qJ(3) * qJD(1) * t105 - g(3)) * t110) * t121 + t134;
t101 = (qJDD(1) * t121 - t116 * t142) * t110;
t146 = t110 * t116;
t130 = -g(3) * t146 + t116 * t149 + t121 * t98;
t140 = t121 ^ 2 * t147;
t144 = qJD(1) * t110;
t139 = t116 * t144;
t94 = pkin(2) * t105 - qJ(3) * t139;
t56 = -pkin(2) * t140 + qJ(3) * t101 - t105 * t94 + t130;
t91 = (t109 * t121 + t111 * t116) * t144;
t153 = -0.2e1 * qJD(3) * t91 - t109 * t56 + t111 * t53;
t114 = sin(qJ(5));
t119 = cos(qJ(5));
t115 = sin(qJ(4));
t120 = cos(qJ(4));
t103 = t105 ^ 2;
t138 = t121 * t144;
t90 = -t109 * t139 + t111 * t138;
t141 = 0.2e1 * qJD(3) * t90 + t109 * t53 + t111 * t56;
t71 = -pkin(3) * t90 - pkin(9) * t91;
t34 = -pkin(3) * t103 + pkin(9) * t104 + t71 * t90 + t141;
t132 = -t112 * g(3) - t110 * t97;
t126 = -t101 * pkin(2) - qJ(3) * t140 + t94 * t139 + qJDD(3) + t132;
t75 = -t109 * t100 + t101 * t111;
t76 = t100 * t111 + t101 * t109;
t41 = (-t105 * t90 - t76) * pkin(9) + (t105 * t91 - t75) * pkin(3) + t126;
t135 = -t115 * t34 + t120 * t41;
t78 = t105 * t120 - t115 * t91;
t59 = qJD(4) * t78 + t104 * t115 + t120 * t76;
t74 = qJDD(4) - t75;
t79 = t105 * t115 + t120 * t91;
t89 = qJD(4) - t90;
t25 = (t78 * t89 - t59) * pkin(10) + (t78 * t79 + t74) * pkin(4) + t135;
t150 = t115 * t41 + t120 * t34;
t58 = -qJD(4) * t79 + t104 * t120 - t115 * t76;
t68 = pkin(4) * t89 - pkin(10) * t79;
t77 = t78 ^ 2;
t27 = -pkin(4) * t77 + pkin(10) * t58 - t68 * t89 + t150;
t151 = t114 * t25 + t119 * t27;
t145 = t110 * t121;
t113 = sin(qJ(6));
t118 = cos(qJ(6));
t33 = -t104 * pkin(3) - t103 * pkin(9) + t91 * t71 - t153;
t125 = -t58 * pkin(4) - t77 * pkin(10) + t79 * t68 + t33;
t62 = -t114 * t79 + t119 * t78;
t63 = t114 * t78 + t119 * t79;
t47 = -pkin(5) * t62 - pkin(11) * t63;
t72 = qJDD(5) + t74;
t84 = qJD(5) + t89;
t83 = t84 ^ 2;
t22 = -pkin(5) * t83 + pkin(11) * t72 + t47 * t62 + t151;
t37 = -qJD(5) * t63 - t114 * t59 + t119 * t58;
t38 = qJD(5) * t62 + t114 * t58 + t119 * t59;
t23 = (-t62 * t84 - t38) * pkin(11) + (t63 * t84 - t37) * pkin(5) + t125;
t49 = -t113 * t63 + t118 * t84;
t31 = qJD(6) * t49 + t113 * t72 + t118 * t38;
t36 = qJDD(6) - t37;
t50 = t113 * t84 + t118 * t63;
t42 = -mrSges(7,1) * t49 + mrSges(7,2) * t50;
t61 = qJD(6) - t62;
t43 = -mrSges(7,2) * t61 + mrSges(7,3) * t49;
t19 = m(7) * (-t113 * t22 + t118 * t23) - t31 * mrSges(7,3) + t36 * mrSges(7,1) - t50 * t42 + t61 * t43;
t30 = -qJD(6) * t50 - t113 * t38 + t118 * t72;
t44 = mrSges(7,1) * t61 - mrSges(7,3) * t50;
t20 = m(7) * (t113 * t23 + t118 * t22) + t30 * mrSges(7,3) - t36 * mrSges(7,2) + t49 * t42 - t61 * t44;
t54 = -mrSges(6,2) * t84 + mrSges(6,3) * t62;
t55 = mrSges(6,1) * t84 - mrSges(6,3) * t63;
t128 = -m(6) * t125 + t37 * mrSges(6,1) - t38 * mrSges(6,2) - t113 * t20 - t118 * t19 + t62 * t54 - t63 * t55;
t66 = -mrSges(5,2) * t89 + mrSges(5,3) * t78;
t67 = mrSges(5,1) * t89 - mrSges(5,3) * t79;
t124 = m(5) * t33 - t58 * mrSges(5,1) + t59 * mrSges(5,2) - t78 * t66 + t79 * t67 - t128;
t70 = -mrSges(4,1) * t90 + mrSges(4,2) * t91;
t80 = -mrSges(4,2) * t105 + mrSges(4,3) * t90;
t14 = m(4) * t153 + t104 * mrSges(4,1) - t76 * mrSges(4,3) + t105 * t80 - t91 * t70 - t124;
t46 = -mrSges(6,1) * t62 + mrSges(6,2) * t63;
t15 = m(6) * t151 - t72 * mrSges(6,2) + t37 * mrSges(6,3) - t113 * t19 + t118 * t20 + t62 * t46 - t84 * t55;
t131 = -t114 * t27 + t119 * t25;
t127 = m(7) * (-pkin(5) * t72 - pkin(11) * t83 + t47 * t63 - t131) - t30 * mrSges(7,1) + t31 * mrSges(7,2) - t49 * t43 + t50 * t44;
t16 = m(6) * t131 + t72 * mrSges(6,1) - t38 * mrSges(6,3) - t63 * t46 + t84 * t54 - t127;
t64 = -mrSges(5,1) * t78 + mrSges(5,2) * t79;
t12 = m(5) * t135 + t74 * mrSges(5,1) - t59 * mrSges(5,3) + t114 * t15 + t119 * t16 - t79 * t64 + t89 * t66;
t13 = m(5) * t150 - t74 * mrSges(5,2) + t58 * mrSges(5,3) - t114 * t16 + t119 * t15 + t78 * t64 - t89 * t67;
t81 = mrSges(4,1) * t105 - mrSges(4,3) * t91;
t7 = m(4) * t141 - t104 * mrSges(4,2) + t75 * mrSges(4,3) - t105 * t81 - t115 * t12 + t120 * t13 + t90 * t70;
t96 = -mrSges(3,2) * t105 + mrSges(3,3) * t138;
t99 = (-mrSges(3,1) * t121 + mrSges(3,2) * t116) * t144;
t5 = m(3) * (-g(3) * t145 + t134) - t100 * mrSges(3,3) + t104 * mrSges(3,1) - t99 * t139 + t105 * t96 + t109 * t7 + t111 * t14;
t95 = mrSges(3,1) * t105 - mrSges(3,3) * t139;
t6 = m(3) * t130 - t104 * mrSges(3,2) + t101 * mrSges(3,3) - t105 * t95 - t109 * t14 + t111 * t7 + t99 * t138;
t129 = m(4) * t126 - t75 * mrSges(4,1) + t76 * mrSges(4,2) + t115 * t13 + t120 * t12 - t90 * t80 + t91 * t81;
t9 = m(3) * t132 + t100 * mrSges(3,2) - t101 * mrSges(3,1) + (t116 * t95 - t121 * t96) * t144 + t129;
t143 = t112 * t9 + t5 * t145 + t6 * t146;
t2 = m(2) * t133 - t123 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t116 * t5 + t121 * t6;
t1 = m(2) * t137 + qJDD(1) * mrSges(2,1) - t123 * mrSges(2,2) - t110 * t9 + (t116 * t6 + t121 * t5) * t112;
t3 = [-m(1) * g(1) - t1 * t117 + t122 * t2, t2, t6, t7, t13, t15, t20; -m(1) * g(2) + t1 * t122 + t117 * t2, t1, t5, t14, t12, t16, t19; (-m(1) - m(2)) * g(3) + t143, -m(2) * g(3) + t143, t9, t129, t124, -t128, t127;];
f_new  = t3;
