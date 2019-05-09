% Calculate vector of cutting forces with Newton-Euler
% S6RRRPRR8
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
% Datum: 2019-05-07 12:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPRR8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR8_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR8_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR8_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR8_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR8_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 12:13:26
% EndTime: 2019-05-07 12:13:51
% DurationCPUTime: 9.73s
% Computational Cost: add. (180555->215), mult. (397825->296), div. (0->0), fcn. (323444->14), ass. (0->116)
t110 = sin(pkin(6));
t116 = sin(qJ(2));
t121 = cos(qJ(2));
t139 = qJD(1) * qJD(2);
t101 = (-qJDD(1) * t121 + t116 * t139) * t110;
t109 = sin(pkin(12));
t111 = cos(pkin(12));
t141 = qJD(1) * t121;
t136 = t110 * t141;
t103 = qJD(3) - t136;
t115 = sin(qJ(3));
t120 = cos(qJ(3));
t112 = cos(pkin(6));
t106 = t112 * qJD(1) + qJD(2);
t104 = t106 ^ 2;
t105 = t112 * qJDD(1) + qJDD(2);
t123 = qJD(1) ^ 2;
t117 = sin(qJ(1));
t122 = cos(qJ(1));
t135 = t117 * g(1) - t122 * g(2);
t151 = pkin(8) * t110;
t96 = qJDD(1) * pkin(1) + t123 * t151 + t135;
t146 = t112 * t96;
t131 = -t122 * g(1) - t117 * g(2);
t97 = -t123 * pkin(1) + qJDD(1) * t151 + t131;
t147 = t116 * t146 + t121 * t97;
t142 = qJD(1) * t110;
t99 = (-pkin(2) * t121 - pkin(9) * t116) * t142;
t64 = -t104 * pkin(2) + t105 * pkin(9) + (-g(3) * t116 + t99 * t141) * t110 + t147;
t100 = (qJDD(1) * t116 + t121 * t139) * t110;
t150 = t112 * g(3);
t65 = t101 * pkin(2) - t100 * pkin(9) - t150 + (-t96 + (pkin(2) * t116 - pkin(9) * t121) * t106 * qJD(1)) * t110;
t132 = -t115 * t64 + t120 * t65;
t137 = t116 * t142;
t89 = t120 * t106 - t115 * t137;
t74 = t89 * qJD(3) + t120 * t100 + t115 * t105;
t90 = t115 * t106 + t120 * t137;
t93 = qJDD(3) + t101;
t33 = (t103 * t89 - t74) * qJ(4) + (t89 * t90 + t93) * pkin(3) + t132;
t148 = t115 * t65 + t120 * t64;
t73 = -t90 * qJD(3) - t115 * t100 + t120 * t105;
t83 = t103 * pkin(3) - t90 * qJ(4);
t88 = t89 ^ 2;
t37 = -t88 * pkin(3) + t73 * qJ(4) - t103 * t83 + t148;
t80 = t109 * t89 + t111 * t90;
t152 = -0.2e1 * qJD(4) * t80 - t109 * t37 + t111 * t33;
t114 = sin(qJ(5));
t119 = cos(qJ(5));
t102 = t103 ^ 2;
t79 = -t109 * t90 + t111 * t89;
t138 = 0.2e1 * qJD(4) * t79 + t109 * t33 + t111 * t37;
t59 = -t79 * pkin(4) - t80 * pkin(10);
t25 = -t102 * pkin(4) + t93 * pkin(10) + t79 * t59 + t138;
t143 = t110 * t121;
t130 = -g(3) * t143 - t116 * t97 + t121 * t146;
t63 = -t105 * pkin(2) - t104 * pkin(9) + t99 * t137 - t130;
t126 = -t73 * pkin(3) - t88 * qJ(4) + t90 * t83 + qJDD(4) + t63;
t55 = -t109 * t74 + t111 * t73;
t56 = t109 * t73 + t111 * t74;
t28 = (-t103 * t79 - t56) * pkin(10) + (t103 * t80 - t55) * pkin(4) + t126;
t149 = t114 * t28 + t119 * t25;
t144 = t110 * t116;
t113 = sin(qJ(6));
t118 = cos(qJ(6));
t133 = -t114 * t25 + t119 * t28;
t67 = t119 * t103 - t114 * t80;
t42 = t67 * qJD(5) + t114 * t93 + t119 * t56;
t54 = qJDD(5) - t55;
t68 = t114 * t103 + t119 * t80;
t78 = qJD(5) - t79;
t19 = (t67 * t78 - t42) * pkin(11) + (t67 * t68 + t54) * pkin(5) + t133;
t41 = -t68 * qJD(5) - t114 * t56 + t119 * t93;
t51 = t78 * pkin(5) - t68 * pkin(11);
t66 = t67 ^ 2;
t20 = -t66 * pkin(5) + t41 * pkin(11) - t78 * t51 + t149;
t46 = -t113 * t68 + t118 * t67;
t31 = t46 * qJD(6) + t113 * t41 + t118 * t42;
t47 = t113 * t67 + t118 * t68;
t38 = -t46 * mrSges(7,1) + t47 * mrSges(7,2);
t75 = qJD(6) + t78;
t43 = -t75 * mrSges(7,2) + t46 * mrSges(7,3);
t52 = qJDD(6) + t54;
t17 = m(7) * (-t113 * t20 + t118 * t19) - t31 * mrSges(7,3) + t52 * mrSges(7,1) - t47 * t38 + t75 * t43;
t30 = -t47 * qJD(6) - t113 * t42 + t118 * t41;
t44 = t75 * mrSges(7,1) - t47 * mrSges(7,3);
t18 = m(7) * (t113 * t19 + t118 * t20) + t30 * mrSges(7,3) - t52 * mrSges(7,2) + t46 * t38 - t75 * t44;
t48 = -t67 * mrSges(6,1) + t68 * mrSges(6,2);
t49 = -t78 * mrSges(6,2) + t67 * mrSges(6,3);
t14 = m(6) * t133 + t54 * mrSges(6,1) - t42 * mrSges(6,3) + t113 * t18 + t118 * t17 - t68 * t48 + t78 * t49;
t50 = t78 * mrSges(6,1) - t68 * mrSges(6,3);
t15 = m(6) * t149 - t54 * mrSges(6,2) + t41 * mrSges(6,3) - t113 * t17 + t118 * t18 + t67 * t48 - t78 * t50;
t69 = -t103 * mrSges(5,2) + t79 * mrSges(5,3);
t70 = t103 * mrSges(5,1) - t80 * mrSges(5,3);
t127 = -m(5) * t126 + t55 * mrSges(5,1) - t56 * mrSges(5,2) - t114 * t15 - t119 * t14 + t79 * t69 - t80 * t70;
t82 = -t103 * mrSges(4,2) + t89 * mrSges(4,3);
t84 = t103 * mrSges(4,1) - t90 * mrSges(4,3);
t124 = m(4) * t63 - t73 * mrSges(4,1) + t74 * mrSges(4,2) - t89 * t82 + t90 * t84 - t127;
t95 = -t106 * mrSges(3,2) + mrSges(3,3) * t136;
t98 = (-mrSges(3,1) * t121 + mrSges(3,2) * t116) * t142;
t10 = m(3) * t130 + t105 * mrSges(3,1) - t100 * mrSges(3,3) + t106 * t95 - t98 * t137 - t124;
t58 = -t79 * mrSges(5,1) + t80 * mrSges(5,2);
t11 = m(5) * t138 - t93 * mrSges(5,2) + t55 * mrSges(5,3) - t103 * t70 - t114 * t14 + t119 * t15 + t79 * t58;
t24 = -t93 * pkin(4) - t102 * pkin(10) + t80 * t59 - t152;
t128 = t30 * mrSges(7,1) + t46 * t43 - m(7) * (-t41 * pkin(5) - t66 * pkin(11) + t68 * t51 + t24) - t31 * mrSges(7,2) - t47 * t44;
t125 = m(6) * t24 - t41 * mrSges(6,1) + t42 * mrSges(6,2) - t67 * t49 + t68 * t50 - t128;
t16 = m(5) * t152 + t93 * mrSges(5,1) - t56 * mrSges(5,3) + t103 * t69 - t80 * t58 - t125;
t81 = -t89 * mrSges(4,1) + t90 * mrSges(4,2);
t7 = m(4) * t132 + t93 * mrSges(4,1) - t74 * mrSges(4,3) + t103 * t82 + t109 * t11 + t111 * t16 - t90 * t81;
t8 = m(4) * t148 - t93 * mrSges(4,2) + t73 * mrSges(4,3) - t103 * t84 - t109 * t16 + t111 * t11 + t89 * t81;
t94 = t106 * mrSges(3,1) - mrSges(3,3) * t137;
t4 = m(3) * (-g(3) * t144 + t147) - t101 * mrSges(3,3) - t105 * mrSges(3,2) + t98 * t136 - t106 * t94 + t120 * t8 - t115 * t7;
t6 = m(3) * (-t110 * t96 - t150) + t100 * mrSges(3,2) + t101 * mrSges(3,1) + t115 * t8 + t120 * t7 + (t116 * t94 - t121 * t95) * t142;
t140 = t10 * t143 + t112 * t6 + t4 * t144;
t2 = m(2) * t131 - t123 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t116 * t10 + t121 * t4;
t1 = m(2) * t135 + qJDD(1) * mrSges(2,1) - t123 * mrSges(2,2) - t110 * t6 + (t121 * t10 + t116 * t4) * t112;
t3 = [-m(1) * g(1) - t117 * t1 + t122 * t2, t2, t4, t8, t11, t15, t18; -m(1) * g(2) + t122 * t1 + t117 * t2, t1, t10, t7, t16, t14, t17; (-m(1) - m(2)) * g(3) + t140, -m(2) * g(3) + t140, t6, t124, -t127, t125, -t128;];
f_new  = t3;
