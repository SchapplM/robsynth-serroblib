% Calculate vector of cutting forces with Newton-Euler
% S6RRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-05-06 13:41
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRPR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 13:30:30
% EndTime: 2019-05-06 13:30:51
% DurationCPUTime: 8.64s
% Computational Cost: add. (134232->213), mult. (356694->297), div. (0->0), fcn. (284719->14), ass. (0->116)
t153 = -2 * qJD(3);
t107 = sin(pkin(11));
t110 = cos(pkin(11));
t114 = sin(qJ(2));
t118 = cos(qJ(2));
t108 = sin(pkin(6));
t143 = qJD(1) * t108;
t89 = (t107 * t114 - t110 * t118) * t143;
t111 = cos(pkin(6));
t101 = qJDD(1) * t111 + qJDD(2);
t102 = qJD(1) * t111 + qJD(2);
t120 = qJD(1) ^ 2;
t115 = sin(qJ(1));
t119 = cos(qJ(1));
t135 = g(1) * t115 - g(2) * t119;
t150 = pkin(8) * t108;
t95 = qJDD(1) * pkin(1) + t120 * t150 + t135;
t148 = t111 * t95;
t131 = -g(1) * t119 - t115 * g(2);
t96 = -pkin(1) * t120 + qJDD(1) * t150 + t131;
t132 = -t114 * t96 + t118 * t148;
t146 = t108 ^ 2 * t120;
t141 = qJD(1) * qJD(2);
t98 = (qJDD(1) * t114 + t118 * t141) * t108;
t53 = t101 * pkin(2) - t98 * qJ(3) + (pkin(2) * t114 * t146 + (qJ(3) * qJD(1) * t102 - g(3)) * t108) * t118 + t132;
t145 = t108 * t114;
t128 = -g(3) * t145 + t114 * t148 + t118 * t96;
t138 = t118 ^ 2 * t146;
t137 = t114 * t143;
t92 = pkin(2) * t102 - qJ(3) * t137;
t99 = (qJDD(1) * t118 - t114 * t141) * t108;
t56 = -pkin(2) * t138 + qJ(3) * t99 - t102 * t92 + t128;
t90 = (t107 * t118 + t110 * t114) * t143;
t152 = -t107 * t56 + t110 * t53 + t153 * t90;
t151 = 2 * qJD(5);
t113 = sin(qJ(4));
t117 = cos(qJ(4));
t100 = t102 ^ 2;
t139 = t107 * t53 + t110 * t56 + t153 * t89;
t72 = pkin(3) * t89 - pkin(9) * t90;
t34 = -pkin(3) * t100 + pkin(9) * t101 - t72 * t89 + t139;
t130 = -t111 * g(3) - t108 * t95;
t123 = -t99 * pkin(2) - qJ(3) * t138 + t137 * t92 + qJDD(3) + t130;
t75 = -t107 * t98 + t110 * t99;
t76 = t107 * t99 + t110 * t98;
t37 = (t102 * t89 - t76) * pkin(9) + (t102 * t90 - t75) * pkin(3) + t123;
t149 = t113 * t37 + t117 * t34;
t144 = t108 * t118;
t112 = sin(qJ(6));
t116 = cos(qJ(6));
t33 = -t101 * pkin(3) - t100 * pkin(9) + t72 * t90 - t152;
t79 = t102 * t113 + t117 * t90;
t58 = -qJD(4) * t79 + t101 * t117 - t113 * t76;
t88 = qJD(4) + t89;
t68 = pkin(4) * t88 - qJ(5) * t79;
t78 = t102 * t117 - t113 * t90;
t77 = t78 ^ 2;
t122 = -t58 * pkin(4) - t77 * qJ(5) + t68 * t79 + qJDD(5) + t33;
t106 = sin(pkin(12));
t109 = cos(pkin(12));
t133 = -t113 * t34 + t117 * t37;
t59 = qJD(4) * t78 + t101 * t113 + t117 * t76;
t74 = qJDD(4) - t75;
t25 = (t78 * t88 - t59) * qJ(5) + (t78 * t79 + t74) * pkin(4) + t133;
t27 = -pkin(4) * t77 + t58 * qJ(5) - t68 * t88 + t149;
t63 = -t106 * t79 + t109 * t78;
t140 = t106 * t25 + t109 * t27 + t151 * t63;
t64 = t106 * t78 + t109 * t79;
t47 = -pkin(5) * t63 - pkin(10) * t64;
t87 = t88 ^ 2;
t22 = -pkin(5) * t87 + pkin(10) * t74 + t63 * t47 + t140;
t41 = -t106 * t59 + t109 * t58;
t42 = t106 * t58 + t109 * t59;
t23 = (-t63 * t88 - t42) * pkin(10) + (t64 * t88 - t41) * pkin(5) + t122;
t49 = -t112 * t64 + t116 * t88;
t31 = t49 * qJD(6) + t112 * t74 + t116 * t42;
t50 = t112 * t88 + t116 * t64;
t38 = -mrSges(7,1) * t49 + mrSges(7,2) * t50;
t40 = qJDD(6) - t41;
t62 = qJD(6) - t63;
t43 = -mrSges(7,2) * t62 + mrSges(7,3) * t49;
t19 = m(7) * (-t112 * t22 + t116 * t23) - t31 * mrSges(7,3) + t40 * mrSges(7,1) - t50 * t38 + t62 * t43;
t30 = -t50 * qJD(6) - t112 * t42 + t116 * t74;
t44 = mrSges(7,1) * t62 - mrSges(7,3) * t50;
t20 = m(7) * (t112 * t23 + t116 * t22) + t30 * mrSges(7,3) - t40 * mrSges(7,2) + t49 * t38 - t62 * t44;
t54 = -mrSges(6,2) * t88 + t63 * mrSges(6,3);
t55 = mrSges(6,1) * t88 - t64 * mrSges(6,3);
t125 = -m(6) * t122 + t41 * mrSges(6,1) - t42 * mrSges(6,2) - t112 * t20 - t116 * t19 + t63 * t54 - t64 * t55;
t67 = -mrSges(5,2) * t88 + mrSges(5,3) * t78;
t69 = mrSges(5,1) * t88 - mrSges(5,3) * t79;
t121 = m(5) * t33 - t58 * mrSges(5,1) + t59 * mrSges(5,2) - t78 * t67 + t79 * t69 - t125;
t71 = mrSges(4,1) * t89 + mrSges(4,2) * t90;
t80 = -mrSges(4,2) * t102 - mrSges(4,3) * t89;
t14 = m(4) * t152 + t101 * mrSges(4,1) - t76 * mrSges(4,3) + t102 * t80 - t90 * t71 - t121;
t46 = -mrSges(6,1) * t63 + mrSges(6,2) * t64;
t15 = m(6) * t140 - mrSges(6,2) * t74 + t41 * mrSges(6,3) - t112 * t19 + t116 * t20 + t63 * t46 - t55 * t88;
t129 = t106 * t27 - t109 * t25;
t124 = m(7) * (-t74 * pkin(5) - t87 * pkin(10) + (t151 + t47) * t64 + t129) - t30 * mrSges(7,1) + t31 * mrSges(7,2) - t49 * t43 + t50 * t44;
t16 = m(6) * (-0.2e1 * qJD(5) * t64 - t129) - t42 * mrSges(6,3) + t74 * mrSges(6,1) - t64 * t46 + t88 * t54 - t124;
t65 = -mrSges(5,1) * t78 + mrSges(5,2) * t79;
t12 = m(5) * t133 + mrSges(5,1) * t74 - t59 * mrSges(5,3) + t106 * t15 + t109 * t16 - t65 * t79 + t67 * t88;
t13 = m(5) * t149 - mrSges(5,2) * t74 + t58 * mrSges(5,3) - t106 * t16 + t109 * t15 + t65 * t78 - t69 * t88;
t81 = mrSges(4,1) * t102 - mrSges(4,3) * t90;
t7 = m(4) * t139 - mrSges(4,2) * t101 + mrSges(4,3) * t75 - t102 * t81 - t113 * t12 + t117 * t13 - t71 * t89;
t136 = t118 * t143;
t94 = -mrSges(3,2) * t102 + mrSges(3,3) * t136;
t97 = (-mrSges(3,1) * t118 + mrSges(3,2) * t114) * t143;
t5 = m(3) * (-g(3) * t144 + t132) - t98 * mrSges(3,3) + t101 * mrSges(3,1) - t97 * t137 + t102 * t94 + t107 * t7 + t110 * t14;
t93 = mrSges(3,1) * t102 - mrSges(3,3) * t137;
t6 = m(3) * t128 - t101 * mrSges(3,2) + t99 * mrSges(3,3) - t102 * t93 - t107 * t14 + t110 * t7 + t136 * t97;
t126 = m(4) * t123 - t75 * mrSges(4,1) + mrSges(4,2) * t76 + t113 * t13 + t117 * t12 + t89 * t80 + t81 * t90;
t9 = m(3) * t130 + t98 * mrSges(3,2) - t99 * mrSges(3,1) + (t114 * t93 - t118 * t94) * t143 + t126;
t142 = t111 * t9 + t144 * t5 + t145 * t6;
t2 = m(2) * t131 - mrSges(2,1) * t120 - qJDD(1) * mrSges(2,2) - t114 * t5 + t118 * t6;
t1 = m(2) * t135 + qJDD(1) * mrSges(2,1) - t120 * mrSges(2,2) - t108 * t9 + (t114 * t6 + t118 * t5) * t111;
t3 = [-m(1) * g(1) - t1 * t115 + t119 * t2, t2, t6, t7, t13, t15, t20; -m(1) * g(2) + t1 * t119 + t115 * t2, t1, t5, t14, t12, t16, t19; (-m(1) - m(2)) * g(3) + t142, -m(2) * g(3) + t142, t9, t126, t121, -t125, t124;];
f_new  = t3;
