% Calculate vector of cutting forces with Newton-Euler
% S6RRRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
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
% Datum: 2019-05-07 23:04
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRPR10_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR10_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR10_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR10_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR10_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR10_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR10_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 22:49:14
% EndTime: 2019-05-07 22:49:32
% DurationCPUTime: 5.46s
% Computational Cost: add. (73053->213), mult. (156820->278), div. (0->0), fcn. (124942->12), ass. (0->111)
t105 = cos(pkin(6));
t100 = t105 * qJDD(1) + qJDD(2);
t109 = sin(qJ(2));
t113 = cos(qJ(2));
t104 = sin(pkin(6));
t138 = t104 * t113;
t115 = qJD(1) ^ 2;
t110 = sin(qJ(1));
t114 = cos(qJ(1));
t129 = t110 * g(1) - t114 * g(2);
t88 = t115 * t104 * pkin(8) + qJDD(1) * pkin(1) + t129;
t140 = t105 * t88;
t126 = -t114 * g(1) - t110 * g(2);
t134 = qJDD(1) * t104;
t89 = -t115 * pkin(1) + pkin(8) * t134 + t126;
t125 = -g(3) * t138 - t109 * t89 + t113 * t140;
t137 = qJD(1) * t104;
t131 = t109 * t137;
t91 = (-pkin(2) * t113 - pkin(9) * t109) * t137;
t101 = t105 * qJD(1) + qJD(2);
t99 = t101 ^ 2;
t52 = -t100 * pkin(2) - t99 * pkin(9) + t91 * t131 - t125;
t108 = sin(qJ(3));
t112 = cos(qJ(3));
t81 = t108 * t101 + t112 * t131;
t136 = qJD(1) * t113;
t92 = (qJD(2) * t136 + qJDD(1) * t109) * t104;
t65 = -t81 * qJD(3) + t112 * t100 - t108 * t92;
t130 = t104 * t136;
t97 = qJD(3) - t130;
t75 = t97 * pkin(3) - t81 * pkin(10);
t80 = t112 * t101 - t108 * t131;
t79 = t80 ^ 2;
t118 = -t65 * pkin(3) - t79 * pkin(10) + t81 * t75 + t52;
t106 = sin(qJ(6));
t111 = cos(qJ(6));
t107 = sin(qJ(4));
t150 = cos(qJ(4));
t70 = t107 * t81 - t150 * t80;
t95 = -qJD(4) - t97;
t148 = t70 * t95;
t151 = -2 * qJD(5);
t66 = t80 * qJD(3) + t108 * t100 + t112 * t92;
t39 = -t70 * qJD(4) + t107 * t65 + t150 * t66;
t71 = t107 * t80 + t150 * t81;
t116 = (-t39 - t148) * qJ(5) + t118 + (-t95 * pkin(4) + t151) * t71;
t141 = t109 * t140 + t113 * t89;
t53 = -t99 * pkin(2) + t100 * pkin(9) + (-g(3) * t109 + t136 * t91) * t104 + t141;
t149 = t105 * g(3);
t93 = -qJD(2) * t131 + t113 * t134;
t54 = -t93 * pkin(2) - t92 * pkin(9) - t149 + (-t88 + (pkin(2) * t109 - pkin(9) * t113) * t101 * qJD(1)) * t104;
t128 = -t108 * t53 + t112 * t54;
t85 = qJDD(3) - t93;
t25 = (t80 * t97 - t66) * pkin(10) + (t80 * t81 + t85) * pkin(3) + t128;
t143 = t108 * t54 + t112 * t53;
t28 = -t79 * pkin(3) + t65 * pkin(10) - t97 * t75 + t143;
t127 = -t107 * t28 + t150 * t25;
t46 = t70 * pkin(4) - t71 * qJ(5);
t84 = qJDD(4) + t85;
t94 = t95 ^ 2;
t21 = -t84 * pkin(4) - t94 * qJ(5) + t71 * t46 + qJDD(5) - t127;
t16 = (t70 * t71 - t84) * pkin(11) + (t39 - t148) * pkin(5) + t21;
t38 = t71 * qJD(4) + t107 * t66 - t150 * t65;
t61 = t71 * pkin(5) + t95 * pkin(11);
t69 = t70 ^ 2;
t19 = (pkin(4) + pkin(11)) * t38 + t116 - t69 * pkin(5) - t71 * t61;
t55 = t106 * t95 + t111 * t70;
t31 = t55 * qJD(6) + t106 * t38 + t111 * t84;
t37 = qJDD(6) + t39;
t56 = t106 * t70 - t111 * t95;
t41 = -t55 * mrSges(7,1) + t56 * mrSges(7,2);
t68 = qJD(6) + t71;
t42 = -t68 * mrSges(7,2) + t55 * mrSges(7,3);
t14 = m(7) * (-t106 * t19 + t111 * t16) - t31 * mrSges(7,3) + t37 * mrSges(7,1) - t56 * t41 + t68 * t42;
t30 = -t56 * qJD(6) - t106 * t84 + t111 * t38;
t43 = t68 * mrSges(7,1) - t56 * mrSges(7,3);
t15 = m(7) * (t106 * t16 + t111 * t19) + t30 * mrSges(7,3) - t37 * mrSges(7,2) + t55 * t41 - t68 * t43;
t58 = t71 * mrSges(6,1) - t95 * mrSges(6,2);
t124 = t106 * t14 - t111 * t15 - m(6) * (t38 * pkin(4) + t116) + t39 * mrSges(6,3) + t71 * t58;
t57 = t70 * mrSges(6,1) + t95 * mrSges(6,3);
t142 = -t95 * mrSges(5,2) + t70 * mrSges(5,3) + t57;
t147 = mrSges(5,1) - mrSges(6,2);
t60 = -t95 * mrSges(5,1) - t71 * mrSges(5,3);
t154 = m(5) * t118 + t39 * mrSges(5,2) - t142 * t70 + t147 * t38 + t71 * t60 - t124;
t73 = -t97 * mrSges(4,2) + t80 * mrSges(4,3);
t74 = t97 * mrSges(4,1) - t81 * mrSges(4,3);
t153 = m(4) * t52 - t65 * mrSges(4,1) + t66 * mrSges(4,2) - t80 * t73 + t81 * t74 + t154;
t146 = -mrSges(5,3) - mrSges(6,1);
t145 = t107 * t25 + t150 * t28;
t48 = -t70 * mrSges(6,2) - t71 * mrSges(6,3);
t144 = -t70 * mrSges(5,1) - t71 * mrSges(5,2) - t48;
t139 = t104 * t109;
t87 = -t101 * mrSges(3,2) + mrSges(3,3) * t130;
t90 = (-mrSges(3,1) * t113 + mrSges(3,2) * t109) * t137;
t10 = m(3) * t125 + t100 * mrSges(3,1) - t92 * mrSges(3,3) + t101 * t87 - t90 * t131 - t153;
t123 = -m(6) * t21 - t106 * t15 - t111 * t14;
t11 = m(5) * t127 + t142 * t95 + t144 * t71 + t146 * t39 + t147 * t84 + t123;
t120 = -t94 * pkin(4) + t84 * qJ(5) - t70 * t46 + t145;
t122 = -t30 * mrSges(7,1) - t55 * t42 + m(7) * (-t38 * pkin(5) - t69 * pkin(11) + (t151 - t61) * t95 + t120) + t31 * mrSges(7,2) + t56 * t43;
t119 = -m(6) * (0.2e1 * qJD(5) * t95 - t120) + t122;
t12 = m(5) * t145 + (t60 - t58) * t95 + (-mrSges(5,2) + mrSges(6,3)) * t84 + t144 * t70 + t146 * t38 + t119;
t72 = -t80 * mrSges(4,1) + t81 * mrSges(4,2);
t7 = m(4) * t128 + t85 * mrSges(4,1) - t66 * mrSges(4,3) + t107 * t12 + t150 * t11 - t81 * t72 + t97 * t73;
t8 = m(4) * t143 - t85 * mrSges(4,2) + t65 * mrSges(4,3) - t107 * t11 + t150 * t12 + t80 * t72 - t97 * t74;
t86 = t101 * mrSges(3,1) - mrSges(3,3) * t131;
t4 = m(3) * (-g(3) * t139 + t141) + t93 * mrSges(3,3) - t100 * mrSges(3,2) + t90 * t130 - t101 * t86 + t112 * t8 - t108 * t7;
t6 = m(3) * (-t104 * t88 - t149) + t92 * mrSges(3,2) - t93 * mrSges(3,1) + t108 * t8 + t112 * t7 + (t109 * t86 - t113 * t87) * t137;
t135 = t10 * t138 + t105 * t6 + t4 * t139;
t2 = m(2) * t126 - t115 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t109 * t10 + t113 * t4;
t1 = m(2) * t129 + qJDD(1) * mrSges(2,1) - t115 * mrSges(2,2) - t104 * t6 + (t113 * t10 + t109 * t4) * t105;
t3 = [-m(1) * g(1) - t110 * t1 + t114 * t2, t2, t4, t8, t12, -t38 * mrSges(6,2) - t70 * t57 - t124, t15; -m(1) * g(2) + t114 * t1 + t110 * t2, t1, t10, t7, t11, t38 * mrSges(6,1) - t84 * mrSges(6,3) + t70 * t48 + t95 * t58 - t119, t14; (-m(1) - m(2)) * g(3) + t135, -m(2) * g(3) + t135, t6, t153, t154, t39 * mrSges(6,1) + t84 * mrSges(6,2) + t71 * t48 - t95 * t57 - t123, t122;];
f_new  = t3;
