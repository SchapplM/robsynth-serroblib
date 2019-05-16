% Calculate vector of cutting forces with Newton-Euler
% S6RRPRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-05-06 15:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRPR10_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR10_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR10_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR10_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR10_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR10_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR10_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 15:37:08
% EndTime: 2019-05-06 15:37:24
% DurationCPUTime: 5.25s
% Computational Cost: add. (63927->212), mult. (146080->280), div. (0->0), fcn. (115527->12), ass. (0->109)
t107 = cos(pkin(6));
t100 = qJDD(1) * t107 + qJDD(2);
t110 = sin(qJ(2));
t113 = cos(qJ(2));
t105 = sin(pkin(6));
t139 = t105 * t113;
t115 = qJD(1) ^ 2;
t111 = sin(qJ(1));
t114 = cos(qJ(1));
t129 = t111 * g(1) - g(2) * t114;
t89 = pkin(8) * t105 * t115 + qJDD(1) * pkin(1) + t129;
t141 = t107 * t89;
t126 = -g(1) * t114 - g(2) * t111;
t136 = qJDD(1) * t105;
t90 = -pkin(1) * t115 + pkin(8) * t136 + t126;
t125 = -g(3) * t139 - t110 * t90 + t113 * t141;
t138 = qJD(1) * t105;
t131 = t110 * t138;
t91 = (-pkin(2) * t113 - qJ(3) * t110) * t138;
t101 = qJD(1) * t107 + qJD(2);
t99 = t101 ^ 2;
t52 = -t100 * pkin(2) - t99 * qJ(3) + t91 * t131 + qJDD(3) - t125;
t104 = sin(pkin(11));
t106 = cos(pkin(11));
t137 = qJD(1) * t113;
t93 = (qJD(2) * t137 + qJDD(1) * t110) * t105;
t73 = t100 * t106 - t104 * t93;
t130 = t105 * t137;
t83 = t101 * t104 + t106 * t131;
t75 = -pkin(3) * t130 - pkin(9) * t83;
t82 = t101 * t106 - t104 * t131;
t81 = t82 ^ 2;
t118 = -t73 * pkin(3) - t81 * pkin(9) + t83 * t75 + t52;
t108 = sin(qJ(6));
t112 = cos(qJ(6));
t109 = sin(qJ(4));
t150 = cos(qJ(4));
t67 = t109 * t83 - t150 * t82;
t96 = -qJD(4) + t130;
t148 = t67 * t96;
t151 = -2 * qJD(5);
t74 = t100 * t104 + t106 * t93;
t41 = -t67 * qJD(4) + t109 * t73 + t150 * t74;
t68 = t109 * t82 + t150 * t83;
t116 = (-t41 - t148) * qJ(5) + t118 + (-t96 * pkin(4) + t151) * t68;
t142 = t110 * t141 + t113 * t90;
t53 = -t99 * pkin(2) + t100 * qJ(3) + (-g(3) * t110 + t91 * t137) * t105 + t142;
t149 = t107 * g(3);
t94 = -qJD(2) * t131 + t113 * t136;
t54 = -t94 * pkin(2) - t149 - t93 * qJ(3) + (-t89 + (pkin(2) * t110 - qJ(3) * t113) * t101 * qJD(1)) * t105;
t128 = -0.2e1 * qJD(3) * t83 - t104 * t53 + t106 * t54;
t25 = (-t82 * t130 - t74) * pkin(9) + (t82 * t83 - t94) * pkin(3) + t128;
t134 = 0.2e1 * qJD(3) * t82 + t104 * t54 + t106 * t53;
t28 = -pkin(3) * t81 + pkin(9) * t73 + t75 * t130 + t134;
t127 = -t109 * t28 + t150 * t25;
t46 = pkin(4) * t67 - qJ(5) * t68;
t86 = qJDD(4) - t94;
t95 = t96 ^ 2;
t21 = -t86 * pkin(4) - t95 * qJ(5) + t68 * t46 + qJDD(5) - t127;
t16 = (t67 * t68 - t86) * pkin(10) + (t41 - t148) * pkin(5) + t21;
t40 = t68 * qJD(4) + t109 * t74 - t150 * t73;
t61 = t68 * pkin(5) + pkin(10) * t96;
t66 = t67 ^ 2;
t19 = t116 + (pkin(4) + pkin(10)) * t40 - t68 * t61 - t66 * pkin(5);
t55 = t108 * t96 + t112 * t67;
t31 = t55 * qJD(6) + t108 * t40 + t112 * t86;
t56 = t108 * t67 - t112 * t96;
t36 = -mrSges(7,1) * t55 + mrSges(7,2) * t56;
t39 = qJDD(6) + t41;
t65 = qJD(6) + t68;
t42 = -mrSges(7,2) * t65 + mrSges(7,3) * t55;
t14 = m(7) * (-t108 * t19 + t112 * t16) - t31 * mrSges(7,3) + t39 * mrSges(7,1) - t56 * t36 + t65 * t42;
t30 = -t56 * qJD(6) - t108 * t86 + t112 * t40;
t43 = mrSges(7,1) * t65 - mrSges(7,3) * t56;
t15 = m(7) * (t108 * t16 + t112 * t19) + t30 * mrSges(7,3) - t39 * mrSges(7,2) + t55 * t36 - t65 * t43;
t58 = t68 * mrSges(6,1) - mrSges(6,2) * t96;
t124 = t108 * t14 - t112 * t15 - m(6) * (t40 * pkin(4) + t116) + t41 * mrSges(6,3) + t68 * t58;
t57 = t67 * mrSges(6,1) + mrSges(6,3) * t96;
t143 = -mrSges(5,2) * t96 + t67 * mrSges(5,3) + t57;
t147 = mrSges(5,1) - mrSges(6,2);
t60 = -mrSges(5,1) * t96 - t68 * mrSges(5,3);
t154 = m(5) * t118 + t41 * mrSges(5,2) - t143 * t67 + t147 * t40 + t68 * t60 - t124;
t71 = mrSges(4,2) * t130 + mrSges(4,3) * t82;
t72 = -mrSges(4,1) * t130 - mrSges(4,3) * t83;
t153 = m(4) * t52 - t73 * mrSges(4,1) + t74 * mrSges(4,2) - t82 * t71 + t83 * t72 + t154;
t146 = -mrSges(5,3) - mrSges(6,1);
t145 = t109 * t25 + t150 * t28;
t48 = -mrSges(6,2) * t67 - mrSges(6,3) * t68;
t144 = -mrSges(5,1) * t67 - mrSges(5,2) * t68 - t48;
t140 = t105 * t110;
t88 = -mrSges(3,2) * t101 + mrSges(3,3) * t130;
t92 = (-mrSges(3,1) * t113 + mrSges(3,2) * t110) * t138;
t11 = m(3) * t125 + t100 * mrSges(3,1) - t93 * mrSges(3,3) + t101 * t88 - t92 * t131 - t153;
t120 = -t95 * pkin(4) + t86 * qJ(5) - t67 * t46 + t145;
t122 = -t30 * mrSges(7,1) - t55 * t42 + m(7) * (-t40 * pkin(5) - t66 * pkin(10) + (t151 - t61) * t96 + t120) + t31 * mrSges(7,2) + t56 * t43;
t119 = -m(6) * (0.2e1 * qJD(5) * t96 - t120) + t122;
t12 = m(5) * t145 + (t60 - t58) * t96 + (-mrSges(5,2) + mrSges(6,3)) * t86 + t144 * t67 + t146 * t40 + t119;
t69 = -mrSges(4,1) * t82 + mrSges(4,2) * t83;
t123 = -m(6) * t21 - t108 * t15 - t112 * t14;
t9 = m(5) * t127 + t143 * t96 + t144 * t68 + t146 * t41 + t147 * t86 + t123;
t7 = m(4) * t128 - t94 * mrSges(4,1) - t74 * mrSges(4,3) + t109 * t12 - t71 * t130 + t150 * t9 - t83 * t69;
t8 = m(4) * t134 + t94 * mrSges(4,2) + t73 * mrSges(4,3) - t109 * t9 + t150 * t12 + t72 * t130 + t82 * t69;
t87 = mrSges(3,1) * t101 - mrSges(3,3) * t131;
t4 = m(3) * (-g(3) * t140 + t142) + t94 * mrSges(3,3) - t100 * mrSges(3,2) + t92 * t130 - t101 * t87 + t106 * t8 - t104 * t7;
t6 = m(3) * (-t105 * t89 - t149) + t93 * mrSges(3,2) - t94 * mrSges(3,1) + t104 * t8 + t106 * t7 + (t110 * t87 - t113 * t88) * t138;
t135 = t107 * t6 + t11 * t139 + t4 * t140;
t2 = m(2) * t126 - t115 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t110 * t11 + t113 * t4;
t1 = m(2) * t129 + qJDD(1) * mrSges(2,1) - t115 * mrSges(2,2) - t105 * t6 + (t11 * t113 + t110 * t4) * t107;
t3 = [-m(1) * g(1) - t1 * t111 + t114 * t2, t2, t4, t8, t12, -t40 * mrSges(6,2) - t67 * t57 - t124, t15; -m(1) * g(2) + t1 * t114 + t111 * t2, t1, t11, t7, t9, t40 * mrSges(6,1) - t86 * mrSges(6,3) + t67 * t48 + t96 * t58 - t119, t14; (-m(1) - m(2)) * g(3) + t135, -m(2) * g(3) + t135, t6, t153, t154, t41 * mrSges(6,1) + t86 * mrSges(6,2) + t68 * t48 - t96 * t57 - t123, t122;];
f_new  = t3;
