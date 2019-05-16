% Calculate vector of cutting forces with Newton-Euler
% S6RRPPPR4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
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
% Datum: 2019-05-06 08:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPPPR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:41:38
% EndTime: 2019-05-06 08:41:42
% DurationCPUTime: 1.46s
% Computational Cost: add. (12351->205), mult. (27322->247), div. (0->0), fcn. (15427->8), ass. (0->97)
t106 = sin(qJ(2));
t109 = cos(qJ(2));
t107 = sin(qJ(1));
t110 = cos(qJ(1));
t133 = t107 * g(1) - t110 * g(2);
t124 = -qJDD(1) * pkin(1) - t133;
t102 = sin(pkin(9));
t103 = cos(pkin(9));
t137 = qJD(1) * qJD(2);
t131 = t109 * t137;
t132 = t106 * t137;
t139 = qJD(1) * t106;
t80 = t106 * qJDD(1) + t131;
t115 = pkin(2) * t132 - 0.2e1 * qJD(3) * t139 + (-t80 - t131) * qJ(3) + t124;
t138 = qJD(1) * t109;
t112 = qJD(1) ^ 2;
t150 = t112 * pkin(7);
t105 = sin(qJ(6));
t108 = cos(qJ(6));
t72 = t102 * qJD(2) + t103 * t138;
t134 = t72 * t139;
t101 = t109 ^ 2;
t81 = t109 * qJDD(1) - t132;
t85 = pkin(3) * t139 - qJD(2) * qJ(4);
t23 = -t85 * t139 + (-pkin(2) - qJ(4)) * t81 + (-pkin(3) * t101 - pkin(7)) * t112 + t115;
t111 = qJD(2) ^ 2;
t126 = -t110 * g(1) - t107 * g(2);
t68 = -t112 * pkin(1) + qJDD(1) * pkin(7) + t126;
t145 = -t109 * g(3) - t106 * t68;
t77 = (-pkin(2) * t109 - qJ(3) * t106) * qJD(1);
t38 = -qJDD(2) * pkin(2) - t111 * qJ(3) + t77 * t139 + qJDD(3) - t145;
t31 = (-t106 * t109 * t112 - qJDD(2)) * qJ(4) + (t80 - t131) * pkin(3) + t38;
t130 = -t102 * t23 + t103 * t31;
t140 = t106 ^ 2 * t112;
t73 = qJD(2) * t103 - t102 * t138;
t47 = pkin(4) * t72 - qJ(5) * t73;
t16 = -t80 * pkin(4) - qJ(5) * t140 + qJDD(5) - t130 + ((2 * qJD(4)) + t47) * t73;
t58 = qJDD(2) * t103 - t102 * t81;
t13 = (-t58 - t134) * pkin(8) + (t72 * t73 - t80) * pkin(5) + t16;
t153 = -2 * qJD(4);
t135 = t102 * t31 + t103 * t23 + t153 * t72;
t152 = 2 * qJD(5);
t120 = -pkin(4) * t140 + t80 * qJ(5) + t139 * t152 - t72 * t47 + t135;
t57 = qJDD(2) * t102 + t103 * t81;
t59 = -pkin(5) * t139 - t73 * pkin(8);
t70 = t72 ^ 2;
t14 = -t70 * pkin(5) + t57 * pkin(8) + t139 * t59 + t120;
t45 = -t105 * t73 + t108 * t72;
t27 = t45 * qJD(6) + t105 * t57 + t108 * t58;
t46 = t105 * t72 + t108 * t73;
t34 = -mrSges(7,1) * t45 + mrSges(7,2) * t46;
t90 = qJD(6) - t139;
t39 = -mrSges(7,2) * t90 + t45 * mrSges(7,3);
t76 = qJDD(6) - t80;
t11 = m(7) * (-t105 * t14 + t108 * t13) - t27 * mrSges(7,3) + t76 * mrSges(7,1) - t46 * t34 + t90 * t39;
t26 = -t46 * qJD(6) - t105 * t58 + t108 * t57;
t40 = mrSges(7,1) * t90 - t46 * mrSges(7,3);
t12 = m(7) * (t105 * t13 + t108 * t14) + t26 * mrSges(7,3) - t76 * mrSges(7,2) + t45 * t34 - t90 * t40;
t56 = -mrSges(6,1) * t139 + t73 * mrSges(6,2);
t123 = m(6) * t120 + t80 * mrSges(6,3) - t105 * t11 + t108 * t12 + t56 * t139;
t48 = mrSges(6,1) * t72 - mrSges(6,3) * t73;
t146 = -mrSges(5,1) * t72 - mrSges(5,2) * t73 - t48;
t147 = -mrSges(5,3) - mrSges(6,2);
t55 = mrSges(5,1) * t139 - t73 * mrSges(5,3);
t6 = m(5) * t135 - t80 * mrSges(5,2) - t139 * t55 + t146 * t72 + t147 * t57 + t123;
t121 = -m(6) * t16 - t105 * t12 - t108 * t11;
t53 = -t72 * mrSges(6,2) + mrSges(6,3) * t139;
t54 = -mrSges(5,2) * t139 - t72 * mrSges(5,3);
t7 = m(5) * t130 + (mrSges(5,1) + mrSges(6,1)) * t80 + (m(5) * t153 + t146) * t73 + t147 * t58 + (t53 + t54) * t139 + t121;
t86 = -mrSges(4,1) * t138 - qJD(2) * mrSges(4,3);
t125 = t102 * t7 - t103 * t6 - m(4) * (-t81 * pkin(2) + t115 - t150) - t86 * t138 + t80 * mrSges(4,3);
t87 = mrSges(4,1) * t139 + qJD(2) * mrSges(4,2);
t142 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t139 - t87;
t149 = mrSges(3,1) - mrSges(4,2);
t84 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t138;
t157 = (t106 * t142 - t109 * t84) * qJD(1) - t149 * t81 + m(3) * (t124 - t150) + t80 * mrSges(3,2) - t125;
t156 = (-t58 + t134) * qJ(5);
t122 = -m(4) * t38 - t102 * t6 - t103 * t7;
t78 = (mrSges(4,2) * t109 - mrSges(4,3) * t106) * qJD(1);
t143 = t78 + (-mrSges(3,1) * t109 + mrSges(3,2) * t106) * qJD(1);
t148 = mrSges(3,3) + mrSges(4,1);
t4 = m(3) * t145 - t148 * t80 + t149 * qJDD(2) + (t84 - t86) * qJD(2) - t143 * t139 + t122;
t144 = -t106 * g(3) + t109 * t68;
t128 = t111 * pkin(2) - qJDD(2) * qJ(3) - t77 * t138 - t144;
t117 = t101 * t112 * qJ(4) - t81 * pkin(3) - qJD(2) * t85 - qJDD(4) + t128;
t136 = qJD(3) * qJD(2);
t116 = -t117 + 0.2e1 * t136;
t95 = -0.2e1 * t136;
t129 = m(7) * (-t70 * pkin(8) + t95 + (-pkin(4) - pkin(5)) * t57 - t156 + (-pkin(4) * t139 + t152 + t59) * t73 + t117) + t27 * mrSges(7,2) - t26 * mrSges(7,1) + t46 * t40 - t45 * t39;
t118 = -t58 * mrSges(6,3) - t73 * t56 - t129 + m(6) * (-0.2e1 * qJD(5) * t73 + t156 + (t139 * t73 + t57) * pkin(4) + t116) + t72 * t53 + t57 * mrSges(6,1);
t114 = m(5) * t116 + t57 * mrSges(5,1) + t58 * mrSges(5,2) + t72 * t54 + t73 * t55 + t118;
t113 = -m(4) * (t95 + t128) + t114;
t9 = t113 + t148 * t81 + (-mrSges(3,2) + mrSges(4,3)) * qJDD(2) - t142 * qJD(2) + m(3) * t144 + t143 * t138;
t151 = t106 * t9 + t109 * t4;
t2 = m(2) * t133 + qJDD(1) * mrSges(2,1) - t112 * mrSges(2,2) - t157;
t1 = m(2) * t126 - t112 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t106 * t4 + t109 * t9;
t3 = [-m(1) * g(1) + t110 * t1 - t107 * t2, t1, t9, t81 * mrSges(4,2) - t139 * t87 - t125, t6, -t57 * mrSges(6,2) - t72 * t48 + t123, t12; -m(1) * g(2) + t107 * t1 + t110 * t2, t2, t4, -t81 * mrSges(4,1) - qJDD(2) * mrSges(4,3) - qJD(2) * t87 - t138 * t78 - t113, t7, t118, t11; (-m(1) - m(2)) * g(3) + t151, -m(2) * g(3) + t151, t157, t80 * mrSges(4,1) + qJDD(2) * mrSges(4,2) + qJD(2) * t86 + t139 * t78 - t122, t114, -t80 * mrSges(6,1) + t58 * mrSges(6,2) - t139 * t53 + t73 * t48 - t121, t129;];
f_new  = t3;
