% Calculate vector of cutting forces with Newton-Euler
% S6RPRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-05-05 16:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRPPR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:47:46
% EndTime: 2019-05-05 16:47:51
% DurationCPUTime: 2.46s
% Computational Cost: add. (27383->193), mult. (67848->240), div. (0->0), fcn. (49636->10), ass. (0->99)
t100 = sin(pkin(9));
t101 = cos(pkin(9));
t103 = sin(qJ(3));
t146 = cos(qJ(3));
t153 = t100 * t103 - t101 * t146;
t108 = qJD(1) ^ 2;
t98 = t101 ^ 2;
t137 = t100 ^ 2 + t98;
t152 = t137 * mrSges(3,3);
t107 = qJD(3) ^ 2;
t131 = qJD(1) * qJD(2);
t127 = -t101 * g(3) - 0.2e1 * t100 * t131;
t135 = pkin(7) * qJDD(1);
t145 = pkin(2) * t108;
t104 = sin(qJ(1));
t106 = cos(qJ(1));
t122 = -t106 * g(1) - t104 * g(2);
t88 = -t108 * pkin(1) + qJDD(1) * qJ(2) + t122;
t53 = (t101 * t145 - t135 - t88) * t100 + t127;
t125 = -t100 * g(3) + (0.2e1 * t131 + t88) * t101;
t62 = t101 * t135 - t98 * t145 + t125;
t140 = -t103 * t62 + t146 * t53;
t86 = t153 * qJD(1);
t117 = t146 * t100 + t101 * t103;
t87 = t117 * qJD(1);
t66 = t86 * pkin(3) - t87 * qJ(4);
t112 = qJDD(3) * pkin(3) + t107 * qJ(4) - t87 * t66 - qJDD(4) + t140;
t136 = cos(pkin(10));
t99 = sin(pkin(10));
t77 = -t136 * qJD(3) + t99 * t87;
t144 = t77 * t86;
t133 = t86 * qJD(3);
t73 = t117 * qJDD(1) - t133;
t61 = t99 * qJDD(3) + t136 * t73;
t151 = (-t61 + t144) * qJ(5) - t112;
t148 = 2 * qJD(5);
t102 = sin(qJ(6));
t105 = cos(qJ(6));
t78 = t99 * qJD(3) + t136 * t87;
t44 = t102 * t77 + t105 * t78;
t60 = -t136 * qJDD(3) + t99 * t73;
t29 = -t44 * qJD(6) - t102 * t61 + t105 * t60;
t43 = -t102 * t78 + t105 * t77;
t30 = t43 * qJD(6) + t102 * t60 + t105 * t61;
t83 = qJD(6) - t86;
t39 = -t83 * mrSges(7,2) + t43 * mrSges(7,3);
t40 = t83 * mrSges(7,1) - t44 * mrSges(7,3);
t59 = -t86 * pkin(5) - t78 * pkin(8);
t76 = t77 ^ 2;
t124 = m(7) * (-t76 * pkin(8) + (-pkin(4) - pkin(5)) * t60 + (-pkin(4) * t86 + t148 + t59) * t78 - t151) + t30 * mrSges(7,2) - t29 * mrSges(7,1) + t44 * t40 - t43 * t39;
t55 = -t77 * mrSges(6,2) + t86 * mrSges(6,3);
t116 = m(6) * (-0.2e1 * qJD(5) * t78 + (t78 * t86 + t60) * pkin(4) + t151) + t77 * t55 + t60 * mrSges(6,1) - t124;
t56 = -t86 * mrSges(5,2) - t77 * mrSges(5,3);
t57 = t86 * mrSges(5,1) - t78 * mrSges(5,3);
t58 = -t86 * mrSges(6,1) + t78 * mrSges(6,2);
t150 = -m(5) * t112 + t60 * mrSges(5,1) + (t57 - t58) * t78 + (mrSges(5,2) - mrSges(6,3)) * t61 + t77 * t56 + t116;
t149 = -2 * qJD(4);
t120 = -t101 * mrSges(3,1) + t100 * mrSges(3,2);
t119 = qJDD(1) * mrSges(3,3) + t108 * t120;
t67 = t86 * mrSges(4,1) + t87 * mrSges(4,2);
t79 = -qJD(3) * mrSges(4,2) - t86 * mrSges(4,3);
t12 = m(4) * t140 + qJDD(3) * mrSges(4,1) - t73 * mrSges(4,3) + qJD(3) * t79 - t87 * t67 - t150;
t139 = t103 * t53 + t146 * t62;
t31 = -t107 * pkin(3) + qJDD(3) * qJ(4) - t86 * t66 + t139;
t128 = t104 * g(1) - t106 * g(2);
t123 = qJDD(2) - t128;
t109 = (-pkin(2) * t101 - pkin(1)) * qJDD(1) + (-t137 * pkin(7) - qJ(2)) * t108 + t123;
t132 = t87 * qJD(3);
t72 = t153 * qJDD(1) + t132;
t33 = (-t73 + t133) * qJ(4) + (t72 + t132) * pkin(3) + t109;
t121 = t136 * t33 - t99 * t31;
t45 = t77 * pkin(4) - t78 * qJ(5);
t85 = t86 ^ 2;
t21 = -t72 * pkin(4) - t85 * qJ(5) + qJDD(5) - t121 + ((2 * qJD(4)) + t45) * t78;
t16 = (-t61 - t144) * pkin(8) + (t77 * t78 - t72) * pkin(5) + t21;
t130 = t136 * t31 + t77 * t149 + t99 * t33;
t114 = -t85 * pkin(4) + t72 * qJ(5) + t86 * t148 - t77 * t45 + t130;
t17 = -t76 * pkin(5) + t60 * pkin(8) + t86 * t59 + t114;
t36 = -t43 * mrSges(7,1) + t44 * mrSges(7,2);
t70 = qJDD(6) - t72;
t14 = m(7) * (-t102 * t17 + t105 * t16) - t30 * mrSges(7,3) + t70 * mrSges(7,1) - t44 * t36 + t83 * t39;
t15 = m(7) * (t102 * t16 + t105 * t17) + t29 * mrSges(7,3) - t70 * mrSges(7,2) + t43 * t36 - t83 * t40;
t115 = -m(6) * t21 - t102 * t15 - t105 * t14;
t46 = t77 * mrSges(6,1) - t78 * mrSges(6,3);
t141 = -t77 * mrSges(5,1) - t78 * mrSges(5,2) - t46;
t142 = -mrSges(5,3) - mrSges(6,2);
t11 = m(5) * t121 + (t56 + t55) * t86 + (m(5) * t149 + t141) * t78 + (mrSges(5,1) + mrSges(6,1)) * t72 + t142 * t61 + t115;
t80 = qJD(3) * mrSges(4,1) - t87 * mrSges(4,3);
t118 = m(6) * t114 + t72 * mrSges(6,3) - t102 * t14 + t105 * t15 + t86 * t58;
t9 = m(5) * t130 - t72 * mrSges(5,2) + t141 * t77 + t142 * t60 - t86 * t57 + t118;
t7 = m(4) * t139 - qJDD(3) * mrSges(4,2) - t72 * mrSges(4,3) - qJD(3) * t80 - t99 * t11 + t136 * t9 - t86 * t67;
t4 = m(3) * t127 + t103 * t7 + t146 * t12 + (-m(3) * t88 - t119) * t100;
t5 = m(3) * t125 + t119 * t101 - t103 * t12 + t146 * t7;
t147 = t100 * t5 + t101 * t4;
t113 = m(4) * t109 + t72 * mrSges(4,1) + t73 * mrSges(4,2) + t136 * t11 + t86 * t79 + t87 * t80 + t99 * t9;
t111 = m(3) * (-qJDD(1) * pkin(1) - t108 * qJ(2) + t123) + t113;
t6 = m(2) * t128 + (-mrSges(2,2) + t152) * t108 + (mrSges(2,1) - t120) * qJDD(1) - t111;
t1 = m(2) * t122 - t108 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t100 * t4 + t101 * t5;
t2 = [-m(1) * g(1) + t106 * t1 - t104 * t6, t1, t5, t7, t9, -t60 * mrSges(6,2) - t77 * t46 + t118, t15; -m(1) * g(2) + t104 * t1 + t106 * t6, t6, t4, t12, t11, -t61 * mrSges(6,3) - t78 * t58 + t116, t14; (-m(1) - m(2)) * g(3) + t147, -m(2) * g(3) + t147, t120 * qJDD(1) - t108 * t152 + t111, t113, t150, -t72 * mrSges(6,1) + t61 * mrSges(6,2) + t78 * t46 - t86 * t55 - t115, t124;];
f_new  = t2;
