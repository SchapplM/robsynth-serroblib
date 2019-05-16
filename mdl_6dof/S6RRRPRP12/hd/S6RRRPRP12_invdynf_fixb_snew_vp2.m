% Calculate vector of cutting forces with Newton-Euler
% S6RRRPRP12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-05-07 09:38
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPRP12_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP12_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP12_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP12_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP12_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP12_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP12_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 09:23:25
% EndTime: 2019-05-07 09:23:34
% DurationCPUTime: 2.61s
% Computational Cost: add. (32183->207), mult. (68996->265), div. (0->0), fcn. (52255->10), ass. (0->102)
t101 = sin(pkin(6));
t105 = sin(qJ(2));
t107 = cos(qJ(2));
t128 = qJD(1) * qJD(2);
t89 = (-qJDD(1) * t107 + t105 * t128) * t101;
t104 = sin(qJ(3));
t149 = cos(qJ(3));
t130 = qJD(1) * t107;
t102 = cos(pkin(6));
t109 = qJD(1) ^ 2;
t106 = sin(qJ(1));
t108 = cos(qJ(1));
t124 = g(1) * t106 - g(2) * t108;
t147 = pkin(8) * t101;
t84 = qJDD(1) * pkin(1) + t109 * t147 + t124;
t135 = t102 * t84;
t122 = -g(1) * t108 - g(2) * t106;
t85 = -pkin(1) * t109 + qJDD(1) * t147 + t122;
t136 = t105 * t135 + t107 * t85;
t131 = qJD(1) * t101;
t87 = (-pkin(2) * t107 - pkin(9) * t105) * t131;
t98 = qJD(1) * t102 + qJD(2);
t96 = t98 ^ 2;
t97 = qJDD(1) * t102 + qJDD(2);
t38 = -t96 * pkin(2) + t97 * pkin(9) + (-g(3) * t105 + t130 * t87) * t101 + t136;
t146 = t102 * g(3);
t88 = (qJDD(1) * t105 + t107 * t128) * t101;
t39 = t89 * pkin(2) - t88 * pkin(9) - t146 + (-t84 + (pkin(2) * t105 - pkin(9) * t107) * t98 * qJD(1)) * t101;
t140 = t104 * t39 + t149 * t38;
t126 = t105 * t131;
t76 = t104 * t126 - t149 * t98;
t77 = t104 * t98 + t126 * t149;
t57 = pkin(3) * t76 - qJ(4) * t77;
t81 = qJDD(3) + t89;
t125 = t101 * t130;
t94 = -qJD(3) + t125;
t93 = t94 ^ 2;
t152 = t93 * pkin(3) - qJ(4) * t81 + 0.2e1 * qJD(4) * t94 + t57 * t76 - t140;
t103 = sin(qJ(5));
t145 = t76 * t94;
t132 = t101 * t107;
t121 = -g(3) * t132 - t105 * t85 + t107 * t135;
t37 = -t97 * pkin(2) - t96 * pkin(9) + t126 * t87 - t121;
t56 = -qJD(3) * t76 + t104 * t97 + t149 * t88;
t110 = (-t56 - t145) * qJ(4) + t37 + (-pkin(3) * t94 - 0.2e1 * qJD(4)) * t77;
t148 = cos(qJ(5));
t123 = -t104 * t38 + t149 * t39;
t25 = -t81 * pkin(3) - t93 * qJ(4) + t57 * t77 + qJDD(4) - t123;
t19 = (t76 * t77 - t81) * pkin(10) + (t56 - t145) * pkin(4) + t25;
t55 = qJD(3) * t77 + t104 * t88 - t149 * t97;
t67 = pkin(4) * t77 + pkin(10) * t94;
t75 = t76 ^ 2;
t23 = -t75 * pkin(4) - t77 * t67 + (pkin(3) + pkin(10)) * t55 + t110;
t141 = t103 * t19 + t148 * t23;
t61 = -t103 * t94 - t148 * t76;
t62 = t103 * t76 - t148 * t94;
t41 = pkin(5) * t61 - qJ(6) * t62;
t74 = qJD(5) + t77;
t49 = -mrSges(7,1) * t74 + mrSges(7,2) * t62;
t53 = qJDD(5) + t56;
t73 = t74 ^ 2;
t127 = m(7) * (-pkin(5) * t73 + qJ(6) * t53 + 0.2e1 * qJD(6) * t74 - t41 * t61 + t141) + t74 * t49 + t53 * mrSges(7,3);
t42 = mrSges(7,1) * t61 - mrSges(7,3) * t62;
t139 = -mrSges(6,1) * t61 - mrSges(6,2) * t62 - t42;
t142 = -mrSges(6,3) - mrSges(7,2);
t30 = qJD(5) * t62 + t103 * t81 - t148 * t55;
t48 = mrSges(6,1) * t74 - mrSges(6,3) * t62;
t12 = m(6) * t141 - t53 * mrSges(6,2) + t139 * t61 + t142 * t30 - t74 * t48 + t127;
t118 = -t103 * t23 + t148 * t19;
t150 = m(7) * (-pkin(5) * t53 - qJ(6) * t73 + t41 * t62 + qJDD(6) - t118);
t31 = -qJD(5) * t61 + t103 * t55 + t148 * t81;
t46 = -mrSges(7,2) * t61 + mrSges(7,3) * t74;
t47 = -mrSges(6,2) * t74 - mrSges(6,3) * t61;
t13 = m(6) * t118 - t150 + (t47 + t46) * t74 + t139 * t62 + (mrSges(6,1) + mrSges(7,1)) * t53 + t142 * t31;
t66 = mrSges(5,1) * t77 - mrSges(5,2) * t94;
t119 = t103 * t13 - t148 * t12 - m(5) * (t55 * pkin(3) + t110) + t56 * mrSges(5,3) + t77 * t66;
t65 = mrSges(5,1) * t76 + mrSges(5,3) * t94;
t137 = mrSges(4,2) * t94 - mrSges(4,3) * t76 - t65;
t144 = mrSges(4,1) - mrSges(5,2);
t64 = -mrSges(4,1) * t94 - mrSges(4,3) * t77;
t151 = m(4) * t37 + t56 * mrSges(4,2) + t137 * t76 + t144 * t55 + t77 * t64 - t119;
t143 = -mrSges(4,3) - mrSges(5,1);
t59 = -mrSges(5,2) * t76 - mrSges(5,3) * t77;
t138 = -mrSges(4,1) * t76 - mrSges(4,2) * t77 - t59;
t133 = t101 * t105;
t114 = -t55 * pkin(4) - t75 * pkin(10) - t67 * t94 - t152;
t116 = -t31 * mrSges(7,3) - t62 * t49 + m(7) * (-0.2e1 * qJD(6) * t62 + (t61 * t74 - t31) * qJ(6) + (t62 * t74 + t30) * pkin(5) + t114) + t30 * mrSges(7,1) + t61 * t46;
t112 = m(6) * t114 + t30 * mrSges(6,1) + mrSges(6,2) * t31 + t61 * t47 + t48 * t62 + t116;
t111 = -m(5) * t152 + t112;
t10 = (t64 - t66) * t94 + (-mrSges(4,2) + mrSges(5,3)) * t81 + t138 * t76 + t143 * t55 + m(4) * t140 + t111;
t82 = mrSges(3,1) * t98 - mrSges(3,3) * t126;
t86 = (-mrSges(3,1) * t107 + mrSges(3,2) * t105) * t131;
t115 = m(5) * t25 + t103 * t12 + t13 * t148;
t9 = m(4) * t123 - t137 * t94 + t138 * t77 + t143 * t56 + t144 * t81 - t115;
t4 = m(3) * (-g(3) * t133 + t136) - t89 * mrSges(3,3) - t97 * mrSges(3,2) + t86 * t125 - t98 * t82 + t149 * t10 - t104 * t9;
t83 = -mrSges(3,2) * t98 + mrSges(3,3) * t125;
t6 = m(3) * (-t101 * t84 - t146) + t88 * mrSges(3,2) + t89 * mrSges(3,1) + t104 * t10 + t149 * t9 + (t105 * t82 - t107 * t83) * t131;
t8 = m(3) * t121 + t97 * mrSges(3,1) - t88 * mrSges(3,3) - t126 * t86 + t98 * t83 - t151;
t129 = t102 * t6 + t132 * t8 + t133 * t4;
t2 = m(2) * t122 - mrSges(2,1) * t109 - qJDD(1) * mrSges(2,2) - t105 * t8 + t107 * t4;
t1 = m(2) * t124 + qJDD(1) * mrSges(2,1) - t109 * mrSges(2,2) - t101 * t6 + (t105 * t4 + t107 * t8) * t102;
t3 = [-m(1) * g(1) - t1 * t106 + t108 * t2, t2, t4, t10, -t55 * mrSges(5,2) - t76 * t65 - t119, t12, -t30 * mrSges(7,2) - t61 * t42 + t127; -m(1) * g(2) + t1 * t108 + t106 * t2, t1, t8, t9, t55 * mrSges(5,1) - t81 * mrSges(5,3) + t76 * t59 + t94 * t66 - t111, t13, t116; (-m(1) - m(2)) * g(3) + t129, -m(2) * g(3) + t129, t6, t151, t56 * mrSges(5,1) + t81 * mrSges(5,2) + t77 * t59 - t94 * t65 + t115, t112, -t53 * mrSges(7,1) + t31 * mrSges(7,2) + t62 * t42 - t74 * t46 + t150;];
f_new  = t3;
