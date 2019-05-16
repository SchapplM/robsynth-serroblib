% Calculate vector of cutting forces with Newton-Euler
% S6RRRPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-05-07 14:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPRR10_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR10_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR10_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR10_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR10_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR10_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR10_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 13:58:09
% EndTime: 2019-05-07 13:58:18
% DurationCPUTime: 3.25s
% Computational Cost: add. (43282->203), mult. (86239->253), div. (0->0), fcn. (58986->10), ass. (0->103)
t113 = sin(qJ(3));
t114 = sin(qJ(2));
t118 = cos(qJ(2));
t151 = cos(qJ(3));
t106 = t118 * qJD(1);
t102 = -t106 + qJD(3);
t112 = sin(qJ(5));
t117 = cos(qJ(5));
t101 = t102 ^ 2;
t140 = qJD(1) * qJD(2);
t137 = t118 * t140;
t138 = t114 * t140;
t121 = qJD(1) ^ 2;
t115 = sin(qJ(1));
t119 = cos(qJ(1));
t136 = t115 * g(1) - t119 * g(2);
t83 = -qJDD(1) * pkin(1) - t121 * pkin(7) - t136;
t94 = t114 * qJDD(1) + t137;
t95 = t118 * qJDD(1) - t138;
t48 = (-t94 - t137) * pkin(8) + (-t95 + t138) * pkin(2) + t83;
t120 = qJD(2) ^ 2;
t132 = -t119 * g(1) - t115 * g(2);
t84 = -t121 * pkin(1) + qJDD(1) * pkin(7) + t132;
t139 = -t114 * g(3) + t118 * t84;
t93 = (-pkin(2) * t118 - pkin(8) * t114) * qJD(1);
t52 = -t120 * pkin(2) + qJDD(2) * pkin(8) + t93 * t106 + t139;
t146 = t113 * t48 + t151 * t52;
t153 = 2 * qJD(4);
t141 = qJD(1) * t114;
t90 = -t151 * qJD(2) + t113 * t141;
t91 = t113 * qJD(2) + t151 * t141;
t69 = t90 * pkin(3) - t91 * qJ(4);
t89 = qJDD(3) - t95;
t129 = -t101 * pkin(3) + t89 * qJ(4) + t102 * t153 - t90 * t69 + t146;
t74 = -t102 * mrSges(5,1) + t91 * mrSges(5,2);
t100 = qJD(5) - t102;
t111 = sin(qJ(6));
t116 = cos(qJ(6));
t143 = t102 * t90;
t133 = -t113 * t52 + t151 * t48;
t32 = -t89 * pkin(3) - t101 * qJ(4) + t91 * t69 + qJDD(4) - t133;
t64 = -t90 * qJD(3) + t113 * qJDD(2) + t151 * t94;
t24 = (-t64 - t143) * pkin(9) + (t90 * t91 - t89) * pkin(4) + t32;
t63 = t91 * qJD(3) - t151 * qJDD(2) + t113 * t94;
t76 = -t102 * pkin(4) - t91 * pkin(9);
t88 = t90 ^ 2;
t27 = -t88 * pkin(4) + t63 * pkin(9) + t102 * t76 + t129;
t135 = -t112 * t27 + t117 * t24;
t66 = -t112 * t91 + t117 * t90;
t39 = t66 * qJD(5) + t112 * t63 + t117 * t64;
t67 = t112 * t90 + t117 * t91;
t87 = qJDD(5) - t89;
t14 = (t100 * t66 - t39) * pkin(10) + (t66 * t67 + t87) * pkin(5) + t135;
t147 = t112 * t24 + t117 * t27;
t38 = -t67 * qJD(5) - t112 * t64 + t117 * t63;
t55 = t100 * pkin(5) - t67 * pkin(10);
t65 = t66 ^ 2;
t15 = -t65 * pkin(5) + t38 * pkin(10) - t100 * t55 + t147;
t44 = -t111 * t67 + t116 * t66;
t21 = t44 * qJD(6) + t111 * t38 + t116 * t39;
t45 = t111 * t66 + t116 * t67;
t35 = -t44 * mrSges(7,1) + t45 * mrSges(7,2);
t96 = qJD(6) + t100;
t40 = -t96 * mrSges(7,2) + t44 * mrSges(7,3);
t81 = qJDD(6) + t87;
t12 = m(7) * (-t111 * t15 + t116 * t14) - t21 * mrSges(7,3) + t81 * mrSges(7,1) - t45 * t35 + t96 * t40;
t20 = -t45 * qJD(6) - t111 * t39 + t116 * t38;
t41 = t96 * mrSges(7,1) - t45 * mrSges(7,3);
t13 = m(7) * (t111 * t14 + t116 * t15) + t20 * mrSges(7,3) - t81 * mrSges(7,2) + t44 * t35 - t96 * t41;
t46 = -t66 * mrSges(6,1) + t67 * mrSges(6,2);
t53 = -t100 * mrSges(6,2) + t66 * mrSges(6,3);
t8 = m(6) * t135 + t87 * mrSges(6,1) - t39 * mrSges(6,3) + t100 * t53 + t111 * t13 + t116 * t12 - t67 * t46;
t54 = t100 * mrSges(6,1) - t67 * mrSges(6,3);
t9 = m(6) * t147 - t87 * mrSges(6,2) + t38 * mrSges(6,3) - t100 * t54 - t111 * t12 + t116 * t13 + t66 * t46;
t130 = m(5) * t129 + t89 * mrSges(5,3) + t102 * t74 - t112 * t8 + t117 * t9;
t70 = t90 * mrSges(5,1) - t91 * mrSges(5,3);
t145 = -t90 * mrSges(4,1) - t91 * mrSges(4,2) - t70;
t148 = -mrSges(4,3) - mrSges(5,2);
t73 = t102 * mrSges(4,1) - t91 * mrSges(4,3);
t5 = m(4) * t146 - t89 * mrSges(4,2) - t102 * t73 + t145 * t90 + t148 * t63 + t130;
t128 = -m(5) * t32 - t112 * t9 - t117 * t8;
t72 = -t102 * mrSges(4,2) - t90 * mrSges(4,3);
t75 = -t90 * mrSges(5,2) + t102 * mrSges(5,3);
t6 = m(4) * t133 + t145 * t91 + (mrSges(4,1) + mrSges(5,1)) * t89 + t148 * t64 + (t72 + t75) * t102 + t128;
t97 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t141;
t98 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t106;
t155 = m(3) * t83 - t95 * mrSges(3,1) + t94 * mrSges(3,2) + (t114 * t97 - t118 * t98) * qJD(1) + t113 * t5 + t151 * t6;
t142 = -t118 * g(3) - t114 * t84;
t51 = -qJDD(2) * pkin(2) - t120 * pkin(8) + t93 * t141 - t142;
t127 = t63 * pkin(3) + t51 + (t143 - t64) * qJ(4);
t150 = pkin(3) * t102;
t123 = -t63 * pkin(4) - t88 * pkin(9) - t127 + (-t150 + t153 + t76) * t91;
t134 = m(7) * (-t38 * pkin(5) - t65 * pkin(10) + t67 * t55 + t123) + t21 * mrSges(7,2) - t20 * mrSges(7,1) + t45 * t41 - t44 * t40;
t126 = m(6) * t123 - t38 * mrSges(6,1) + t39 * mrSges(6,2) - t66 * t53 + t67 * t54 + t134;
t125 = m(5) * ((-(2 * qJD(4)) + t150) * t91 + t127) + t63 * mrSges(5,1) + t90 * t75 - t126;
t154 = m(4) * t51 + t63 * mrSges(4,1) + (t73 - t74) * t91 + (mrSges(4,2) - mrSges(5,3)) * t64 + t90 * t72 + t125;
t92 = (-mrSges(3,1) * t118 + mrSges(3,2) * t114) * qJD(1);
t11 = m(3) * t142 + qJDD(2) * mrSges(3,1) - t94 * mrSges(3,3) + qJD(2) * t98 - t92 * t141 - t154;
t4 = m(3) * t139 - qJDD(2) * mrSges(3,2) + t95 * mrSges(3,3) - qJD(2) * t97 + t92 * t106 - t113 * t6 + t151 * t5;
t152 = t118 * t11 + t114 * t4;
t2 = m(2) * t136 + qJDD(1) * mrSges(2,1) - t121 * mrSges(2,2) - t155;
t1 = m(2) * t132 - t121 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t114 * t11 + t118 * t4;
t3 = [-m(1) * g(1) + t119 * t1 - t115 * t2, t1, t4, t5, -t63 * mrSges(5,2) - t90 * t70 + t130, t9, t13; -m(1) * g(2) + t115 * t1 + t119 * t2, t2, t11, t6, -t64 * mrSges(5,3) - t91 * t74 + t125, t8, t12; (-m(1) - m(2)) * g(3) + t152, -m(2) * g(3) + t152, t155, t154, -t89 * mrSges(5,1) + t64 * mrSges(5,2) - t102 * t75 + t91 * t70 - t128, t126, t134;];
f_new  = t3;
