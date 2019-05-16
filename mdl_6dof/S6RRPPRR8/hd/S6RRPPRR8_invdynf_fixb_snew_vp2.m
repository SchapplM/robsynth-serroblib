% Calculate vector of cutting forces with Newton-Euler
% S6RRPPRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-05-06 11:23
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPPRR8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR8_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR8_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR8_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR8_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR8_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 11:17:42
% EndTime: 2019-05-06 11:17:49
% DurationCPUTime: 3.08s
% Computational Cost: add. (37507->201), mult. (82705->254), div. (0->0), fcn. (55521->10), ass. (0->101)
t156 = -2 * qJD(4);
t117 = cos(qJ(2));
t104 = t117 * qJD(1);
t110 = sin(pkin(10));
t113 = sin(qJ(2));
t141 = qJD(1) * t113;
t144 = cos(pkin(10));
t88 = -t144 * qJD(2) + t110 * t141;
t138 = t88 * t104;
t119 = qJD(2) ^ 2;
t120 = qJD(1) ^ 2;
t114 = sin(qJ(1));
t118 = cos(qJ(1));
t132 = -t118 * g(1) - t114 * g(2);
t84 = -t120 * pkin(1) + qJDD(1) * pkin(7) + t132;
t145 = -t117 * g(3) - t113 * t84;
t92 = (-pkin(2) * t117 - qJ(3) * t113) * qJD(1);
t51 = -qJDD(2) * pkin(2) - t119 * qJ(3) + t92 * t141 + qJDD(3) - t145;
t140 = qJD(1) * qJD(2);
t136 = t117 * t140;
t94 = t113 * qJDD(1) + t136;
t73 = -t144 * qJDD(2) + t110 * t94;
t74 = t110 * qJDD(2) + t144 * t94;
t89 = t110 * qJD(2) + t144 * t141;
t155 = t51 + (-t138 - t74) * qJ(4) + (-t104 * t89 + t73) * pkin(3) + t89 * t156;
t112 = sin(qJ(5));
t116 = cos(qJ(5));
t152 = -2 * qJD(3);
t101 = t113 * t140;
t135 = t114 * g(1) - t118 * g(2);
t83 = -qJDD(1) * pkin(1) - t120 * pkin(7) - t135;
t95 = t117 * qJDD(1) - t101;
t48 = (-t94 - t136) * qJ(3) + (-t95 + t101) * pkin(2) + t83;
t137 = -t113 * g(3) + t117 * t84;
t52 = -t119 * pkin(2) + qJDD(2) * qJ(3) + t92 * t104 + t137;
t139 = t110 * t48 + t144 * t52 + t88 * t152;
t142 = t117 ^ 2 * t120;
t63 = t88 * pkin(3) - t89 * qJ(4);
t126 = -pkin(3) * t142 - t95 * qJ(4) + t104 * t156 - t88 * t63 + t139;
t111 = sin(qJ(6));
t115 = cos(qJ(6));
t131 = -t110 * t52 + t144 * t48;
t30 = t95 * pkin(3) - qJ(4) * t142 + qJDD(4) - t131 + ((2 * qJD(3)) + t63) * t89;
t20 = (-t74 + t138) * pkin(8) + (t88 * t89 + t95) * pkin(4) + t30;
t75 = pkin(4) * t104 - t89 * pkin(8);
t86 = t88 ^ 2;
t22 = -t86 * pkin(4) + t73 * pkin(8) - t75 * t104 + t126;
t134 = -t112 * t22 + t116 * t20;
t61 = -t112 * t89 + t116 * t88;
t41 = t61 * qJD(5) + t112 * t73 + t116 * t74;
t62 = t112 * t88 + t116 * t89;
t91 = qJDD(5) + t95;
t99 = t104 + qJD(5);
t14 = (t61 * t99 - t41) * pkin(9) + (t61 * t62 + t91) * pkin(5) + t134;
t148 = t112 * t20 + t116 * t22;
t40 = -t62 * qJD(5) - t112 * t74 + t116 * t73;
t57 = t99 * pkin(5) - t62 * pkin(9);
t60 = t61 ^ 2;
t15 = -t60 * pkin(5) + t40 * pkin(9) - t99 * t57 + t148;
t44 = -t111 * t62 + t115 * t61;
t26 = t44 * qJD(6) + t111 * t40 + t115 * t41;
t45 = t111 * t61 + t115 * t62;
t33 = -t44 * mrSges(7,1) + t45 * mrSges(7,2);
t98 = qJD(6) + t99;
t36 = -t98 * mrSges(7,2) + t44 * mrSges(7,3);
t85 = qJDD(6) + t91;
t12 = m(7) * (-t111 * t15 + t115 * t14) - t26 * mrSges(7,3) + t85 * mrSges(7,1) - t45 * t33 + t98 * t36;
t25 = -t45 * qJD(6) - t111 * t41 + t115 * t40;
t37 = t98 * mrSges(7,1) - t45 * mrSges(7,3);
t13 = m(7) * (t111 * t14 + t115 * t15) + t25 * mrSges(7,3) - t85 * mrSges(7,2) + t44 * t33 - t98 * t37;
t46 = -t61 * mrSges(6,1) + t62 * mrSges(6,2);
t53 = -t99 * mrSges(6,2) + t61 * mrSges(6,3);
t8 = m(6) * t134 + t91 * mrSges(6,1) - t41 * mrSges(6,3) + t111 * t13 + t115 * t12 - t62 * t46 + t99 * t53;
t54 = t99 * mrSges(6,1) - t62 * mrSges(6,3);
t9 = m(6) * t148 - t91 * mrSges(6,2) + t40 * mrSges(6,3) - t111 * t12 + t115 * t13 + t61 * t46 - t99 * t54;
t130 = m(5) * t126 - t95 * mrSges(5,3) - t112 * t8 + t116 * t9;
t72 = mrSges(5,1) * t104 + t89 * mrSges(5,2);
t146 = -mrSges(4,1) * t104 - t89 * mrSges(4,3) - t72;
t64 = t88 * mrSges(5,1) - t89 * mrSges(5,3);
t147 = -t88 * mrSges(4,1) - t89 * mrSges(4,2) - t64;
t149 = -mrSges(4,3) - mrSges(5,2);
t5 = m(4) * t139 + t95 * mrSges(4,2) + t146 * t104 + t147 * t88 + t149 * t73 + t130;
t128 = -m(5) * t30 - t112 * t9 - t116 * t8;
t69 = -t88 * mrSges(5,2) - mrSges(5,3) * t104;
t70 = mrSges(4,2) * t104 - t88 * mrSges(4,3);
t6 = m(4) * t131 + (-mrSges(4,1) - mrSges(5,1)) * t95 + (m(4) * t152 + t147) * t89 + t149 * t74 + (-t69 - t70) * t104 + t128;
t96 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t141;
t97 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t104;
t154 = m(3) * t83 - t95 * mrSges(3,1) + t94 * mrSges(3,2) + (t113 * t96 - t117 * t97) * qJD(1) + t110 * t5 + t144 * t6;
t121 = -t73 * pkin(4) - t86 * pkin(8) + t89 * t75 - t155;
t133 = m(7) * (-t40 * pkin(5) - t60 * pkin(9) + t62 * t57 + t121) + t26 * mrSges(7,2) - t25 * mrSges(7,1) + t45 * t37 - t44 * t36;
t127 = m(6) * t121 - t40 * mrSges(6,1) + t41 * mrSges(6,2) - t61 * t53 + t62 * t54 + t133;
t125 = m(5) * t155 + t73 * mrSges(5,1) + t88 * t69 - t127;
t153 = m(4) * t51 + t73 * mrSges(4,1) + t146 * t89 + (mrSges(4,2) - mrSges(5,3)) * t74 + t88 * t70 + t125;
t93 = (-mrSges(3,1) * t117 + mrSges(3,2) * t113) * qJD(1);
t11 = m(3) * t145 + qJDD(2) * mrSges(3,1) - t94 * mrSges(3,3) + qJD(2) * t97 - t93 * t141 - t153;
t4 = m(3) * t137 - qJDD(2) * mrSges(3,2) + t95 * mrSges(3,3) - qJD(2) * t96 + t93 * t104 - t110 * t6 + t144 * t5;
t151 = t117 * t11 + t113 * t4;
t2 = m(2) * t135 + qJDD(1) * mrSges(2,1) - t120 * mrSges(2,2) - t154;
t1 = m(2) * t132 - t120 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t113 * t11 + t117 * t4;
t3 = [-m(1) * g(1) + t118 * t1 - t114 * t2, t1, t4, t5, -t73 * mrSges(5,2) - t72 * t104 - t88 * t64 + t130, t9, t13; -m(1) * g(2) + t114 * t1 + t118 * t2, t2, t11, t6, -t74 * mrSges(5,3) - t89 * t72 + t125, t8, t12; (-m(1) - m(2)) * g(3) + t151, -m(2) * g(3) + t151, t154, t153, t95 * mrSges(5,1) + t74 * mrSges(5,2) + t69 * t104 + t89 * t64 - t128, t127, t133;];
f_new  = t3;
