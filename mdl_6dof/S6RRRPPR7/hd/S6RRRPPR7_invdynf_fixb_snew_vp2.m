% Calculate vector of cutting forces with Newton-Euler
% S6RRRPPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-05-07 05:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPPR7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR7_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR7_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR7_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR7_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR7_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 05:51:51
% EndTime: 2019-05-07 05:52:00
% DurationCPUTime: 3.20s
% Computational Cost: add. (41301->202), mult. (84045->255), div. (0->0), fcn. (56810->10), ass. (0->101)
t114 = sin(qJ(3));
t115 = sin(qJ(2));
t118 = cos(qJ(2));
t151 = cos(qJ(3));
t141 = t118 * qJD(1);
t102 = -qJD(3) + t141;
t110 = sin(pkin(10));
t111 = cos(pkin(10));
t101 = t102 ^ 2;
t140 = qJD(1) * qJD(2);
t104 = t115 * t140;
t137 = t118 * t140;
t121 = qJD(1) ^ 2;
t116 = sin(qJ(1));
t119 = cos(qJ(1));
t136 = t116 * g(1) - t119 * g(2);
t84 = -qJDD(1) * pkin(1) - t121 * pkin(7) - t136;
t95 = t115 * qJDD(1) + t137;
t96 = t118 * qJDD(1) - t104;
t48 = (-t95 - t137) * pkin(8) + (-t96 + t104) * pkin(2) + t84;
t120 = qJD(2) ^ 2;
t132 = -t119 * g(1) - t116 * g(2);
t85 = -t121 * pkin(1) + qJDD(1) * pkin(7) + t132;
t138 = -t115 * g(3) + t118 * t85;
t94 = (-pkin(2) * t118 - pkin(8) * t115) * qJD(1);
t52 = -t120 * pkin(2) + qJDD(2) * pkin(8) + t94 * t141 + t138;
t147 = t114 * t48 + t151 * t52;
t153 = -2 * qJD(4);
t142 = qJD(1) * t115;
t91 = -t151 * qJD(2) + t114 * t142;
t92 = t114 * qJD(2) + t151 * t142;
t71 = t91 * pkin(3) - t92 * qJ(4);
t90 = -qJDD(3) + t96;
t129 = -t101 * pkin(3) - t90 * qJ(4) + t102 * t153 - t91 * t71 + t147;
t77 = t102 * mrSges(5,1) + t92 * mrSges(5,2);
t113 = sin(qJ(6));
t117 = cos(qJ(6));
t100 = qJD(6) + t102;
t144 = t102 * t91;
t133 = -t114 * t52 + t151 * t48;
t32 = t90 * pkin(3) - t101 * qJ(4) + t92 * t71 + qJDD(4) - t133;
t66 = -t91 * qJD(3) + t114 * qJDD(2) + t151 * t95;
t20 = (-t66 + t144) * qJ(5) + (t91 * t92 + t90) * pkin(4) + t32;
t65 = t92 * qJD(3) - t151 * qJDD(2) + t114 * t95;
t75 = t102 * pkin(4) - t92 * qJ(5);
t89 = t91 ^ 2;
t24 = -t89 * pkin(4) + t65 * qJ(5) - t102 * t75 + t129;
t69 = t110 * t91 + t111 * t92;
t135 = -0.2e1 * qJD(5) * t69 - t110 * t24 + t111 * t20;
t41 = t110 * t65 + t111 * t66;
t68 = -t110 * t92 + t111 * t91;
t14 = (t102 * t68 - t41) * pkin(9) + (t68 * t69 + t90) * pkin(5) + t135;
t139 = 0.2e1 * qJD(5) * t68 + t110 * t20 + t111 * t24;
t40 = -t110 * t66 + t111 * t65;
t55 = t102 * pkin(5) - t69 * pkin(9);
t67 = t68 ^ 2;
t15 = -t67 * pkin(5) + t40 * pkin(9) - t102 * t55 + t139;
t44 = -t113 * t69 + t117 * t68;
t28 = t44 * qJD(6) + t113 * t40 + t117 * t41;
t45 = t113 * t68 + t117 * t69;
t35 = -t44 * mrSges(7,1) + t45 * mrSges(7,2);
t36 = -t100 * mrSges(7,2) + t44 * mrSges(7,3);
t88 = qJDD(6) + t90;
t12 = m(7) * (-t113 * t15 + t117 * t14) - t28 * mrSges(7,3) + t88 * mrSges(7,1) - t45 * t35 + t100 * t36;
t27 = -t45 * qJD(6) - t113 * t41 + t117 * t40;
t37 = t100 * mrSges(7,1) - t45 * mrSges(7,3);
t13 = m(7) * (t113 * t14 + t117 * t15) + t27 * mrSges(7,3) - t88 * mrSges(7,2) + t44 * t35 - t100 * t37;
t46 = -t68 * mrSges(6,1) + t69 * mrSges(6,2);
t53 = -t102 * mrSges(6,2) + t68 * mrSges(6,3);
t8 = m(6) * t135 + t90 * mrSges(6,1) - t41 * mrSges(6,3) + t102 * t53 + t113 * t13 + t117 * t12 - t69 * t46;
t54 = t102 * mrSges(6,1) - t69 * mrSges(6,3);
t9 = m(6) * t139 - t90 * mrSges(6,2) + t40 * mrSges(6,3) - t102 * t54 - t113 * t12 + t117 * t13 + t68 * t46;
t130 = m(5) * t129 - t90 * mrSges(5,3) - t102 * t77 - t110 * t8 + t111 * t9;
t72 = t91 * mrSges(5,1) - t92 * mrSges(5,3);
t146 = -t91 * mrSges(4,1) - t92 * mrSges(4,2) - t72;
t148 = -mrSges(4,3) - mrSges(5,2);
t76 = -t102 * mrSges(4,1) - t92 * mrSges(4,3);
t5 = m(4) * t147 + t90 * mrSges(4,2) + t102 * t76 + t146 * t91 + t148 * t65 + t130;
t128 = -m(5) * t32 - t110 * t9 - t111 * t8;
t74 = t102 * mrSges(4,2) - t91 * mrSges(4,3);
t78 = -t91 * mrSges(5,2) - t102 * mrSges(5,3);
t6 = m(4) * t133 + t146 * t92 + (-mrSges(4,1) - mrSges(5,1)) * t90 + t148 * t66 + (-t74 - t78) * t102 + t128;
t97 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t142;
t98 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t141;
t155 = m(3) * t84 - t96 * mrSges(3,1) + t95 * mrSges(3,2) + (t115 * t97 - t118 * t98) * qJD(1) + t114 * t5 + t151 * t6;
t143 = -t118 * g(3) - t115 * t85;
t51 = -qJDD(2) * pkin(2) - t120 * pkin(8) + t94 * t142 - t143;
t127 = t65 * pkin(3) + t51 + (-t144 - t66) * qJ(4);
t150 = pkin(3) * t102;
t122 = -t65 * pkin(4) - t89 * qJ(5) + qJDD(5) - t127 + ((2 * qJD(4)) + t150 + t75) * t92;
t134 = m(7) * (-t40 * pkin(5) - t67 * pkin(9) + t69 * t55 + t122) + t28 * mrSges(7,2) - t27 * mrSges(7,1) + t45 * t37 - t44 * t36;
t126 = m(6) * t122 - t40 * mrSges(6,1) + t41 * mrSges(6,2) - t68 * t53 + t69 * t54 + t134;
t125 = m(5) * ((t153 - t150) * t92 + t127) + t65 * mrSges(5,1) + t91 * t78 - t126;
t154 = m(4) * t51 + t65 * mrSges(4,1) + (t76 - t77) * t92 + (mrSges(4,2) - mrSges(5,3)) * t66 + t91 * t74 + t125;
t93 = (-mrSges(3,1) * t118 + mrSges(3,2) * t115) * qJD(1);
t11 = m(3) * t143 + qJDD(2) * mrSges(3,1) - t95 * mrSges(3,3) + qJD(2) * t98 - t93 * t142 - t154;
t4 = m(3) * t138 - qJDD(2) * mrSges(3,2) + t96 * mrSges(3,3) - qJD(2) * t97 - t114 * t6 + t93 * t141 + t151 * t5;
t152 = t118 * t11 + t115 * t4;
t2 = m(2) * t136 + qJDD(1) * mrSges(2,1) - t121 * mrSges(2,2) - t155;
t1 = m(2) * t132 - t121 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t115 * t11 + t118 * t4;
t3 = [-m(1) * g(1) + t119 * t1 - t116 * t2, t1, t4, t5, -t65 * mrSges(5,2) - t91 * t72 + t130, t9, t13; -m(1) * g(2) + t116 * t1 + t119 * t2, t2, t11, t6, -t66 * mrSges(5,3) - t92 * t77 + t125, t8, t12; (-m(1) - m(2)) * g(3) + t152, -m(2) * g(3) + t152, t155, t154, t90 * mrSges(5,1) + t66 * mrSges(5,2) + t102 * t78 + t92 * t72 - t128, t126, t134;];
f_new  = t3;
