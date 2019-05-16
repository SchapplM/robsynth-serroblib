% Calculate vector of cutting forces with Newton-Euler
% S6RRRPPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-05-07 06:43
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPPR9_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR9_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR9_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR9_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR9_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR9_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR9_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 06:31:53
% EndTime: 2019-05-07 06:32:05
% DurationCPUTime: 5.01s
% Computational Cost: add. (66448->214), mult. (145476->282), div. (0->0), fcn. (113754->12), ass. (0->110)
t108 = sin(pkin(6));
t112 = sin(qJ(2));
t115 = cos(qJ(2));
t135 = qJD(1) * qJD(2);
t96 = (-qJDD(1) * t115 + t112 * t135) * t108;
t137 = qJD(1) * t115;
t132 = t108 * t137;
t101 = qJD(3) - t132;
t100 = t101 ^ 2;
t111 = sin(qJ(3));
t153 = cos(qJ(3));
t109 = cos(pkin(6));
t104 = t109 * qJD(1) + qJD(2);
t102 = t104 ^ 2;
t103 = t109 * qJDD(1) + qJDD(2);
t117 = qJD(1) ^ 2;
t113 = sin(qJ(1));
t116 = cos(qJ(1));
t130 = t113 * g(1) - t116 * g(2);
t152 = pkin(8) * t108;
t91 = qJDD(1) * pkin(1) + t117 * t152 + t130;
t142 = t109 * t91;
t128 = -t116 * g(1) - t113 * g(2);
t92 = -t117 * pkin(1) + qJDD(1) * t152 + t128;
t143 = t112 * t142 + t115 * t92;
t138 = qJD(1) * t108;
t94 = (-pkin(2) * t115 - pkin(9) * t112) * t138;
t47 = -t102 * pkin(2) + t103 * pkin(9) + (-g(3) * t112 + t94 * t137) * t108 + t143;
t151 = t109 * g(3);
t95 = (qJDD(1) * t112 + t115 * t135) * t108;
t48 = t96 * pkin(2) - t95 * pkin(9) - t151 + (-t91 + (pkin(2) * t112 - pkin(9) * t115) * t104 * qJD(1)) * t108;
t147 = -t111 * t47 + t153 * t48;
t133 = t112 * t138;
t84 = -t153 * t104 + t111 * t133;
t85 = t111 * t104 + t153 * t133;
t69 = t84 * pkin(3) - t85 * qJ(4);
t88 = qJDD(3) + t96;
t120 = t88 * pkin(3) + t100 * qJ(4) - t85 * t69 - qJDD(4) + t147;
t107 = sin(pkin(11));
t141 = cos(pkin(11));
t74 = -t141 * t101 + t107 * t85;
t150 = t74 * t84;
t68 = -t84 * qJD(3) + t111 * t103 + t153 * t95;
t54 = t107 * t88 + t141 * t68;
t157 = (-t54 + t150) * qJ(5) - t120;
t154 = 2 * qJD(5);
t110 = sin(qJ(6));
t114 = cos(qJ(6));
t75 = t107 * t101 + t141 * t85;
t51 = t110 * t74 + t114 * t75;
t53 = t107 * t68 - t141 * t88;
t32 = -t51 * qJD(6) - t110 * t54 + t114 * t53;
t50 = -t110 * t75 + t114 * t74;
t33 = t50 * qJD(6) + t110 * t53 + t114 * t54;
t81 = qJD(6) - t84;
t39 = -t81 * mrSges(7,2) + t50 * mrSges(7,3);
t40 = t81 * mrSges(7,1) - t51 * mrSges(7,3);
t63 = -t84 * pkin(5) - t75 * pkin(10);
t73 = t74 ^ 2;
t129 = m(7) * (-t73 * pkin(10) + (-pkin(4) - pkin(5)) * t53 + (-pkin(4) * t84 + t154 + t63) * t75 - t157) + t33 * mrSges(7,2) - t32 * mrSges(7,1) + t51 * t40 - t50 * t39;
t59 = -t74 * mrSges(6,2) + t84 * mrSges(6,3);
t123 = m(6) * (-0.2e1 * qJD(5) * t75 + (t75 * t84 + t53) * pkin(4) + t157) + t74 * t59 + t53 * mrSges(6,1) - t129;
t60 = -t84 * mrSges(5,2) - t74 * mrSges(5,3);
t61 = t84 * mrSges(5,1) - t75 * mrSges(5,3);
t62 = -t84 * mrSges(6,1) + t75 * mrSges(6,2);
t156 = -m(5) * t120 + t53 * mrSges(5,1) + (t61 - t62) * t75 + (mrSges(5,2) - mrSges(6,3)) * t54 + t74 * t60 + t123;
t155 = -2 * qJD(4);
t148 = -mrSges(5,3) - mrSges(6,2);
t146 = t111 * t48 + t153 * t47;
t56 = t74 * mrSges(6,1) - t75 * mrSges(6,3);
t145 = -t74 * mrSges(5,1) - t75 * mrSges(5,2) - t56;
t140 = t108 * t112;
t139 = t108 * t115;
t70 = t84 * mrSges(4,1) + t85 * mrSges(4,2);
t76 = -t101 * mrSges(4,2) - t84 * mrSges(4,3);
t12 = m(4) * t147 + t88 * mrSges(4,1) - t68 * mrSges(4,3) + t101 * t76 - t85 * t70 - t156;
t89 = t104 * mrSges(3,1) - mrSges(3,3) * t133;
t27 = -t100 * pkin(3) + t88 * qJ(4) - t84 * t69 + t146;
t126 = -g(3) * t139 - t112 * t92 + t115 * t142;
t46 = -t103 * pkin(2) - t102 * pkin(9) + t94 * t133 - t126;
t67 = t85 * qJD(3) - t153 * t103 + t111 * t95;
t29 = (t101 * t84 - t68) * qJ(4) + (t101 * t85 + t67) * pkin(3) + t46;
t134 = t107 * t29 + t141 * t27 + t74 * t155;
t55 = t74 * pkin(4) - t75 * qJ(5);
t83 = t84 ^ 2;
t121 = -t83 * pkin(4) + t67 * qJ(5) + t84 * t154 - t74 * t55 + t134;
t127 = -t107 * t27 + t141 * t29;
t21 = -t67 * pkin(4) - t83 * qJ(5) + qJDD(5) - t127 + ((2 * qJD(4)) + t55) * t75;
t16 = (-t54 - t150) * pkin(10) + (t74 * t75 - t67) * pkin(5) + t21;
t17 = -t73 * pkin(5) + t53 * pkin(10) + t84 * t63 + t121;
t36 = -t50 * mrSges(7,1) + t51 * mrSges(7,2);
t66 = qJDD(6) - t67;
t14 = m(7) * (-t110 * t17 + t114 * t16) - t33 * mrSges(7,3) + t66 * mrSges(7,1) - t51 * t36 + t81 * t39;
t15 = m(7) * (t110 * t16 + t114 * t17) + t32 * mrSges(7,3) - t66 * mrSges(7,2) + t50 * t36 - t81 * t40;
t125 = m(6) * t121 + t67 * mrSges(6,3) - t110 * t14 + t114 * t15 + t84 * t62;
t10 = m(5) * t134 - t67 * mrSges(5,2) + t145 * t74 + t148 * t53 - t84 * t61 + t125;
t122 = -m(6) * t21 - t110 * t15 - t114 * t14;
t11 = m(5) * t127 + (t60 + t59) * t84 + (m(5) * t155 + t145) * t75 + (mrSges(5,1) + mrSges(6,1)) * t67 + t148 * t54 + t122;
t77 = t101 * mrSges(4,1) - t85 * mrSges(4,3);
t9 = m(4) * t146 - t88 * mrSges(4,2) - t67 * mrSges(4,3) + t141 * t10 - t101 * t77 - t107 * t11 - t84 * t70;
t93 = (-mrSges(3,1) * t115 + mrSges(3,2) * t112) * t138;
t4 = m(3) * (-g(3) * t140 + t143) - t96 * mrSges(3,3) - t103 * mrSges(3,2) + t93 * t132 - t104 * t89 + t153 * t9 - t111 * t12;
t90 = -t104 * mrSges(3,2) + mrSges(3,3) * t132;
t6 = m(3) * (-t108 * t91 - t151) + t95 * mrSges(3,2) + t96 * mrSges(3,1) + t111 * t9 + t153 * t12 + (t112 * t89 - t115 * t90) * t138;
t118 = m(4) * t46 + t67 * mrSges(4,1) + t68 * mrSges(4,2) + t107 * t10 + t141 * t11 + t84 * t76 + t85 * t77;
t8 = m(3) * t126 + t103 * mrSges(3,1) - t95 * mrSges(3,3) + t104 * t90 - t93 * t133 - t118;
t136 = t109 * t6 + t8 * t139 + t4 * t140;
t2 = m(2) * t128 - t117 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t112 * t8 + t115 * t4;
t1 = m(2) * t130 + qJDD(1) * mrSges(2,1) - t117 * mrSges(2,2) - t108 * t6 + (t112 * t4 + t115 * t8) * t109;
t3 = [-m(1) * g(1) - t113 * t1 + t116 * t2, t2, t4, t9, t10, -t53 * mrSges(6,2) - t74 * t56 + t125, t15; -m(1) * g(2) + t116 * t1 + t113 * t2, t1, t8, t12, t11, -t54 * mrSges(6,3) - t75 * t62 + t123, t14; (-m(1) - m(2)) * g(3) + t136, -m(2) * g(3) + t136, t6, t118, t156, -t67 * mrSges(6,1) + t54 * mrSges(6,2) + t75 * t56 - t84 * t59 - t122, t129;];
f_new  = t3;
