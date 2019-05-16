% Calculate vector of cutting forces with Newton-Euler
% S6RRPRRP14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-05-06 19:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRRP14_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP14_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP14_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP14_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP14_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP14_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP14_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 19:17:22
% EndTime: 2019-05-06 19:17:29
% DurationCPUTime: 2.65s
% Computational Cost: add. (29699->213), mult. (66163->269), div. (0->0), fcn. (47508->10), ass. (0->105)
t158 = -2 * qJD(3);
t105 = sin(qJ(2));
t101 = sin(pkin(6));
t134 = qJD(1) * t101;
t125 = t105 * t134;
t102 = cos(pkin(6));
t96 = t102 * qJD(1) + qJD(2);
t157 = (pkin(2) * t96 + t158) * t125;
t108 = cos(qJ(2));
t136 = t101 * t105;
t110 = qJD(1) ^ 2;
t106 = sin(qJ(1));
t109 = cos(qJ(1));
t124 = t106 * g(1) - t109 * g(2);
t78 = t110 * t101 * pkin(8) + qJDD(1) * pkin(1) + t124;
t139 = t102 * t78;
t122 = -t109 * g(1) - t106 * g(2);
t131 = qJDD(1) * t101;
t79 = -t110 * pkin(1) + pkin(8) * t131 + t122;
t120 = -g(3) * t136 + t105 * t139 + t108 * t79;
t133 = qJD(1) * t108;
t126 = t101 * t133;
t80 = (-pkin(2) * t108 - qJ(3) * t105) * t134;
t94 = t96 ^ 2;
t95 = t102 * qJDD(1) + qJDD(2);
t156 = t94 * pkin(2) - t95 * qJ(3) - t80 * t126 + t96 * t158 - t120;
t104 = sin(qJ(4));
t107 = cos(qJ(4));
t138 = t110 * t101 ^ 2;
t128 = t108 ^ 2 * t138;
t151 = t102 * g(3);
t154 = -pkin(2) - pkin(9);
t83 = pkin(3) * t125 - t96 * pkin(9);
t84 = (qJD(2) * t133 + qJDD(1) * t105) * t101;
t85 = -qJD(2) * t125 + t108 * t131;
t29 = -pkin(3) * t128 - t151 - t84 * qJ(3) + t154 * t85 + (-t78 + (-qJ(3) * t108 * t96 - t105 * t83) * qJD(1)) * t101 + t157;
t135 = t101 * t108;
t142 = g(3) * t135 + t105 * t79;
t117 = -t94 * qJ(3) + t80 * t125 + qJDD(3) + t142;
t31 = t84 * pkin(3) + t154 * t95 + (-pkin(3) * t96 * t134 - pkin(9) * t105 * t138 - t139) * t108 + t117;
t123 = -t104 * t29 + t107 * t31;
t68 = -t104 * t96 - t107 * t126;
t69 = -t104 * t126 + t107 * t96;
t55 = -t68 * pkin(4) - t69 * pkin(10);
t73 = qJDD(4) + t84;
t90 = qJD(4) + t125;
t88 = t90 ^ 2;
t21 = -t73 * pkin(4) - t88 * pkin(10) + t69 * t55 - t123;
t103 = sin(qJ(5));
t152 = cos(qJ(5));
t53 = t68 * qJD(4) - t104 * t85 + t107 * t95;
t57 = t103 * t90 + t152 * t69;
t33 = t57 * qJD(5) + t103 * t53 - t152 * t73;
t56 = t103 * t69 - t152 * t90;
t34 = -t56 * qJD(5) + t103 * t73 + t152 * t53;
t66 = qJD(5) - t68;
t44 = -t56 * mrSges(7,2) + t66 * mrSges(7,3);
t129 = m(7) * (-0.2e1 * qJD(6) * t57 + (t56 * t66 - t34) * qJ(6) + (t57 * t66 + t33) * pkin(5) + t21) + t33 * mrSges(7,1) + t56 * t44;
t45 = -t66 * mrSges(6,2) - t56 * mrSges(6,3);
t46 = t66 * mrSges(6,1) - t57 * mrSges(6,3);
t47 = -t66 * mrSges(7,1) + t57 * mrSges(7,2);
t155 = m(6) * t21 + t33 * mrSges(6,1) + (t46 - t47) * t57 + (mrSges(6,2) - mrSges(7,3)) * t34 + t56 * t45 + t129;
t145 = t104 * t31 + t107 * t29;
t22 = -t88 * pkin(4) + t73 * pkin(10) + t68 * t55 + t145;
t111 = t85 * pkin(3) - pkin(9) * t128 + t96 * t83 - t156;
t52 = -t69 * qJD(4) - t104 * t95 - t107 * t85;
t24 = (-t68 * t90 - t53) * pkin(10) + (t69 * t90 - t52) * pkin(4) + t111;
t118 = -t103 * t22 + t152 * t24;
t40 = t56 * pkin(5) - t57 * qJ(6);
t50 = qJDD(5) - t52;
t65 = t66 ^ 2;
t153 = m(7) * (-t50 * pkin(5) - t65 * qJ(6) + t57 * t40 + qJDD(6) - t118);
t150 = mrSges(3,1) - mrSges(4,2);
t148 = mrSges(3,3) + mrSges(4,1);
t147 = -mrSges(6,3) - mrSges(7,2);
t146 = t103 * t24 + t152 * t22;
t41 = t56 * mrSges(7,1) - t57 * mrSges(7,3);
t144 = -t56 * mrSges(6,1) - t57 * mrSges(6,2) - t41;
t77 = mrSges(4,1) * t125 + t96 * mrSges(4,2);
t141 = t96 * mrSges(3,1) - mrSges(3,3) * t125 - t77;
t81 = (mrSges(4,2) * t108 - mrSges(4,3) * t105) * t134;
t140 = t81 + (-mrSges(3,1) * t108 + mrSges(3,2) * t105) * t134;
t130 = m(7) * (-t65 * pkin(5) + t50 * qJ(6) + 0.2e1 * qJD(6) * t66 - t56 * t40 + t146) + t66 * t47 + t50 * mrSges(7,3);
t13 = m(6) * t146 - t50 * mrSges(6,2) + t144 * t56 + t147 * t33 - t66 * t46 + t130;
t15 = m(6) * t118 - t153 + (t45 + t44) * t66 + t144 * t57 + (mrSges(6,1) + mrSges(7,1)) * t50 + t147 * t34;
t54 = -t68 * mrSges(5,1) + t69 * mrSges(5,2);
t59 = t90 * mrSges(5,1) - t69 * mrSges(5,3);
t10 = m(5) * t145 - t73 * mrSges(5,2) + t52 * mrSges(5,3) - t103 * t15 + t152 * t13 + t68 * t54 - t90 * t59;
t58 = -t90 * mrSges(5,2) + t68 * mrSges(5,3);
t11 = m(5) * t123 + t73 * mrSges(5,1) - t53 * mrSges(5,3) - t69 * t54 + t90 * t58 - t155;
t121 = -t101 * t78 - t151;
t76 = -mrSges(4,1) * t126 - t96 * mrSges(4,3);
t119 = t107 * t10 - t104 * t11 + m(4) * (-t85 * pkin(2) + (-t96 * t126 - t84) * qJ(3) + t121 + t157) + t76 * t126 - t84 * mrSges(4,3);
t75 = -t96 * mrSges(3,2) + mrSges(3,3) * t126;
t5 = m(3) * t121 + t84 * mrSges(3,2) - t150 * t85 + (t141 * t105 - t108 * t75) * t134 + t119;
t127 = t108 * t139;
t116 = -m(4) * (-t95 * pkin(2) + t117 - t127) - t104 * t10 - t107 * t11;
t6 = m(3) * (t127 - t142) + (t75 - t76) * t96 + t150 * t95 - t148 * t84 - t140 * t125 + t116;
t114 = m(5) * t111 - t52 * mrSges(5,1) + t53 * mrSges(5,2) + t103 * t13 + t152 * t15 - t68 * t58 + t69 * t59;
t112 = -m(4) * t156 + t114;
t8 = m(3) * t120 - t141 * t96 + (-mrSges(3,2) + mrSges(4,3)) * t95 + t148 * t85 + t140 * t126 + t112;
t132 = t102 * t5 + t6 * t135 + t8 * t136;
t2 = m(2) * t122 - t110 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t105 * t6 + t108 * t8;
t1 = m(2) * t124 + qJDD(1) * mrSges(2,1) - t110 * mrSges(2,2) - t101 * t5 + (t105 * t8 + t108 * t6) * t102;
t3 = [-m(1) * g(1) - t106 * t1 + t109 * t2, t2, t8, t85 * mrSges(4,2) - t77 * t125 + t119, t10, t13, -t33 * mrSges(7,2) - t56 * t41 + t130; -m(1) * g(2) + t109 * t1 + t106 * t2, t1, t6, -t85 * mrSges(4,1) - t95 * mrSges(4,3) - t81 * t126 - t96 * t77 - t112, t11, t15, -t34 * mrSges(7,3) - t57 * t47 + t129; (-m(1) - m(2)) * g(3) + t132, -m(2) * g(3) + t132, t5, t84 * mrSges(4,1) + t95 * mrSges(4,2) + t81 * t125 + t96 * t76 - t116, t114, t155, -t50 * mrSges(7,1) + t34 * mrSges(7,2) + t57 * t41 - t66 * t44 + t153;];
f_new  = t3;
