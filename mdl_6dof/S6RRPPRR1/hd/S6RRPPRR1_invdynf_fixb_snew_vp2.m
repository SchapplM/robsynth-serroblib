% Calculate vector of cutting forces with Newton-Euler
% S6RRPPRR1
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
% Datum: 2019-05-06 09:35
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPPRR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:31:01
% EndTime: 2019-05-06 09:31:08
% DurationCPUTime: 2.69s
% Computational Cost: add. (31432->203), mult. (73313->259), div. (0->0), fcn. (51018->10), ass. (0->100)
t112 = sin(qJ(2));
t116 = cos(qJ(2));
t139 = qJD(1) * t112;
t119 = qJD(1) ^ 2;
t140 = t116 ^ 2 * t119;
t113 = sin(qJ(1));
t117 = cos(qJ(1));
t143 = t113 * g(1) - t117 * g(2);
t86 = -qJDD(1) * pkin(1) - t119 * pkin(7) - t143;
t137 = qJD(1) * qJD(2);
t92 = t116 * qJDD(1) - t112 * t137;
t93 = qJD(2) * pkin(2) - qJ(3) * t139;
t126 = -t92 * pkin(2) - qJ(3) * t140 + t93 * t139 + qJDD(3) + t86;
t108 = sin(pkin(10));
t138 = qJD(1) * t116;
t142 = cos(pkin(10));
t83 = t108 * t139 - t142 * t138;
t141 = qJD(2) * t83;
t91 = t112 * qJDD(1) + t116 * t137;
t69 = t108 * t91 - t142 * t92;
t70 = t108 * t92 + t142 * t91;
t125 = t69 * pkin(3) + t126 + (t141 - t70) * qJ(4);
t110 = sin(qJ(6));
t114 = cos(qJ(6));
t145 = pkin(3) * qJD(2);
t150 = 2 * qJD(4);
t84 = (t108 * t116 + t142 * t112) * qJD(1);
t76 = -qJD(2) * pkin(4) - t84 * pkin(8);
t82 = t83 ^ 2;
t122 = -t69 * pkin(4) - t82 * pkin(8) - t125 + (-t145 + t150 + t76) * t84;
t103 = -qJD(2) + qJD(5);
t101 = t103 ^ 2;
t102 = -qJDD(2) + qJDD(5);
t111 = sin(qJ(5));
t115 = cos(qJ(5));
t118 = qJD(2) ^ 2;
t134 = -t117 * g(1) - t113 * g(2);
t87 = -t119 * pkin(1) + qJDD(1) * pkin(7) + t134;
t144 = t112 * t87;
t46 = qJDD(2) * pkin(2) - t91 * qJ(3) - t144 + (pkin(2) * t112 * t119 + qJ(3) * t137 - g(3)) * t116;
t135 = -t112 * g(3) + t116 * t87;
t47 = -pkin(2) * t140 + t92 * qJ(3) - qJD(2) * t93 + t135;
t133 = -t108 * t47 + t142 * t46;
t62 = t83 * pkin(3) - t84 * qJ(4);
t26 = -qJDD(2) * pkin(3) - t118 * qJ(4) + qJDD(4) - t133 + ((2 * qJD(3)) + t62) * t84;
t20 = (-t70 - t141) * pkin(8) + (t83 * t84 - qJDD(2)) * pkin(4) + t26;
t151 = -2 * qJD(3);
t136 = t108 * t46 + t142 * t47 + t83 * t151;
t127 = -t118 * pkin(3) + qJDD(2) * qJ(4) + qJD(2) * t150 - t83 * t62 + t136;
t22 = -t82 * pkin(4) + t69 * pkin(8) + qJD(2) * t76 + t127;
t147 = t111 * t20 + t115 * t22;
t57 = -t111 * t84 + t115 * t83;
t58 = t111 * t83 + t115 * t84;
t42 = -t57 * pkin(5) - t58 * pkin(9);
t17 = -t101 * pkin(5) + t102 * pkin(9) + t57 * t42 + t147;
t34 = -t58 * qJD(5) - t111 * t70 + t115 * t69;
t35 = t57 * qJD(5) + t111 * t69 + t115 * t70;
t18 = (-t103 * t57 - t35) * pkin(9) + (t103 * t58 - t34) * pkin(5) + t122;
t52 = t114 * t103 - t110 * t58;
t28 = t52 * qJD(6) + t110 * t102 + t114 * t35;
t33 = qJDD(6) - t34;
t53 = t110 * t103 + t114 * t58;
t36 = -t52 * mrSges(7,1) + t53 * mrSges(7,2);
t56 = qJD(6) - t57;
t37 = -t56 * mrSges(7,2) + t52 * mrSges(7,3);
t14 = m(7) * (-t110 * t17 + t114 * t18) - t28 * mrSges(7,3) + t33 * mrSges(7,1) - t53 * t36 + t56 * t37;
t27 = -t53 * qJD(6) + t114 * t102 - t110 * t35;
t38 = t56 * mrSges(7,1) - t53 * mrSges(7,3);
t15 = m(7) * (t110 * t18 + t114 * t17) + t27 * mrSges(7,3) - t33 * mrSges(7,2) + t52 * t36 - t56 * t38;
t54 = -t103 * mrSges(6,2) + t57 * mrSges(6,3);
t55 = t103 * mrSges(6,1) - t58 * mrSges(6,3);
t129 = m(6) * t122 - t34 * mrSges(6,1) + t35 * mrSges(6,2) + t110 * t15 + t114 * t14 - t57 * t54 + t58 * t55;
t74 = -qJD(2) * mrSges(5,1) + t84 * mrSges(5,2);
t75 = -t83 * mrSges(5,2) + qJD(2) * mrSges(5,3);
t123 = t70 * mrSges(5,3) + t84 * t74 + t129 - m(5) * ((-(2 * qJD(4)) + t145) * t84 + t125) - t83 * t75 - t69 * mrSges(5,1);
t72 = -qJD(2) * mrSges(4,2) - t83 * mrSges(4,3);
t73 = qJD(2) * mrSges(4,1) - t84 * mrSges(4,3);
t121 = m(4) * t126 + t69 * mrSges(4,1) + t70 * mrSges(4,2) + t83 * t72 + t84 * t73 - t123;
t94 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t139;
t95 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t138;
t152 = m(3) * t86 - t92 * mrSges(3,1) + t91 * mrSges(3,2) + (t112 * t94 - t116 * t95) * qJD(1) + t121;
t41 = -t57 * mrSges(6,1) + t58 * mrSges(6,2);
t10 = m(6) * t147 - t102 * mrSges(6,2) + t34 * mrSges(6,3) - t103 * t55 - t110 * t14 + t114 * t15 + t57 * t41;
t132 = -t111 * t22 + t115 * t20;
t124 = m(7) * (-t102 * pkin(5) - t101 * pkin(9) + t58 * t42 - t132) - t27 * mrSges(7,1) + t28 * mrSges(7,2) - t52 * t37 + t53 * t38;
t11 = m(6) * t132 + t102 * mrSges(6,1) - t35 * mrSges(6,3) + t103 * t54 - t58 * t41 - t124;
t130 = m(5) * t127 + qJDD(2) * mrSges(5,3) + qJD(2) * t74 + t115 * t10 - t111 * t11;
t63 = t83 * mrSges(5,1) - t84 * mrSges(5,3);
t146 = -t83 * mrSges(4,1) - t84 * mrSges(4,2) - t63;
t148 = -mrSges(4,3) - mrSges(5,2);
t6 = m(4) * t136 - qJDD(2) * mrSges(4,2) - qJD(2) * t73 + t146 * t83 + t148 * t69 + t130;
t128 = -m(5) * t26 - t111 * t10 - t115 * t11;
t7 = m(4) * t133 + (m(4) * t151 + t146) * t84 + t148 * t70 + (mrSges(4,1) + mrSges(5,1)) * qJDD(2) + (t72 + t75) * qJD(2) + t128;
t90 = (-mrSges(3,1) * t116 + mrSges(3,2) * t112) * qJD(1);
t4 = m(3) * (-t116 * g(3) - t144) - t91 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t90 * t139 + qJD(2) * t95 + t108 * t6 + t142 * t7;
t5 = m(3) * t135 - qJDD(2) * mrSges(3,2) + t92 * mrSges(3,3) - qJD(2) * t94 - t108 * t7 + t90 * t138 + t142 * t6;
t149 = t112 * t5 + t116 * t4;
t8 = m(2) * t143 + qJDD(1) * mrSges(2,1) - t119 * mrSges(2,2) - t152;
t1 = m(2) * t134 - t119 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t112 * t4 + t116 * t5;
t2 = [-m(1) * g(1) + t117 * t1 - t113 * t8, t1, t5, t6, -t69 * mrSges(5,2) - t83 * t63 + t130, t10, t15; -m(1) * g(2) + t113 * t1 + t117 * t8, t8, t4, t7, -t123, t11, t14; (-m(1) - m(2)) * g(3) + t149, -m(2) * g(3) + t149, t152, t121, -qJDD(2) * mrSges(5,1) + t70 * mrSges(5,2) - qJD(2) * t75 + t84 * t63 - t128, t129, t124;];
f_new  = t2;
