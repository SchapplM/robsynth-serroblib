% Calculate vector of cutting forces with Newton-Euler
% S6RRRPPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-05-07 07:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPPR10_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR10_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR10_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR10_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR10_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR10_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR10_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 06:58:26
% EndTime: 2019-05-07 06:58:37
% DurationCPUTime: 4.99s
% Computational Cost: add. (66378->211), mult. (145044->281), div. (0->0), fcn. (112288->12), ass. (0->109)
t107 = sin(pkin(6));
t112 = sin(qJ(2));
t115 = cos(qJ(2));
t136 = qJD(1) * qJD(2);
t94 = (-qJDD(1) * t115 + t112 * t136) * t107;
t111 = sin(qJ(3));
t153 = cos(qJ(3));
t109 = cos(pkin(6));
t103 = t109 * qJD(1) + qJD(2);
t101 = t103 ^ 2;
t102 = t109 * qJDD(1) + qJDD(2);
t138 = qJD(1) * t115;
t117 = qJD(1) ^ 2;
t113 = sin(qJ(1));
t116 = cos(qJ(1));
t132 = t113 * g(1) - t116 * g(2);
t152 = pkin(8) * t107;
t89 = qJDD(1) * pkin(1) + t117 * t152 + t132;
t143 = t109 * t89;
t129 = -t116 * g(1) - t113 * g(2);
t90 = -t117 * pkin(1) + qJDD(1) * t152 + t129;
t144 = t112 * t143 + t115 * t90;
t139 = qJD(1) * t107;
t92 = (-pkin(2) * t115 - pkin(9) * t112) * t139;
t43 = -t101 * pkin(2) + t102 * pkin(9) + (-g(3) * t112 + t92 * t138) * t107 + t144;
t151 = t109 * g(3);
t93 = (qJDD(1) * t112 + t115 * t136) * t107;
t44 = t94 * pkin(2) - t93 * pkin(9) - t151 + (-t89 + (pkin(2) * t112 - pkin(9) * t115) * t103 * qJD(1)) * t107;
t147 = t111 * t44 + t153 * t43;
t134 = t112 * t139;
t81 = -t153 * t103 + t111 * t134;
t82 = t111 * t103 + t153 * t134;
t61 = t81 * pkin(3) - t82 * qJ(4);
t86 = qJDD(3) + t94;
t133 = t107 * t138;
t99 = -qJD(3) + t133;
t98 = t99 ^ 2;
t155 = t98 * pkin(3) - t86 * qJ(4) + 0.2e1 * qJD(4) * t99 + t81 * t61 - t147;
t106 = sin(pkin(11));
t108 = cos(pkin(11));
t110 = sin(qJ(6));
t114 = cos(qJ(6));
t150 = t81 * t99;
t130 = -t111 * t43 + t153 * t44;
t28 = -t86 * pkin(3) - t98 * qJ(4) + t82 * t61 + qJDD(4) - t130;
t60 = -t81 * qJD(3) + t111 * t102 + t153 * t93;
t22 = (t81 * t82 - t86) * qJ(5) + (t60 - t150) * pkin(4) + t28;
t140 = t107 * t115;
t128 = -g(3) * t140 - t112 * t90 + t115 * t143;
t42 = -t102 * pkin(2) - t101 * pkin(9) + t92 * t134 - t128;
t118 = (-t60 - t150) * qJ(4) + t42 + (-t99 * pkin(3) - 0.2e1 * qJD(4)) * t82;
t59 = t82 * qJD(3) - t153 * t102 + t111 * t93;
t72 = t82 * pkin(4) + t99 * qJ(5);
t80 = t81 ^ 2;
t26 = -t80 * pkin(4) - t82 * t72 + (pkin(3) + qJ(5)) * t59 + t118;
t69 = t106 * t81 - t108 * t99;
t131 = -0.2e1 * qJD(5) * t69 - t106 * t26 + t108 * t22;
t49 = t106 * t59 + t108 * t86;
t68 = t106 * t99 + t108 * t81;
t16 = (t68 * t82 - t49) * pkin(10) + (t68 * t69 + t60) * pkin(5) + t131;
t135 = 0.2e1 * qJD(5) * t68 + t106 * t22 + t108 * t26;
t48 = -t106 * t86 + t108 * t59;
t55 = t82 * pkin(5) - t69 * pkin(10);
t67 = t68 ^ 2;
t17 = -t67 * pkin(5) + t48 * pkin(10) - t82 * t55 + t135;
t45 = -t110 * t69 + t114 * t68;
t33 = t45 * qJD(6) + t110 * t48 + t114 * t49;
t46 = t110 * t68 + t114 * t69;
t35 = -t45 * mrSges(7,1) + t46 * mrSges(7,2);
t79 = qJD(6) + t82;
t36 = -t79 * mrSges(7,2) + t45 * mrSges(7,3);
t57 = qJDD(6) + t60;
t14 = m(7) * (-t110 * t17 + t114 * t16) - t33 * mrSges(7,3) + t57 * mrSges(7,1) - t46 * t35 + t79 * t36;
t32 = -t46 * qJD(6) - t110 * t49 + t114 * t48;
t37 = t79 * mrSges(7,1) - t46 * mrSges(7,3);
t15 = m(7) * (t110 * t16 + t114 * t17) + t32 * mrSges(7,3) - t57 * mrSges(7,2) + t45 * t35 - t79 * t37;
t50 = -t68 * mrSges(6,1) + t69 * mrSges(6,2);
t53 = -t82 * mrSges(6,2) + t68 * mrSges(6,3);
t11 = m(6) * t131 + t60 * mrSges(6,1) - t49 * mrSges(6,3) + t110 * t15 + t114 * t14 - t69 * t50 + t82 * t53;
t54 = t82 * mrSges(6,1) - t69 * mrSges(6,3);
t12 = m(6) * t135 - t60 * mrSges(6,2) + t48 * mrSges(6,3) - t110 * t14 + t114 * t15 + t68 * t50 - t82 * t54;
t74 = t82 * mrSges(5,1) - t99 * mrSges(5,2);
t126 = t106 * t11 - t108 * t12 - m(5) * (t59 * pkin(3) + t118) + t60 * mrSges(5,3) + t82 * t74;
t73 = t81 * mrSges(5,1) + t99 * mrSges(5,3);
t145 = t99 * mrSges(4,2) - t81 * mrSges(4,3) - t73;
t149 = mrSges(4,1) - mrSges(5,2);
t71 = -t99 * mrSges(4,1) - t82 * mrSges(4,3);
t154 = m(4) * t42 + t60 * mrSges(4,2) + t145 * t81 + t149 * t59 + t82 * t71 - t126;
t148 = -mrSges(4,3) - mrSges(5,1);
t63 = -t81 * mrSges(5,2) - t82 * mrSges(5,3);
t146 = -t81 * mrSges(4,1) - t82 * mrSges(4,2) - t63;
t141 = t107 * t112;
t121 = -t59 * pkin(4) - t80 * qJ(5) - t99 * t72 + qJDD(5) - t155;
t123 = -t32 * mrSges(7,1) - t45 * t36 + m(7) * (-t48 * pkin(5) - t67 * pkin(10) + t69 * t55 + t121) + t33 * mrSges(7,2) + t46 * t37;
t120 = m(6) * t121 - t48 * mrSges(6,1) + t49 * mrSges(6,2) - t68 * t53 + t69 * t54 + t123;
t119 = -m(5) * t155 + t120;
t13 = (t71 - t74) * t99 + (-mrSges(4,2) + mrSges(5,3)) * t86 + t146 * t81 + t148 * t59 + m(4) * t147 + t119;
t87 = t103 * mrSges(3,1) - mrSges(3,3) * t134;
t124 = -m(5) * t28 - t106 * t12 - t108 * t11;
t9 = m(4) * t130 - t145 * t99 + t146 * t82 + t148 * t60 + t149 * t86 + t124;
t91 = (-mrSges(3,1) * t115 + mrSges(3,2) * t112) * t139;
t4 = m(3) * (-g(3) * t141 + t144) - t94 * mrSges(3,3) - t102 * mrSges(3,2) + t91 * t133 - t103 * t87 + t153 * t13 - t111 * t9;
t88 = -t103 * mrSges(3,2) + mrSges(3,3) * t133;
t6 = m(3) * (-t107 * t89 - t151) + t93 * mrSges(3,2) + t94 * mrSges(3,1) + t111 * t13 + t153 * t9 + (t112 * t87 - t115 * t88) * t139;
t8 = m(3) * t128 + t102 * mrSges(3,1) - t93 * mrSges(3,3) + t103 * t88 - t91 * t134 - t154;
t137 = t109 * t6 + t8 * t140 + t4 * t141;
t2 = m(2) * t129 - t117 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t112 * t8 + t115 * t4;
t1 = m(2) * t132 + qJDD(1) * mrSges(2,1) - t117 * mrSges(2,2) - t107 * t6 + (t112 * t4 + t115 * t8) * t109;
t3 = [-m(1) * g(1) - t113 * t1 + t116 * t2, t2, t4, t13, -t59 * mrSges(5,2) - t81 * t73 - t126, t12, t15; -m(1) * g(2) + t116 * t1 + t113 * t2, t1, t8, t9, t59 * mrSges(5,1) - t86 * mrSges(5,3) + t81 * t63 + t99 * t74 - t119, t11, t14; (-m(1) - m(2)) * g(3) + t137, -m(2) * g(3) + t137, t6, t154, t60 * mrSges(5,1) + t86 * mrSges(5,2) + t82 * t63 - t99 * t73 - t124, t120, t123;];
f_new  = t3;
