% Calculate vector of cutting forces with Newton-Euler
% S6RRRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-07 20:14
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRPR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 20:06:48
% EndTime: 2019-05-07 20:07:05
% DurationCPUTime: 6.62s
% Computational Cost: add. (126297->208), mult. (255363->275), div. (0->0), fcn. (188286->12), ass. (0->106)
t111 = sin(qJ(2));
t116 = cos(qJ(2));
t118 = qJD(1) ^ 2;
t106 = sin(pkin(11));
t107 = cos(pkin(11));
t108 = sin(qJ(6));
t113 = cos(qJ(6));
t109 = sin(qJ(4));
t114 = cos(qJ(4));
t104 = qJD(2) + qJD(3);
t105 = t116 ^ 2;
t112 = sin(qJ(1));
t117 = cos(qJ(1));
t132 = t112 * g(1) - t117 * g(2);
t126 = -qJDD(1) * pkin(1) - t132;
t137 = qJD(1) * t111;
t135 = qJD(1) * qJD(2);
t96 = qJDD(1) * t116 - t111 * t135;
t99 = qJD(2) * pkin(2) - pkin(8) * t137;
t123 = -t96 * pkin(2) + t99 * t137 + (-pkin(8) * t105 - pkin(7)) * t118 + t126;
t110 = sin(qJ(3));
t115 = cos(qJ(3));
t89 = (t110 * t116 + t111 * t115) * qJD(1);
t95 = qJDD(1) * t111 + t116 * t135;
t70 = -t89 * qJD(3) - t110 * t95 + t115 * t96;
t136 = qJD(1) * t116;
t88 = -t110 * t137 + t115 * t136;
t71 = qJD(3) * t88 + t110 * t96 + t115 * t95;
t35 = (-t104 * t88 - t71) * pkin(9) + (t104 * t89 - t70) * pkin(3) + t123;
t102 = t104 ^ 2;
t103 = qJDD(2) + qJDD(3);
t128 = -g(1) * t117 - g(2) * t112;
t91 = -pkin(1) * t118 + qJDD(1) * pkin(7) + t128;
t138 = t111 * t91;
t141 = pkin(2) * t118;
t63 = qJDD(2) * pkin(2) - t95 * pkin(8) - t138 + (pkin(8) * t135 + t111 * t141 - g(3)) * t116;
t133 = -g(3) * t111 + t116 * t91;
t64 = pkin(8) * t96 - qJD(2) * t99 - t105 * t141 + t133;
t139 = t110 * t63 + t115 * t64;
t78 = -pkin(3) * t88 - pkin(9) * t89;
t41 = -pkin(3) * t102 + pkin(9) * t103 + t78 * t88 + t139;
t131 = -t109 * t41 + t114 * t35;
t80 = t104 * t114 - t109 * t89;
t49 = qJD(4) * t80 + t103 * t109 + t114 * t71;
t68 = qJDD(4) - t70;
t81 = t104 * t109 + t114 * t89;
t87 = qJD(4) - t88;
t23 = (t80 * t87 - t49) * qJ(5) + (t80 * t81 + t68) * pkin(4) + t131;
t140 = t109 * t35 + t114 * t41;
t48 = -qJD(4) * t81 + t103 * t114 - t109 * t71;
t74 = pkin(4) * t87 - qJ(5) * t81;
t79 = t80 ^ 2;
t25 = -pkin(4) * t79 + qJ(5) * t48 - t74 * t87 + t140;
t58 = t106 * t80 + t107 * t81;
t129 = -0.2e1 * qJD(5) * t58 - t106 * t25 + t107 * t23;
t38 = t106 * t48 + t107 * t49;
t57 = -t106 * t81 + t107 * t80;
t17 = (t57 * t87 - t38) * pkin(10) + (t57 * t58 + t68) * pkin(5) + t129;
t134 = 0.2e1 * qJD(5) * t57 + t106 * t23 + t107 * t25;
t37 = -t106 * t49 + t107 * t48;
t52 = pkin(5) * t87 - pkin(10) * t58;
t56 = t57 ^ 2;
t18 = -pkin(5) * t56 + pkin(10) * t37 - t52 * t87 + t134;
t45 = -t108 * t58 + t113 * t57;
t28 = qJD(6) * t45 + t108 * t37 + t113 * t38;
t46 = t108 * t57 + t113 * t58;
t32 = -mrSges(7,1) * t45 + mrSges(7,2) * t46;
t85 = qJD(6) + t87;
t42 = -mrSges(7,2) * t85 + mrSges(7,3) * t45;
t65 = qJDD(6) + t68;
t14 = m(7) * (-t108 * t18 + t113 * t17) - t28 * mrSges(7,3) + t65 * mrSges(7,1) - t46 * t32 + t85 * t42;
t27 = -qJD(6) * t46 - t108 * t38 + t113 * t37;
t43 = mrSges(7,1) * t85 - mrSges(7,3) * t46;
t15 = m(7) * (t108 * t17 + t113 * t18) + t27 * mrSges(7,3) - t65 * mrSges(7,2) + t45 * t32 - t85 * t43;
t47 = -mrSges(6,1) * t57 + mrSges(6,2) * t58;
t50 = -mrSges(6,2) * t87 + mrSges(6,3) * t57;
t12 = m(6) * t129 + t68 * mrSges(6,1) - t38 * mrSges(6,3) + t108 * t15 + t113 * t14 - t58 * t47 + t87 * t50;
t51 = mrSges(6,1) * t87 - mrSges(6,3) * t58;
t13 = m(6) * t134 - t68 * mrSges(6,2) + t37 * mrSges(6,3) - t108 * t14 + t113 * t15 + t57 * t47 - t87 * t51;
t62 = -mrSges(5,1) * t80 + mrSges(5,2) * t81;
t73 = -mrSges(5,2) * t87 + mrSges(5,3) * t80;
t10 = m(5) * t131 + t68 * mrSges(5,1) - t49 * mrSges(5,3) + t106 * t13 + t107 * t12 - t81 * t62 + t87 * t73;
t75 = mrSges(5,1) * t87 - mrSges(5,3) * t81;
t11 = m(5) * t140 - t68 * mrSges(5,2) + t48 * mrSges(5,3) - t106 * t12 + t107 * t13 + t80 * t62 - t87 * t75;
t82 = -mrSges(4,2) * t104 + mrSges(4,3) * t88;
t83 = mrSges(4,1) * t104 - mrSges(4,3) * t89;
t124 = -m(4) * t123 + t70 * mrSges(4,1) - t71 * mrSges(4,2) - t114 * t10 - t109 * t11 + t88 * t82 - t89 * t83;
t97 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t137;
t98 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t136;
t143 = (t111 * t97 - t116 * t98) * qJD(1) + m(3) * (-t118 * pkin(7) + t126) - t96 * mrSges(3,1) + t95 * mrSges(3,2) - t124;
t130 = -t110 * t64 + t115 * t63;
t40 = -pkin(3) * t103 - pkin(9) * t102 + t89 * t78 - t130;
t120 = -pkin(4) * t48 - qJ(5) * t79 + t81 * t74 + qJDD(5) + t40;
t125 = t27 * mrSges(7,1) + t45 * t42 - m(7) * (-pkin(5) * t37 - pkin(10) * t56 + t52 * t58 + t120) - t28 * mrSges(7,2) - t46 * t43;
t122 = -m(6) * t120 + t37 * mrSges(6,1) - t38 * mrSges(6,2) + t57 * t50 - t58 * t51 + t125;
t119 = m(5) * t40 - t48 * mrSges(5,1) + t49 * mrSges(5,2) - t80 * t73 + t81 * t75 - t122;
t77 = -mrSges(4,1) * t88 + mrSges(4,2) * t89;
t16 = m(4) * t130 + t103 * mrSges(4,1) - t71 * mrSges(4,3) + t104 * t82 - t89 * t77 - t119;
t7 = m(4) * t139 - t103 * mrSges(4,2) + t70 * mrSges(4,3) - t109 * t10 - t104 * t83 + t114 * t11 + t88 * t77;
t94 = (-mrSges(3,1) * t116 + mrSges(3,2) * t111) * qJD(1);
t4 = m(3) * (-t116 * g(3) - t138) - t95 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t94 * t137 + qJD(2) * t98 + t110 * t7 + t115 * t16;
t5 = m(3) * t133 - qJDD(2) * mrSges(3,2) + t96 * mrSges(3,3) - qJD(2) * t97 - t110 * t16 + t115 * t7 + t94 * t136;
t142 = t111 * t5 + t116 * t4;
t6 = m(2) * t132 + qJDD(1) * mrSges(2,1) - t118 * mrSges(2,2) - t143;
t1 = m(2) * t128 - t118 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t111 * t4 + t116 * t5;
t2 = [-m(1) * g(1) + t1 * t117 - t112 * t6, t1, t5, t7, t11, t13, t15; -m(1) * g(2) + t1 * t112 + t117 * t6, t6, t4, t16, t10, t12, t14; (-m(1) - m(2)) * g(3) + t142, -m(2) * g(3) + t142, t143, -t124, t119, -t122, -t125;];
f_new  = t2;
