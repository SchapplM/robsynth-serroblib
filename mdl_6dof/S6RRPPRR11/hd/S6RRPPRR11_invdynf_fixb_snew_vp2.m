% Calculate vector of cutting forces with Newton-Euler
% S6RRPPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-05-06 12:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPPRR11_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR11_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR11_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR11_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR11_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR11_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR11_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 11:58:57
% EndTime: 2019-05-06 11:59:07
% DurationCPUTime: 4.60s
% Computational Cost: add. (56035->215), mult. (132389->285), div. (0->0), fcn. (98080->12), ass. (0->112)
t158 = -2 * qJD(3);
t108 = cos(pkin(6));
t100 = t108 * qJD(1) + qJD(2);
t111 = sin(qJ(2));
t106 = sin(pkin(6));
t143 = qJD(1) * t106;
t134 = t111 * t143;
t157 = (pkin(2) * t100 + t158) * t134;
t115 = cos(qJ(2));
t145 = t106 * t111;
t117 = qJD(1) ^ 2;
t112 = sin(qJ(1));
t116 = cos(qJ(1));
t133 = t112 * g(1) - t116 * g(2);
t84 = t117 * t106 * pkin(8) + qJDD(1) * pkin(1) + t133;
t148 = t108 * t84;
t131 = -t116 * g(1) - t112 * g(2);
t140 = qJDD(1) * t106;
t85 = -t117 * pkin(1) + pkin(8) * t140 + t131;
t128 = -g(3) * t145 + t111 * t148 + t115 * t85;
t142 = qJD(1) * t115;
t135 = t106 * t142;
t86 = (-pkin(2) * t115 - qJ(3) * t111) * t143;
t98 = t100 ^ 2;
t99 = t108 * qJDD(1) + qJDD(2);
t156 = t98 * pkin(2) - t99 * qJ(3) + t100 * t158 - t86 * t135 - t128;
t155 = t108 * g(3);
t154 = mrSges(3,1) - mrSges(4,2);
t153 = mrSges(3,3) + mrSges(4,1);
t152 = -pkin(2) - qJ(4);
t110 = sin(qJ(5));
t114 = cos(qJ(5));
t105 = sin(pkin(11));
t107 = cos(pkin(11));
t146 = t106 ^ 2 * t117;
t136 = t115 ^ 2 * t146;
t81 = pkin(3) * t134 - t100 * qJ(4);
t89 = (qJD(2) * t142 + qJDD(1) * t111) * t106;
t90 = -qJD(2) * t134 + t115 * t140;
t32 = -pkin(3) * t136 - t155 - t89 * qJ(3) + t152 * t90 + (-t84 + (-qJ(3) * t100 * t115 - t111 * t81) * qJD(1)) * t106 + t157;
t144 = t106 * t115;
t147 = g(3) * t144 + t111 * t85;
t126 = -t98 * qJ(3) + t86 * t134 + qJDD(3) + t147;
t35 = t89 * pkin(3) + t152 * t99 + (-pkin(3) * t100 * t143 - qJ(4) * t111 * t146 - t148) * t115 + t126;
t74 = t107 * t100 - t105 * t135;
t132 = -0.2e1 * qJD(4) * t74 - t105 * t32 + t107 * t35;
t63 = -t105 * t90 + t107 * t99;
t73 = -t105 * t100 - t107 * t135;
t22 = (t73 * t134 - t63) * pkin(9) + (t73 * t74 + t89) * pkin(4) + t132;
t138 = 0.2e1 * qJD(4) * t73 + t105 * t35 + t107 * t32;
t62 = -t105 * t99 - t107 * t90;
t64 = pkin(4) * t134 - t74 * pkin(9);
t72 = t73 ^ 2;
t24 = -pkin(4) * t72 + t62 * pkin(9) - t64 * t134 + t138;
t151 = t110 * t22 + t114 * t24;
t83 = mrSges(4,1) * t134 + t100 * mrSges(4,2);
t150 = t100 * mrSges(3,1) - mrSges(3,3) * t134 - t83;
t87 = (mrSges(4,2) * t115 - mrSges(4,3) * t111) * t143;
t149 = t87 + (-mrSges(3,1) * t115 + mrSges(3,2) * t111) * t143;
t119 = t90 * pkin(3) - qJ(4) * t136 + t100 * t81 + qJDD(4) - t156;
t109 = sin(qJ(6));
t113 = cos(qJ(6));
t118 = -t62 * pkin(4) - t72 * pkin(9) + t74 * t64 + t119;
t56 = -t110 * t74 + t114 * t73;
t57 = t110 * t73 + t114 * t74;
t46 = -t56 * pkin(5) - t57 * pkin(10);
t78 = qJDD(5) + t89;
t94 = qJD(5) + t134;
t92 = t94 ^ 2;
t19 = -pkin(5) * t92 + pkin(10) * t78 + t56 * t46 + t151;
t39 = -t57 * qJD(5) - t110 * t63 + t114 * t62;
t40 = t56 * qJD(5) + t110 * t62 + t114 * t63;
t20 = t118 + (t57 * t94 - t39) * pkin(5) + (-t56 * t94 - t40) * pkin(10);
t49 = -t109 * t57 + t113 * t94;
t28 = t49 * qJD(6) + t109 * t78 + t113 * t40;
t50 = t109 * t94 + t113 * t57;
t36 = -mrSges(7,1) * t49 + mrSges(7,2) * t50;
t38 = qJDD(6) - t39;
t55 = qJD(6) - t56;
t41 = -mrSges(7,2) * t55 + mrSges(7,3) * t49;
t16 = m(7) * (-t109 * t19 + t113 * t20) - t28 * mrSges(7,3) + t38 * mrSges(7,1) - t50 * t36 + t55 * t41;
t27 = -t50 * qJD(6) - t109 * t40 + t113 * t78;
t42 = mrSges(7,1) * t55 - mrSges(7,3) * t50;
t17 = m(7) * (t109 * t20 + t113 * t19) + t27 * mrSges(7,3) - t38 * mrSges(7,2) + t49 * t36 - t55 * t42;
t51 = -t94 * mrSges(6,2) + t56 * mrSges(6,3);
t52 = t94 * mrSges(6,1) - t57 * mrSges(6,3);
t123 = m(6) * t118 - t39 * mrSges(6,1) + t40 * mrSges(6,2) + t109 * t17 + t113 * t16 - t56 * t51 + t57 * t52;
t60 = -mrSges(5,2) * t134 + t73 * mrSges(5,3);
t61 = mrSges(5,1) * t134 - t74 * mrSges(5,3);
t121 = m(5) * t119 - t62 * mrSges(5,1) + t63 * mrSges(5,2) - t73 * t60 + t74 * t61 + t123;
t120 = -m(4) * t156 + t121;
t11 = -t150 * t100 + m(3) * t128 + t120 + (-mrSges(3,2) + mrSges(4,3)) * t99 + t153 * t90 + t149 * t135;
t130 = -t106 * t84 - t155;
t45 = -t56 * mrSges(6,1) + t57 * mrSges(6,2);
t12 = m(6) * t151 - t78 * mrSges(6,2) + t39 * mrSges(6,3) - t109 * t16 + t113 * t17 + t56 * t45 - t94 * t52;
t129 = -t110 * t24 + t114 * t22;
t122 = m(7) * (-pkin(5) * t78 - pkin(10) * t92 + t57 * t46 - t129) - t27 * mrSges(7,1) + t28 * mrSges(7,2) - t49 * t41 + t50 * t42;
t13 = m(6) * t129 + t78 * mrSges(6,1) - t40 * mrSges(6,3) - t57 * t45 + t94 * t51 - t122;
t58 = -t73 * mrSges(5,1) + t74 * mrSges(5,2);
t8 = m(5) * t132 + t89 * mrSges(5,1) - t63 * mrSges(5,3) + t110 * t12 + t114 * t13 + t60 * t134 - t74 * t58;
t82 = -mrSges(4,1) * t135 - t100 * mrSges(4,3);
t9 = m(5) * t138 - t89 * mrSges(5,2) + t62 * mrSges(5,3) - t110 * t13 + t114 * t12 - t61 * t134 + t73 * t58;
t127 = -t105 * t8 + t107 * t9 + m(4) * (-t90 * pkin(2) + (-t100 * t135 - t89) * qJ(3) + t130 + t157) + t82 * t135 - t89 * mrSges(4,3);
t80 = -t100 * mrSges(3,2) + mrSges(3,3) * t135;
t5 = m(3) * t130 + t89 * mrSges(3,2) - t154 * t90 + (t150 * t111 - t115 * t80) * t143 + t127;
t137 = t115 * t148;
t125 = -m(4) * (-t99 * pkin(2) + t126 - t137) - t105 * t9 - t107 * t8;
t6 = m(3) * (t137 - t147) + t154 * t99 - t153 * t89 + (t80 - t82) * t100 - t149 * t134 + t125;
t139 = t108 * t5 + t11 * t145 + t6 * t144;
t2 = m(2) * t131 - t117 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t115 * t11 - t111 * t6;
t1 = m(2) * t133 + qJDD(1) * mrSges(2,1) - t117 * mrSges(2,2) - t106 * t5 + (t11 * t111 + t115 * t6) * t108;
t3 = [-m(1) * g(1) - t1 * t112 + t116 * t2, t2, t11, t90 * mrSges(4,2) - t83 * t134 + t127, t9, t12, t17; -m(1) * g(2) + t1 * t116 + t112 * t2, t1, t6, -t90 * mrSges(4,1) - t99 * mrSges(4,3) - t100 * t83 - t87 * t135 - t120, t8, t13, t16; (-m(1) - m(2)) * g(3) + t139, -m(2) * g(3) + t139, t5, t89 * mrSges(4,1) + t99 * mrSges(4,2) + t100 * t82 + t87 * t134 - t125, t121, t123, t122;];
f_new  = t3;
