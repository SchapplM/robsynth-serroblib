% Calculate vector of cutting forces with Newton-Euler
% S6RRRPPR6
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
% Datum: 2019-05-07 05:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPPR6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR6_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR6_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR6_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR6_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR6_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 05:26:40
% EndTime: 2019-05-07 05:26:53
% DurationCPUTime: 5.34s
% Computational Cost: add. (67191->215), mult. (149153->282), div. (0->0), fcn. (117622->12), ass. (0->113)
t108 = sin(qJ(2));
t112 = cos(qJ(2));
t104 = sin(pkin(6));
t138 = t104 * t112;
t105 = cos(pkin(6));
t114 = qJD(1) ^ 2;
t109 = sin(qJ(1));
t113 = cos(qJ(1));
t129 = t109 * g(1) - t113 * g(2);
t151 = pkin(8) * t104;
t89 = qJDD(1) * pkin(1) + t114 * t151 + t129;
t141 = t105 * t89;
t127 = -t113 * g(1) - t109 * g(2);
t90 = -t114 * pkin(1) + qJDD(1) * t151 + t127;
t125 = -g(3) * t138 - t108 * t90 + t112 * t141;
t137 = qJD(1) * t104;
t131 = t108 * t137;
t92 = (-pkin(2) * t112 - pkin(9) * t108) * t137;
t100 = t105 * qJD(1) + qJD(2);
t98 = t100 ^ 2;
t99 = t105 * qJDD(1) + qJDD(2);
t52 = -t99 * pkin(2) - t98 * pkin(9) + t92 * t131 - t125;
t107 = sin(qJ(3));
t111 = cos(qJ(3));
t83 = t107 * t100 + t111 * t131;
t135 = qJD(1) * qJD(2);
t93 = (qJDD(1) * t108 + t112 * t135) * t104;
t65 = -t83 * qJD(3) - t107 * t93 + t111 * t99;
t136 = qJD(1) * t112;
t130 = t104 * t136;
t96 = -qJD(3) + t130;
t76 = -t96 * pkin(3) - t83 * qJ(4);
t82 = t111 * t100 - t107 * t131;
t81 = t82 ^ 2;
t117 = -t65 * pkin(3) - t81 * qJ(4) + t83 * t76 + qJDD(4) + t52;
t106 = sin(qJ(6));
t110 = cos(qJ(6));
t103 = sin(pkin(11));
t140 = cos(pkin(11));
t72 = t103 * t83 - t140 * t82;
t149 = t72 * t96;
t152 = -2 * qJD(5);
t66 = t82 * qJD(3) + t107 * t99 + t111 * t93;
t43 = t103 * t65 + t140 * t66;
t73 = t103 * t82 + t140 * t83;
t115 = (-t43 - t149) * qJ(5) + t117 + (-t96 * pkin(4) + t152) * t73;
t142 = t108 * t141 + t112 * t90;
t53 = -t98 * pkin(2) + t99 * pkin(9) + (-g(3) * t108 + t92 * t136) * t104 + t142;
t150 = t105 * g(3);
t94 = (-qJDD(1) * t112 + t108 * t135) * t104;
t54 = t94 * pkin(2) - t93 * pkin(9) - t150 + (-t89 + (pkin(2) * t108 - pkin(9) * t112) * t100 * qJD(1)) * t104;
t128 = -t107 * t53 + t111 * t54;
t86 = qJDD(3) + t94;
t25 = (-t82 * t96 - t66) * qJ(4) + (t82 * t83 + t86) * pkin(3) + t128;
t144 = t107 * t54 + t111 * t53;
t28 = -t81 * pkin(3) + t65 * qJ(4) + t96 * t76 + t144;
t126 = -t103 * t28 + t140 * t25;
t46 = t72 * pkin(4) - t73 * qJ(5);
t157 = (2 * qJD(4)) + t46;
t95 = t96 ^ 2;
t21 = -t86 * pkin(4) - t95 * qJ(5) + t157 * t73 + qJDD(5) - t126;
t16 = (t72 * t73 - t86) * pkin(10) + (t43 - t149) * pkin(5) + t21;
t42 = t103 * t66 - t140 * t65;
t61 = t73 * pkin(5) + t96 * pkin(10);
t71 = t72 ^ 2;
t19 = t115 + (pkin(4) + pkin(10)) * t42 - t73 * t61 - t71 * pkin(5);
t55 = t106 * t96 + t110 * t72;
t33 = t55 * qJD(6) + t106 * t42 + t110 * t86;
t56 = t106 * t72 - t110 * t96;
t36 = -t55 * mrSges(7,1) + t56 * mrSges(7,2);
t70 = qJD(6) + t73;
t37 = -t70 * mrSges(7,2) + t55 * mrSges(7,3);
t41 = qJDD(6) + t43;
t14 = m(7) * (-t106 * t19 + t110 * t16) - t33 * mrSges(7,3) + t41 * mrSges(7,1) - t56 * t36 + t70 * t37;
t32 = -t56 * qJD(6) - t106 * t86 + t110 * t42;
t38 = t70 * mrSges(7,1) - t56 * mrSges(7,3);
t15 = m(7) * (t106 * t16 + t110 * t19) + t32 * mrSges(7,3) - t41 * mrSges(7,2) + t55 * t36 - t70 * t38;
t58 = t73 * mrSges(6,1) - t96 * mrSges(6,2);
t124 = t106 * t14 - t110 * t15 - m(6) * (t42 * pkin(4) + t115) + t43 * mrSges(6,3) + t73 * t58;
t57 = t72 * mrSges(6,1) + t96 * mrSges(6,3);
t143 = -t96 * mrSges(5,2) + t72 * mrSges(5,3) + t57;
t148 = mrSges(5,1) - mrSges(6,2);
t60 = -t96 * mrSges(5,1) - t73 * mrSges(5,3);
t158 = m(5) * t117 + t43 * mrSges(5,2) - t143 * t72 + t148 * t42 + t73 * t60 - t124;
t75 = t96 * mrSges(4,2) + t82 * mrSges(4,3);
t77 = -t96 * mrSges(4,1) - t83 * mrSges(4,3);
t156 = m(4) * t52 - t65 * mrSges(4,1) + t66 * mrSges(4,2) - t82 * t75 + t83 * t77 + t158;
t154 = -2 * qJD(4);
t147 = -mrSges(5,3) - mrSges(6,1);
t146 = t103 * t25 + t140 * t28;
t48 = -t72 * mrSges(6,2) - t73 * mrSges(6,3);
t145 = -t72 * mrSges(5,1) - t73 * mrSges(5,2) - t48;
t139 = t104 * t108;
t88 = -t100 * mrSges(3,2) + mrSges(3,3) * t130;
t91 = (-mrSges(3,1) * t112 + mrSges(3,2) * t108) * t137;
t11 = m(3) * t125 + t99 * mrSges(3,1) - t93 * mrSges(3,3) + t100 * t88 - t91 * t131 - t156;
t123 = t95 * pkin(4) - t86 * qJ(5) - t146;
t68 = t72 * t154;
t120 = -t32 * mrSges(7,1) - t55 * t37 + m(7) * (-t42 * pkin(5) - t71 * pkin(10) - t72 * t46 + t68 + (t152 - t61) * t96 - t123) + t33 * mrSges(7,2) + t56 * t38;
t118 = -m(6) * (0.2e1 * qJD(5) * t96 + t157 * t72 + t123) + t120;
t12 = m(5) * (t68 + t146) + (t60 - t58) * t96 + (-mrSges(5,2) + mrSges(6,3)) * t86 + t145 * t72 + t147 * t42 + t118;
t74 = -t82 * mrSges(4,1) + t83 * mrSges(4,2);
t121 = -m(6) * t21 - t106 * t15 - t110 * t14;
t9 = m(5) * t126 + t143 * t96 + t148 * t86 + (m(5) * t154 + t145) * t73 + t147 * t43 + t121;
t7 = m(4) * t128 + t86 * mrSges(4,1) - t66 * mrSges(4,3) + t103 * t12 + t140 * t9 - t83 * t74 - t96 * t75;
t8 = m(4) * t144 - t86 * mrSges(4,2) + t65 * mrSges(4,3) - t103 * t9 + t140 * t12 + t82 * t74 + t96 * t77;
t87 = t100 * mrSges(3,1) - mrSges(3,3) * t131;
t4 = m(3) * (-g(3) * t139 + t142) - t94 * mrSges(3,3) - t99 * mrSges(3,2) + t91 * t130 - t100 * t87 + t111 * t8 - t107 * t7;
t6 = m(3) * (-t104 * t89 - t150) + t93 * mrSges(3,2) + t94 * mrSges(3,1) + t107 * t8 + t111 * t7 + (t108 * t87 - t112 * t88) * t137;
t134 = t105 * t6 + t11 * t138 + t4 * t139;
t2 = m(2) * t127 - t114 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t108 * t11 + t112 * t4;
t1 = m(2) * t129 + qJDD(1) * mrSges(2,1) - t114 * mrSges(2,2) - t104 * t6 + (t108 * t4 + t112 * t11) * t105;
t3 = [-m(1) * g(1) - t109 * t1 + t113 * t2, t2, t4, t8, t12, -t42 * mrSges(6,2) - t72 * t57 - t124, t15; -m(1) * g(2) + t113 * t1 + t109 * t2, t1, t11, t7, t9, t42 * mrSges(6,1) - t86 * mrSges(6,3) + t72 * t48 + t96 * t58 - t118, t14; (-m(1) - m(2)) * g(3) + t134, -m(2) * g(3) + t134, t6, t156, t158, t43 * mrSges(6,1) + t86 * mrSges(6,2) + t73 * t48 - t96 * t57 - t121, t120;];
f_new  = t3;
