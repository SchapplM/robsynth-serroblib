% Calculate vector of cutting forces with Newton-Euler
% S6RRRPPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
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
% Datum: 2019-05-07 04:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPPP1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPP1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPP1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPP1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPP1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPP1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPP1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 03:57:15
% EndTime: 2019-05-07 03:57:25
% DurationCPUTime: 3.87s
% Computational Cost: add. (47712->221), mult. (100730->268), div. (0->0), fcn. (73264->10), ass. (0->110)
t146 = cos(qJ(2));
t124 = t146 * qJD(1);
t116 = t124 - qJD(3);
t97 = sin(pkin(6));
t112 = t116 * t97;
t136 = cos(pkin(10));
t137 = cos(pkin(6));
t117 = t137 * t136;
t101 = cos(qJ(3));
t99 = sin(qJ(2));
t135 = qJD(1) * t99;
t98 = sin(qJ(3));
t85 = t101 * qJD(2) - t98 * t135;
t86 = t98 * qJD(2) + t101 * t135;
t96 = sin(pkin(10));
t59 = t136 * t112 - t85 * t117 + t96 * t86;
t109 = t85 * t137 - t112;
t60 = t109 * t96 + t136 * t86;
t34 = t59 * pkin(4) - t60 * qJ(5);
t155 = (2 * qJD(4)) + t34;
t123 = qJ(4) * t137;
t122 = qJD(2) * t124;
t129 = qJD(2) * t135;
t104 = qJD(1) ^ 2;
t100 = sin(qJ(1));
t102 = cos(qJ(1));
t128 = t100 * g(1) - t102 * g(2);
t81 = -qJDD(1) * pkin(1) - t104 * pkin(8) - t128;
t89 = t99 * qJDD(1) + t122;
t90 = t146 * qJDD(1) - t129;
t54 = (-t122 - t89) * pkin(9) + (-t90 + t129) * pkin(2) + t81;
t103 = qJD(2) ^ 2;
t120 = -t102 * g(1) - t100 * g(2);
t82 = -t104 * pkin(1) + qJDD(1) * pkin(8) + t120;
t130 = -t99 * g(3) + t146 * t82;
t88 = (-t146 * pkin(2) - pkin(9) * t99) * qJD(1);
t63 = -t103 * pkin(2) + qJDD(2) * pkin(9) + t88 * t124 + t130;
t127 = t101 * t54 - t98 * t63;
t66 = t85 * qJD(3) + t98 * qJDD(2) + t101 * t89;
t138 = qJ(4) * t97;
t67 = -t85 * pkin(3) - t86 * t138;
t68 = t109 * qJ(4);
t84 = qJDD(3) - t90;
t23 = t84 * pkin(3) - t116 * t68 - t66 * t123 - t86 * t67 + t127;
t139 = -t146 * g(3) - t99 * t82;
t62 = -qJDD(2) * pkin(2) - t103 * pkin(9) + t88 * t135 - t139;
t65 = -t86 * qJD(3) + t101 * qJDD(2) - t98 * t89;
t73 = -t116 * pkin(3) - t86 * t123;
t27 = -t65 * pkin(3) - t66 * t138 - t85 * t68 + t86 * t73 + t62;
t121 = t137 * t27 - t97 * t23 + qJDD(4);
t71 = t137 * t116 + t97 * t85;
t145 = t59 * t71;
t149 = -2 * qJD(5);
t113 = t137 * t65 + t84 * t97;
t41 = t113 * t96 + t136 * t66;
t107 = (-t41 - t145) * qJ(5) + t121 + (-t71 * pkin(4) + t149) * t60;
t125 = t97 * t136;
t40 = -t65 * t117 - t84 * t125 + t96 * t66;
t48 = t60 * mrSges(6,1) - t71 * mrSges(6,2);
t154 = m(6) * (t40 * pkin(4) + t107) - t60 * t48 - t41 * mrSges(6,3);
t126 = t96 * t137;
t140 = t101 * t63 + t98 * t54;
t24 = t113 * qJ(4) + t116 * t73 + t85 * t67 + t140;
t132 = t97 * t96 * t27 + t23 * t126 + t136 * t24;
t53 = t137 * t84 - t97 * t65;
t70 = t71 ^ 2;
t110 = t70 * pkin(4) - t53 * qJ(5) - t132;
t44 = t60 * pkin(5) + t71 * qJ(6);
t45 = t60 * mrSges(7,1) + t71 * mrSges(7,3);
t151 = -2 * qJD(4);
t56 = t59 * t151;
t58 = t59 ^ 2;
t131 = t71 * t45 - m(7) * (-t40 * pkin(5) - t58 * qJ(6) - t59 * t34 + qJDD(6) + t56 + (t149 - t44) * t71 - t110) - t53 * mrSges(7,2);
t118 = m(6) * (0.2e1 * qJD(5) * t71 + t155 * t59 + t110) + t131;
t36 = -t59 * mrSges(6,2) - t60 * mrSges(6,3);
t142 = -t59 * mrSges(5,1) - t60 * mrSges(5,2) - t36;
t143 = -mrSges(5,3) - mrSges(6,1);
t33 = -t60 * mrSges(7,2) + t59 * mrSges(7,3);
t43 = -t71 * mrSges(5,1) - t60 * mrSges(5,3);
t10 = m(5) * (t56 + t132) + (t43 - t48) * t71 + (-mrSges(5,2) + mrSges(6,3)) * t53 + (-t33 + t142) * t59 + (-mrSges(7,1) + t143) * t40 - t118;
t148 = 2 * qJD(6);
t47 = -t59 * mrSges(7,1) - t71 * mrSges(7,2);
t133 = m(7) * (-t58 * pkin(5) + t59 * t148 - t60 * t44 + (pkin(4) + qJ(6)) * t40 + t107) + t59 * t47 + t40 * mrSges(7,3);
t42 = t71 * mrSges(5,2) - t59 * mrSges(5,3);
t46 = t59 * mrSges(6,1) + t71 * mrSges(6,3);
t11 = m(5) * t121 + (t43 - t45) * t60 + (t42 - t46) * t59 + (mrSges(5,2) - mrSges(7,2)) * t41 + (mrSges(5,1) - mrSges(6,2)) * t40 + t133 + t154;
t74 = t116 * mrSges(4,2) + t85 * mrSges(4,3);
t75 = -t116 * mrSges(4,1) - t86 * mrSges(4,3);
t106 = t23 * t117 + t27 * t125 - t96 * t24;
t18 = -t53 * pkin(4) - t70 * qJ(5) + t155 * t60 + qJDD(5) - t106;
t134 = m(7) * (t71 * t148 + (t59 * t60 - t53) * qJ(6) + (t41 - t145) * pkin(5) + t18) + t60 * t33 + t41 * mrSges(7,1);
t119 = m(6) * t18 + t134;
t141 = -t46 + t47;
t144 = mrSges(6,2) - mrSges(7,3);
t9 = m(5) * t106 + (m(5) * t151 + t142) * t60 + t143 * t41 + (-t42 - t141) * t71 + (mrSges(5,1) - t144) * t53 - t119;
t153 = m(4) * t62 - t65 * mrSges(4,1) + t66 * mrSges(4,2) + t137 * t11 + (t10 * t96 + t136 * t9) * t97 - t85 * t74 + t86 * t75;
t69 = -t85 * mrSges(4,1) + t86 * mrSges(4,2);
t7 = m(4) * t127 + t84 * mrSges(4,1) - t66 * mrSges(4,3) + t10 * t126 - t97 * t11 - t116 * t74 + t9 * t117 - t86 * t69;
t8 = m(4) * t140 - t84 * mrSges(4,2) + t65 * mrSges(4,3) + t136 * t10 + t116 * t75 + t85 * t69 - t96 * t9;
t92 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t135;
t93 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t124;
t152 = m(3) * t81 - t90 * mrSges(3,1) + t89 * mrSges(3,2) - (t146 * t93 - t92 * t99) * qJD(1) + t101 * t7 + t98 * t8;
t87 = (-mrSges(3,1) * t146 + mrSges(3,2) * t99) * qJD(1);
t4 = m(3) * t130 - qJDD(2) * mrSges(3,2) + t90 * mrSges(3,3) - qJD(2) * t92 + t101 * t8 + t87 * t124 - t98 * t7;
t6 = m(3) * t139 + qJDD(2) * mrSges(3,1) - t89 * mrSges(3,3) + qJD(2) * t93 - t87 * t135 - t153;
t147 = t146 * t6 + t99 * t4;
t111 = -t41 * mrSges(7,2) - t60 * t45 + t133;
t2 = m(2) * t128 + qJDD(1) * mrSges(2,1) - t104 * mrSges(2,2) - t152;
t1 = m(2) * t120 - t104 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t146 * t4 - t99 * t6;
t3 = [-m(1) * g(1) + t102 * t1 - t100 * t2, t1, t4, t8, t10, -t40 * mrSges(6,2) - t59 * t46 + t111 + t154, t111; -m(1) * g(2) + t100 * t1 + t102 * t2, t2, t6, t7, t9, -t53 * mrSges(6,3) + t71 * t48 + (t33 + t36) * t59 + (mrSges(6,1) + mrSges(7,1)) * t40 + t118, -t53 * mrSges(7,3) + t71 * t47 + t134; (-m(1) - m(2)) * g(3) + t147, -m(2) * g(3) + t147, t152, t153, t11, t41 * mrSges(6,1) + t141 * t71 + t144 * t53 + t60 * t36 + t119, -t40 * mrSges(7,1) - t59 * t33 - t131;];
f_new  = t3;
