% Calculate vector of cutting forces with Newton-Euler
% S6RRPRRP11
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-05-06 18:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRRP11_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP11_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP11_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP11_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP11_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP11_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP11_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:52:49
% EndTime: 2019-05-06 18:52:55
% DurationCPUTime: 1.91s
% Computational Cost: add. (19510->201), mult. (40090->243), div. (0->0), fcn. (24519->8), ass. (0->93)
t152 = -2 * qJD(3);
t107 = qJD(1) ^ 2;
t106 = qJD(2) ^ 2;
t100 = sin(qJ(2));
t104 = cos(qJ(2));
t101 = sin(qJ(1));
t105 = cos(qJ(1));
t119 = -t105 * g(1) - t101 * g(2);
t68 = -t107 * pkin(1) + qJDD(1) * pkin(7) + t119;
t124 = -t100 * g(3) + t104 * t68;
t131 = qJD(1) * t104;
t76 = (-pkin(2) * t104 - qJ(3) * t100) * qJD(1);
t149 = t106 * pkin(2) - qJDD(2) * qJ(3) + (qJD(2) * t152) - t76 * t131 - t124;
t130 = qJD(1) * qJD(2);
t121 = t100 * t130;
t80 = t104 * qJDD(1) - t121;
t92 = t100 * qJD(1);
t86 = pkin(3) * t92 - (qJD(2) * pkin(8));
t97 = t104 ^ 2;
t111 = -t97 * t107 * pkin(8) + t80 * pkin(3) + qJD(2) * t86 - t149;
t103 = cos(qJ(4));
t99 = sin(qJ(4));
t75 = t103 * qJD(2) - t99 * t131;
t53 = -t75 * qJD(4) - t99 * qJDD(2) - t103 * t80;
t89 = t92 + qJD(4);
t61 = t89 * pkin(4) - t75 * pkin(9);
t74 = -t99 * qJD(2) - t103 * t131;
t72 = t74 ^ 2;
t109 = -t53 * pkin(4) - t72 * pkin(9) + t75 * t61 + t111;
t102 = cos(qJ(5));
t54 = t74 * qJD(4) + t103 * qJDD(2) - t99 * t80;
t98 = sin(qJ(5));
t57 = t102 * t75 + t98 * t74;
t26 = -t57 * qJD(5) + t102 * t53 - t98 * t54;
t56 = t102 * t74 - t98 * t75;
t27 = t56 * qJD(5) + t102 * t54 + t98 * t53;
t87 = qJD(5) + t89;
t47 = t87 * pkin(5) - t57 * qJ(6);
t48 = t87 * mrSges(7,1) - t57 * mrSges(7,3);
t55 = t56 ^ 2;
t126 = m(7) * (-t26 * pkin(5) - t55 * qJ(6) + t57 * t47 + qJDD(6) + t109) + t27 * mrSges(7,2) + t57 * t48;
t45 = -t87 * mrSges(7,2) + t56 * mrSges(7,3);
t46 = -t87 * mrSges(6,2) + t56 * mrSges(6,3);
t49 = t87 * mrSges(6,1) - t57 * mrSges(6,3);
t112 = m(6) * t109 + t27 * mrSges(6,2) - (t45 + t46) * t56 - (mrSges(6,1) + mrSges(7,1)) * t26 + t57 * t49 + t126;
t59 = -t89 * mrSges(5,2) + t74 * mrSges(5,3);
t60 = t89 * mrSges(5,1) - t75 * mrSges(5,3);
t108 = m(5) * t111 - t53 * mrSges(5,1) + t54 * mrSges(5,2) - t74 * t59 + t75 * t60 + t112;
t151 = -m(4) * t149 + t108;
t125 = t101 * g(1) - t105 * g(2);
t116 = -qJDD(1) * pkin(1) - t125;
t120 = t104 * t130;
t79 = t100 * qJDD(1) + t120;
t110 = pkin(2) * t121 + t92 * t152 + (-t79 - t120) * qJ(3) + t116;
t144 = t107 * pkin(7);
t30 = -t86 * t92 + (-pkin(2) - pkin(8)) * t80 + (-pkin(3) * t97 - pkin(7)) * t107 + t110;
t135 = -t104 * g(3) - t100 * t68;
t44 = -qJDD(2) * pkin(2) - t106 * qJ(3) + t76 * t92 + qJDD(3) - t135;
t35 = (-t100 * t104 * t107 - qJDD(2)) * pkin(8) + (t79 - t120) * pkin(3) + t44;
t122 = t103 * t35 - t99 * t30;
t58 = -t74 * mrSges(5,1) + t75 * mrSges(5,2);
t73 = qJDD(4) + t79;
t18 = (t74 * t89 - t54) * pkin(9) + (t74 * t75 + t73) * pkin(4) + t122;
t137 = t103 * t30 + t99 * t35;
t20 = -t72 * pkin(4) + t53 * pkin(9) - t89 * t61 + t137;
t123 = t102 * t18 - t98 * t20;
t69 = qJDD(5) + t73;
t128 = m(7) * (-0.2e1 * qJD(6) * t57 + (t56 * t87 - t27) * qJ(6) + (t56 * t57 + t69) * pkin(5) + t123) + t87 * t45 + t69 * mrSges(7,1);
t38 = -t56 * mrSges(7,1) + t57 * mrSges(7,2);
t39 = -t56 * mrSges(6,1) + t57 * mrSges(6,2);
t8 = m(6) * t123 + t69 * mrSges(6,1) + t87 * t46 + (-t39 - t38) * t57 + (-mrSges(6,3) - mrSges(7,3)) * t27 + t128;
t138 = t102 * t20 + t98 * t18;
t127 = m(7) * (-t55 * pkin(5) + t26 * qJ(6) + 0.2e1 * qJD(6) * t56 - t87 * t47 + t138) + t26 * mrSges(7,3) + t56 * t38;
t9 = m(6) * t138 + t26 * mrSges(6,3) + t56 * t39 + (-t49 - t48) * t87 + (-mrSges(6,2) - mrSges(7,2)) * t69 + t127;
t6 = m(5) * t122 + t73 * mrSges(5,1) - t54 * mrSges(5,3) + t102 * t8 - t75 * t58 + t89 * t59 + t98 * t9;
t7 = m(5) * t137 - t73 * mrSges(5,2) + t53 * mrSges(5,3) + t102 * t9 + t74 * t58 - t89 * t60 - t98 * t8;
t84 = -mrSges(4,1) * t131 - (qJD(2) * mrSges(4,3));
t117 = -t103 * t7 + t99 * t6 - m(4) * (-t80 * pkin(2) + t110 - t144) - t84 * t131 + t79 * mrSges(4,3);
t85 = mrSges(4,1) * t92 + (qJD(2) * mrSges(4,2));
t133 = (qJD(2) * mrSges(3,1)) - mrSges(3,3) * t92 - t85;
t141 = mrSges(3,1) - mrSges(4,2);
t83 = -(qJD(2) * mrSges(3,2)) + mrSges(3,3) * t131;
t150 = (t133 * t100 - t104 * t83) * qJD(1) - t141 * t80 + m(3) * (t116 - t144) + t79 * mrSges(3,2) - t117;
t77 = (mrSges(4,2) * t104 - mrSges(4,3) * t100) * qJD(1);
t134 = t77 + (-mrSges(3,1) * t104 + mrSges(3,2) * t100) * qJD(1);
t139 = -mrSges(3,3) - mrSges(4,1);
t11 = -t139 * t80 + (-mrSges(3,2) + mrSges(4,3)) * qJDD(2) - t133 * qJD(2) + m(3) * t124 + t134 * t131 + t151;
t115 = -m(4) * t44 - t103 * t6 - t99 * t7;
t4 = m(3) * t135 + t139 * t79 + t141 * qJDD(2) + (t83 - t84) * qJD(2) - t134 * t92 + t115;
t145 = t100 * t11 + t104 * t4;
t2 = m(2) * t125 + qJDD(1) * mrSges(2,1) - t107 * mrSges(2,2) - t150;
t1 = m(2) * t119 - t107 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t100 * t4 + t104 * t11;
t3 = [-m(1) * g(1) + t105 * t1 - t101 * t2, t1, t11, t80 * mrSges(4,2) - t85 * t92 - t117, t7, t9, -t69 * mrSges(7,2) - t87 * t48 + t127; -m(1) * g(2) + t101 * t1 + t105 * t2, t2, t4, -t80 * mrSges(4,1) - qJDD(2) * mrSges(4,3) - qJD(2) * t85 - t77 * t131 - t151, t6, t8, -t27 * mrSges(7,3) - t57 * t38 + t128; (-m(1) - m(2)) * g(3) + t145, -m(2) * g(3) + t145, t150, t79 * mrSges(4,1) + qJDD(2) * mrSges(4,2) + qJD(2) * t84 + t77 * t92 - t115, t108, t112, -t26 * mrSges(7,1) - t56 * t45 + t126;];
f_new  = t3;
