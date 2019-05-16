% Calculate vector of cutting forces with Newton-Euler
% S6RPRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-05-05 21:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRPR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:55:21
% EndTime: 2019-05-05 21:55:29
% DurationCPUTime: 3.63s
% Computational Cost: add. (50859->180), mult. (108760->240), div. (0->0), fcn. (75749->12), ass. (0->95)
t101 = qJD(1) ^ 2;
t100 = cos(qJ(1));
t96 = sin(qJ(1));
t115 = g(1) * t96 - g(2) * t100;
t72 = qJDD(1) * pkin(1) + t115;
t111 = -g(1) * t100 - g(2) * t96;
t74 = -pkin(1) * t101 + t111;
t90 = sin(pkin(10));
t92 = cos(pkin(10));
t112 = t92 * t72 - t74 * t90;
t108 = -qJDD(1) * pkin(2) - t112;
t95 = sin(qJ(3));
t121 = qJD(1) * t95;
t119 = qJD(1) * qJD(3);
t99 = cos(qJ(3));
t76 = qJDD(1) * t99 - t119 * t95;
t79 = qJD(3) * pkin(3) - pkin(8) * t121;
t87 = t99 ^ 2;
t105 = -t76 * pkin(3) + t79 * t121 + (-pkin(8) * t87 - pkin(7)) * t101 + t108;
t94 = sin(qJ(4));
t98 = cos(qJ(4));
t70 = (t94 * t99 + t95 * t98) * qJD(1);
t116 = t99 * t119;
t75 = qJDD(1) * t95 + t116;
t48 = -qJD(4) * t70 - t75 * t94 + t76 * t98;
t86 = qJD(3) + qJD(4);
t63 = pkin(4) * t86 - qJ(5) * t70;
t69 = (-t94 * t95 + t98 * t99) * qJD(1);
t65 = t69 ^ 2;
t103 = -t48 * pkin(4) - t65 * qJ(5) + t63 * t70 + qJDD(5) + t105;
t125 = 2 * qJD(5);
t122 = t72 * t90 + t74 * t92;
t60 = -pkin(2) * t101 + qJDD(1) * pkin(7) + t122;
t88 = -g(3) + qJDD(2);
t113 = -t95 * t60 + t88 * t99;
t41 = (-t75 + t116) * pkin(8) + (t101 * t95 * t99 + qJDD(3)) * pkin(3) + t113;
t123 = t60 * t99 + t88 * t95;
t42 = -pkin(3) * t101 * t87 + pkin(8) * t76 - qJD(3) * t79 + t123;
t114 = t41 * t98 - t94 * t42;
t49 = t69 * qJD(4) + t75 * t98 + t76 * t94;
t85 = qJDD(3) + qJDD(4);
t21 = (t69 * t86 - t49) * qJ(5) + (t69 * t70 + t85) * pkin(4) + t114;
t124 = t41 * t94 + t42 * t98;
t23 = -t65 * pkin(4) + t48 * qJ(5) - t63 * t86 + t124;
t89 = sin(pkin(11));
t91 = cos(pkin(11));
t57 = t69 * t91 - t70 * t89;
t117 = t125 * t57 + t21 * t89 + t23 * t91;
t58 = t69 * t89 + t70 * t91;
t40 = -pkin(5) * t57 - pkin(9) * t58;
t84 = t86 ^ 2;
t18 = -pkin(5) * t84 + pkin(9) * t85 + t57 * t40 + t117;
t30 = t48 * t91 - t49 * t89;
t31 = t48 * t89 + t49 * t91;
t19 = (-t57 * t86 - t31) * pkin(9) + (t58 * t86 - t30) * pkin(5) + t103;
t93 = sin(qJ(6));
t97 = cos(qJ(6));
t46 = -t58 * t93 + t86 * t97;
t27 = t46 * qJD(6) + t31 * t97 + t85 * t93;
t29 = qJDD(6) - t30;
t47 = t58 * t97 + t86 * t93;
t32 = -mrSges(7,1) * t46 + mrSges(7,2) * t47;
t54 = qJD(6) - t57;
t33 = -mrSges(7,2) * t54 + mrSges(7,3) * t46;
t15 = m(7) * (-t18 * t93 + t19 * t97) - t27 * mrSges(7,3) + t29 * mrSges(7,1) - t47 * t32 + t54 * t33;
t26 = -t47 * qJD(6) - t31 * t93 + t85 * t97;
t34 = mrSges(7,1) * t54 - mrSges(7,3) * t47;
t16 = m(7) * (t18 * t97 + t19 * t93) + t26 * mrSges(7,3) - t29 * mrSges(7,2) + t46 * t32 - t54 * t34;
t50 = -mrSges(6,2) * t86 + t57 * mrSges(6,3);
t51 = mrSges(6,1) * t86 - t58 * mrSges(6,3);
t107 = -m(6) * t103 + t30 * mrSges(6,1) - t31 * mrSges(6,2) - t15 * t97 - t16 * t93 + t57 * t50 - t58 * t51;
t62 = -mrSges(5,2) * t86 + t69 * mrSges(5,3);
t64 = mrSges(5,1) * t86 - mrSges(5,3) * t70;
t104 = -m(5) * t105 + t48 * mrSges(5,1) - t49 * mrSges(5,2) + t69 * t62 - t64 * t70 + t107;
t77 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t121;
t120 = qJD(1) * t99;
t78 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t120;
t126 = (t77 * t95 - t78 * t99) * qJD(1) + m(4) * (-t101 * pkin(7) + t108) - t76 * mrSges(4,1) + t75 * mrSges(4,2) - t104;
t73 = (-mrSges(4,1) * t99 + mrSges(4,2) * t95) * qJD(1);
t39 = -mrSges(6,1) * t57 + mrSges(6,2) * t58;
t11 = m(6) * t117 - mrSges(6,2) * t85 + t30 * mrSges(6,3) - t15 * t93 + t16 * t97 + t57 * t39 - t51 * t86;
t110 = -t91 * t21 + t89 * t23;
t106 = m(7) * (-t85 * pkin(5) - t84 * pkin(9) + (t125 + t40) * t58 + t110) - t26 * mrSges(7,1) + t27 * mrSges(7,2) - t46 * t33 + t47 * t34;
t12 = m(6) * (-0.2e1 * qJD(5) * t58 - t110) - t31 * mrSges(6,3) + t85 * mrSges(6,1) - t58 * t39 + t86 * t50 - t106;
t61 = -t69 * mrSges(5,1) + mrSges(5,2) * t70;
t8 = m(5) * t114 + mrSges(5,1) * t85 - t49 * mrSges(5,3) + t11 * t89 + t12 * t91 - t61 * t70 + t62 * t86;
t9 = m(5) * t124 - mrSges(5,2) * t85 + t48 * mrSges(5,3) + t11 * t91 - t12 * t89 + t69 * t61 - t64 * t86;
t6 = m(4) * t113 + qJDD(3) * mrSges(4,1) - t75 * mrSges(4,3) + qJD(3) * t78 - t121 * t73 + t98 * t8 + t94 * t9;
t7 = m(4) * t123 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t76 - qJD(3) * t77 + t120 * t73 - t8 * t94 + t9 * t98;
t118 = m(3) * t88 + t6 * t99 + t7 * t95;
t10 = m(3) * t112 + qJDD(1) * mrSges(3,1) - t101 * mrSges(3,2) - t126;
t3 = m(3) * t122 - mrSges(3,1) * t101 - qJDD(1) * mrSges(3,2) - t6 * t95 + t7 * t99;
t2 = m(2) * t111 - mrSges(2,1) * t101 - qJDD(1) * mrSges(2,2) - t10 * t90 + t3 * t92;
t1 = m(2) * t115 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t101 + t10 * t92 + t3 * t90;
t4 = [-m(1) * g(1) - t1 * t96 + t100 * t2, t2, t3, t7, t9, t11, t16; -m(1) * g(2) + t1 * t100 + t2 * t96, t1, t10, t6, t8, t12, t15; (-m(1) - m(2)) * g(3) + t118, -m(2) * g(3) + t118, t118, t126, -t104, -t107, t106;];
f_new  = t4;
