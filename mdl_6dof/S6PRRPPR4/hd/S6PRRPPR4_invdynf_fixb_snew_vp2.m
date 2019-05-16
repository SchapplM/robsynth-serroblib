% Calculate vector of cutting forces with Newton-Euler
% S6PRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
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
% Datum: 2019-05-05 03:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRPPR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:08:25
% EndTime: 2019-05-05 03:08:30
% DurationCPUTime: 2.01s
% Computational Cost: add. (24067->176), mult. (50178->226), div. (0->0), fcn. (33480->12), ass. (0->94)
t101 = cos(qJ(2));
t92 = sin(pkin(10));
t94 = cos(pkin(10));
t78 = g(1) * t92 - g(2) * t94;
t95 = cos(pkin(6));
t133 = t78 * t95;
t79 = -g(1) * t94 - g(2) * t92;
t90 = -g(3) + qJDD(1);
t93 = sin(pkin(6));
t98 = sin(qJ(2));
t140 = (t90 * t93 + t133) * t101 - t98 * t79;
t102 = qJD(3) ^ 2;
t97 = sin(qJ(3));
t123 = qJD(2) * t97;
t100 = cos(qJ(3));
t103 = qJD(2) ^ 2;
t132 = t93 * t98;
t119 = t101 * t79 + t90 * t132 + t98 * t133;
t40 = -pkin(2) * t103 + qJDD(2) * pkin(8) + t119;
t61 = -t78 * t93 + t90 * t95;
t129 = t100 * t61 - t97 * t40;
t74 = (-pkin(3) * t100 - qJ(4) * t97) * qJD(2);
t106 = qJDD(3) * pkin(3) + qJ(4) * t102 - t74 * t123 - qJDD(4) + t129;
t122 = qJD(2) * t100;
t124 = cos(pkin(11));
t91 = sin(pkin(11));
t68 = -t124 * qJD(3) + t91 * t123;
t118 = t68 * t122;
t121 = qJD(2) * qJD(3);
t116 = t100 * t121;
t76 = qJDD(2) * t97 + t116;
t59 = t91 * qJDD(3) + t124 * t76;
t139 = -(t59 + t118) * qJ(5) - t106;
t136 = -2 * qJD(4);
t128 = t100 * t40 + t97 * t61;
t27 = -pkin(3) * t102 + qJDD(3) * qJ(4) + t74 * t122 + t128;
t39 = -qJDD(2) * pkin(2) - pkin(8) * t103 - t140;
t85 = t97 * t121;
t77 = qJDD(2) * t100 - t85;
t31 = (-t76 - t116) * qJ(4) + (-t77 + t85) * pkin(3) + t39;
t120 = t124 * t27 + t68 * t136 + t91 * t31;
t125 = t103 * t100 ^ 2;
t135 = -2 * qJD(5);
t69 = t91 * qJD(3) + t124 * t123;
t48 = pkin(4) * t68 - qJ(5) * t69;
t107 = -pkin(4) * t125 - t77 * qJ(5) + t122 * t135 - t68 * t48 + t120;
t112 = t124 * t31 - t91 * t27;
t20 = t77 * pkin(4) - qJ(5) * t125 + qJDD(5) - t112 + ((2 * qJD(4)) + t48) * t69;
t16 = (-t59 + t118) * pkin(9) + (t68 * t69 + t77) * pkin(5) + t20;
t58 = -t124 * qJDD(3) + t76 * t91;
t60 = pkin(5) * t122 - pkin(9) * t69;
t66 = t68 ^ 2;
t17 = -pkin(5) * t66 + t58 * pkin(9) - t60 * t122 + t107;
t96 = sin(qJ(6));
t99 = cos(qJ(6));
t46 = t68 * t99 - t69 * t96;
t33 = t46 * qJD(6) + t58 * t96 + t59 * t99;
t47 = t68 * t96 + t69 * t99;
t36 = -mrSges(7,1) * t46 + mrSges(7,2) * t47;
t83 = qJD(6) + t122;
t41 = -mrSges(7,2) * t83 + t46 * mrSges(7,3);
t71 = qJDD(6) + t77;
t14 = m(7) * (t16 * t99 - t17 * t96) - t33 * mrSges(7,3) + t71 * mrSges(7,1) - t47 * t36 + t83 * t41;
t32 = -t47 * qJD(6) + t58 * t99 - t59 * t96;
t42 = mrSges(7,1) * t83 - t47 * mrSges(7,3);
t15 = m(7) * (t16 * t96 + t17 * t99) + t32 * mrSges(7,3) - t71 * mrSges(7,2) + t46 * t36 - t83 * t42;
t110 = m(6) * t107 - t77 * mrSges(6,3) - t96 * t14 + t99 * t15;
t57 = mrSges(6,1) * t122 + mrSges(6,2) * t69;
t126 = -mrSges(5,1) * t122 - mrSges(5,3) * t69 - t57;
t49 = mrSges(6,1) * t68 - mrSges(6,3) * t69;
t127 = -mrSges(5,1) * t68 - mrSges(5,2) * t69 - t49;
t130 = -mrSges(5,3) - mrSges(6,2);
t10 = m(5) * t120 + t77 * mrSges(5,2) + t126 * t122 + t127 * t68 + t130 * t58 + t110;
t109 = -m(6) * t20 - t99 * t14 - t96 * t15;
t54 = -mrSges(6,2) * t68 - mrSges(6,3) * t122;
t55 = mrSges(5,2) * t122 - mrSges(5,3) * t68;
t11 = m(5) * t112 + (-mrSges(5,1) - mrSges(6,1)) * t77 + (m(5) * t136 + t127) * t69 + t130 * t59 + (-t54 - t55) * t122 + t109;
t80 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t123;
t81 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t122;
t138 = m(4) * t39 - t77 * mrSges(4,1) + t76 * mrSges(4,2) - (t100 * t81 - t80 * t97) * qJD(2) + t91 * t10 + t124 * t11;
t115 = m(7) * (-pkin(9) * t66 + (-pkin(4) - pkin(5)) * t58 + (pkin(4) * t122 + (2 * qJD(5)) + t60) * t69 - t139) + t33 * mrSges(7,2) - t32 * mrSges(7,1) + t47 * t42 - t46 * t41;
t108 = m(6) * (t69 * t135 + (-t69 * t122 + t58) * pkin(4) + t139) + t68 * t54 + t58 * mrSges(6,1) - t115;
t137 = -m(5) * t106 + t58 * mrSges(5,1) + t126 * t69 + (mrSges(5,2) - mrSges(6,3)) * t59 + t68 * t55 + t108;
t8 = m(3) * t140 + qJDD(2) * mrSges(3,1) - t103 * mrSges(3,2) - t138;
t134 = t101 * t8;
t75 = (-mrSges(4,1) * t100 + mrSges(4,2) * t97) * qJD(2);
t12 = m(4) * t129 + qJDD(3) * mrSges(4,1) - t76 * mrSges(4,3) + qJD(3) * t81 - t75 * t123 - t137;
t9 = m(4) * t128 - qJDD(3) * mrSges(4,2) + t77 * mrSges(4,3) - qJD(3) * t80 + t124 * t10 - t91 * t11 + t75 * t122;
t4 = m(3) * t119 - t103 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t100 * t9 - t97 * t12;
t6 = m(3) * t61 + t100 * t12 + t9 * t97;
t117 = m(2) * t90 + t4 * t132 + t93 * t134 + t95 * t6;
t2 = m(2) * t79 + t101 * t4 - t8 * t98;
t1 = m(2) * t78 - t6 * t93 + (t4 * t98 + t134) * t95;
t3 = [-m(1) * g(1) - t1 * t92 + t2 * t94, t2, t4, t9, t10, -t58 * mrSges(6,2) - t57 * t122 - t68 * t49 + t110, t15; -m(1) * g(2) + t1 * t94 + t2 * t92, t1, t8, t12, t11, -t59 * mrSges(6,3) - t69 * t57 + t108, t14; -m(1) * g(3) + t117, t117, t6, t138, t137, t77 * mrSges(6,1) + t59 * mrSges(6,2) + t54 * t122 + t69 * t49 - t109, t115;];
f_new  = t3;
