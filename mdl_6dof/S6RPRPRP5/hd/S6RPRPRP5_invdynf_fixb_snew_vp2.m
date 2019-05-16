% Calculate vector of cutting forces with Newton-Euler
% S6RPRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-05-05 17:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRPRP5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:46:57
% EndTime: 2019-05-05 17:47:06
% DurationCPUTime: 3.27s
% Computational Cost: add. (37747->188), mult. (92725->240), div. (0->0), fcn. (68922->10), ass. (0->96)
t137 = cos(qJ(3));
t95 = sin(pkin(9));
t97 = cos(pkin(9));
t99 = sin(qJ(3));
t141 = -t97 * t137 + t95 * t99;
t103 = qJD(1) ^ 2;
t93 = t97 ^ 2;
t129 = t95 ^ 2 + t93;
t140 = t129 * mrSges(3,3);
t115 = -t97 * mrSges(3,1) + t95 * mrSges(3,2);
t113 = qJDD(1) * mrSges(3,3) + t103 * t115;
t125 = qJD(1) * qJD(2);
t121 = -t97 * g(3) - 0.2e1 * t95 * t125;
t102 = qJD(3) ^ 2;
t128 = pkin(7) * qJDD(1);
t135 = pkin(2) * t103;
t100 = sin(qJ(1));
t101 = cos(qJ(1));
t114 = -t101 * g(1) - t100 * g(2);
t85 = -t103 * pkin(1) + qJDD(1) * qJ(2) + t114;
t54 = (t97 * t135 - t128 - t85) * t95 + t121;
t119 = -t95 * g(3) + (0.2e1 * t125 + t85) * t97;
t60 = t97 * t128 - t93 * t135 + t119;
t117 = t137 * t54 - t99 * t60;
t83 = t141 * qJD(1);
t111 = t137 * t95 + t97 * t99;
t84 = t111 * qJD(1);
t64 = t83 * pkin(3) - t84 * qJ(4);
t28 = -qJDD(3) * pkin(3) - t102 * qJ(4) + t84 * t64 + qJDD(4) - t117;
t94 = sin(pkin(10));
t96 = cos(pkin(10));
t76 = t94 * qJD(3) + t96 * t84;
t57 = t83 * pkin(4) - t76 * pkin(8);
t127 = t83 * qJD(3);
t70 = t111 * qJDD(1) - t127;
t58 = t96 * qJDD(3) - t94 * t70;
t75 = t96 * qJD(3) - t94 * t84;
t74 = t75 ^ 2;
t105 = -t58 * pkin(4) - t74 * pkin(8) + t76 * t57 + t28;
t136 = cos(qJ(5));
t98 = sin(qJ(5));
t47 = t136 * t76 + t98 * t75;
t59 = t94 * qJDD(3) + t96 * t70;
t29 = t47 * qJD(5) - t136 * t58 + t98 * t59;
t46 = -t136 * t75 + t98 * t76;
t30 = -t46 * qJD(5) + t136 * t59 + t98 * t58;
t81 = qJD(5) + t83;
t41 = -t46 * mrSges(7,2) + t81 * mrSges(7,3);
t44 = -t81 * mrSges(7,1) + t47 * mrSges(7,2);
t110 = t30 * mrSges(7,3) + t47 * t44 - m(7) * (-0.2e1 * qJD(6) * t47 + (t46 * t81 - t30) * qJ(6) + (t47 * t81 + t29) * pkin(5) + t105) - t29 * mrSges(7,1) - t46 * t41;
t42 = -t81 * mrSges(6,2) - t46 * mrSges(6,3);
t43 = t81 * mrSges(6,1) - t47 * mrSges(6,3);
t107 = m(6) * t105 + t29 * mrSges(6,1) + t30 * mrSges(6,2) + t46 * t42 + t47 * t43 - t110;
t55 = -t83 * mrSges(5,2) + t75 * mrSges(5,3);
t56 = t83 * mrSges(5,1) - t76 * mrSges(5,3);
t104 = m(5) * t28 - t58 * mrSges(5,1) + t59 * mrSges(5,2) - t75 * t55 + t76 * t56 + t107;
t65 = t83 * mrSges(4,1) + t84 * mrSges(4,2);
t77 = -qJD(3) * mrSges(4,2) - t83 * mrSges(4,3);
t14 = m(4) * t117 + qJDD(3) * mrSges(4,1) - t70 * mrSges(4,3) + qJD(3) * t77 - t84 * t65 - t104;
t130 = t137 * t60 + t99 * t54;
t31 = -t102 * pkin(3) + qJDD(3) * qJ(4) - t83 * t64 + t130;
t120 = t100 * g(1) - t101 * g(2);
t116 = qJDD(2) - t120;
t106 = (-pkin(2) * t97 - pkin(1)) * qJDD(1) + (-t129 * pkin(7) - qJ(2)) * t103 + t116;
t126 = t84 * qJD(3);
t69 = t141 * qJDD(1) + t126;
t34 = (-t70 + t127) * qJ(4) + (t69 + t126) * pkin(3) + t106;
t118 = -0.2e1 * qJD(4) * t76 - t94 * t31 + t96 * t34;
t20 = (t75 * t83 - t59) * pkin(8) + (t75 * t76 + t69) * pkin(4) + t118;
t123 = 0.2e1 * qJD(4) * t75 + t96 * t31 + t94 * t34;
t22 = -t74 * pkin(4) + t58 * pkin(8) - t83 * t57 + t123;
t132 = t136 * t22 + t98 * t20;
t37 = t46 * pkin(5) - t47 * qJ(6);
t67 = qJDD(5) + t69;
t80 = t81 ^ 2;
t124 = m(7) * (-t80 * pkin(5) + t67 * qJ(6) + 0.2e1 * qJD(6) * t81 - t46 * t37 + t132) + t81 * t44 + t67 * mrSges(7,3);
t38 = t46 * mrSges(7,1) - t47 * mrSges(7,3);
t131 = -t46 * mrSges(6,1) - t47 * mrSges(6,2) - t38;
t133 = -mrSges(6,3) - mrSges(7,2);
t12 = m(6) * t132 - t67 * mrSges(6,2) + t131 * t46 + t133 * t29 - t81 * t43 + t124;
t112 = t136 * t20 - t98 * t22;
t138 = m(7) * (-t67 * pkin(5) - t80 * qJ(6) + t47 * t37 + qJDD(6) - t112);
t13 = m(6) * t112 - t138 + (t42 + t41) * t81 + (mrSges(6,1) + mrSges(7,1)) * t67 + t131 * t47 + t133 * t30;
t48 = -t75 * mrSges(5,1) + t76 * mrSges(5,2);
t10 = m(5) * t118 + t69 * mrSges(5,1) - t59 * mrSges(5,3) + t98 * t12 + t136 * t13 - t76 * t48 + t83 * t55;
t11 = m(5) * t123 - t69 * mrSges(5,2) + t58 * mrSges(5,3) + t136 * t12 - t98 * t13 + t75 * t48 - t83 * t56;
t78 = qJD(3) * mrSges(4,1) - t84 * mrSges(4,3);
t7 = m(4) * t130 - qJDD(3) * mrSges(4,2) - t69 * mrSges(4,3) - qJD(3) * t78 - t94 * t10 + t96 * t11 - t83 * t65;
t4 = m(3) * t121 + t99 * t7 + t137 * t14 + (-m(3) * t85 - t113) * t95;
t5 = m(3) * t119 + t113 * t97 + t137 * t7 - t99 * t14;
t139 = t97 * t4 + t95 * t5;
t109 = m(4) * t106 + t69 * mrSges(4,1) + t70 * mrSges(4,2) + t96 * t10 + t94 * t11 + t83 * t77 + t84 * t78;
t108 = m(3) * (-qJDD(1) * pkin(1) - t103 * qJ(2) + t116) + t109;
t6 = m(2) * t120 + (-mrSges(2,2) + t140) * t103 + (mrSges(2,1) - t115) * qJDD(1) - t108;
t1 = m(2) * t114 - t103 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t95 * t4 + t97 * t5;
t2 = [-m(1) * g(1) + t101 * t1 - t100 * t6, t1, t5, t7, t11, t12, -t29 * mrSges(7,2) - t46 * t38 + t124; -m(1) * g(2) + t100 * t1 + t101 * t6, t6, t4, t14, t10, t13, -t110; (-m(1) - m(2)) * g(3) + t139, -m(2) * g(3) + t139, t115 * qJDD(1) - t103 * t140 + t108, t109, t104, t107, -t67 * mrSges(7,1) + t30 * mrSges(7,2) + t47 * t38 - t81 * t41 + t138;];
f_new  = t2;
