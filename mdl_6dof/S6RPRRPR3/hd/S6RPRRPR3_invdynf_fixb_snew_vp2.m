% Calculate vector of cutting forces with Newton-Euler
% S6RPRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-05-05 22:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRPR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:15:11
% EndTime: 2019-05-05 22:15:16
% DurationCPUTime: 1.55s
% Computational Cost: add. (19553->177), mult. (37196->219), div. (0->0), fcn. (23031->10), ass. (0->90)
t93 = sin(qJ(3));
t119 = qJD(1) * t93;
t94 = sin(qJ(1));
t97 = cos(qJ(1));
t115 = t94 * g(1) - t97 * g(2);
t69 = qJDD(1) * pkin(1) + t115;
t109 = -t97 * g(1) - t94 * g(2);
t99 = qJD(1) ^ 2;
t71 = -t99 * pkin(1) + t109;
t89 = sin(pkin(10));
t90 = cos(pkin(10));
t120 = t89 * t69 + t90 * t71;
t42 = -t99 * pkin(2) + qJDD(1) * pkin(7) + t120;
t88 = -g(3) + qJDD(2);
t96 = cos(qJ(3));
t124 = -t93 * t42 + t96 * t88;
t72 = (-pkin(3) * t96 - pkin(8) * t93) * qJD(1);
t98 = qJD(3) ^ 2;
t102 = qJDD(3) * pkin(3) + t98 * pkin(8) - t72 * t119 + t124;
t129 = cos(qJ(4));
t92 = sin(qJ(4));
t67 = -qJD(3) * t129 + t119 * t92;
t118 = t96 * qJD(1);
t80 = qJD(4) - t118;
t128 = t67 * t80;
t117 = qJD(1) * qJD(3);
t113 = t96 * t117;
t73 = t93 * qJDD(1) + t113;
t46 = -t67 * qJD(4) + t92 * qJDD(3) + t129 * t73;
t133 = (-t46 + t128) * qJ(5) - t102;
t111 = t90 * t69 - t89 * t71;
t41 = -qJDD(1) * pkin(2) - t99 * pkin(7) - t111;
t114 = t93 * t117;
t74 = t96 * qJDD(1) - t114;
t27 = (-t73 - t113) * pkin(8) + (-t74 + t114) * pkin(3) + t41;
t123 = t96 * t42 + t93 * t88;
t33 = -t98 * pkin(3) + qJDD(3) * pkin(8) + t118 * t72 + t123;
t125 = t129 * t33 + t92 * t27;
t130 = 2 * qJD(5);
t68 = t92 * qJD(3) + t119 * t129;
t50 = t67 * pkin(4) - t68 * qJ(5);
t66 = qJDD(4) - t74;
t79 = t80 ^ 2;
t105 = -t79 * pkin(4) + t66 * qJ(5) + t130 * t80 - t67 * t50 + t125;
t108 = t129 * t27 - t92 * t33;
t19 = -t66 * pkin(4) - t79 * qJ(5) + t68 * t50 + qJDD(5) - t108;
t14 = (-t46 - t128) * pkin(9) + (t67 * t68 - t66) * pkin(5) + t19;
t45 = t68 * qJD(4) - qJDD(3) * t129 + t92 * t73;
t57 = -t80 * pkin(5) - t68 * pkin(9);
t65 = t67 ^ 2;
t15 = -t65 * pkin(5) + t45 * pkin(9) + t80 * t57 + t105;
t91 = sin(qJ(6));
t95 = cos(qJ(6));
t47 = t95 * t67 - t91 * t68;
t25 = t47 * qJD(6) + t91 * t45 + t95 * t46;
t48 = t91 * t67 + t95 * t68;
t34 = -t47 * mrSges(7,1) + t48 * mrSges(7,2);
t78 = qJD(6) - t80;
t35 = -t78 * mrSges(7,2) + t47 * mrSges(7,3);
t64 = qJDD(6) - t66;
t12 = m(7) * (t95 * t14 - t91 * t15) - t25 * mrSges(7,3) + t64 * mrSges(7,1) - t48 * t34 + t78 * t35;
t24 = -t48 * qJD(6) + t95 * t45 - t91 * t46;
t36 = t78 * mrSges(7,1) - t48 * mrSges(7,3);
t13 = m(7) * (t91 * t14 + t95 * t15) + t24 * mrSges(7,3) - t64 * mrSges(7,2) + t47 * t34 - t78 * t36;
t55 = -t80 * mrSges(6,1) + t68 * mrSges(6,2);
t106 = m(6) * t105 + t66 * mrSges(6,3) - t91 * t12 + t95 * t13 + t80 * t55;
t51 = t67 * mrSges(6,1) - t68 * mrSges(6,3);
t122 = -t67 * mrSges(5,1) - t68 * mrSges(5,2) - t51;
t126 = -mrSges(5,3) - mrSges(6,2);
t54 = t80 * mrSges(5,1) - t68 * mrSges(5,3);
t7 = m(5) * t125 - t66 * mrSges(5,2) + t122 * t67 + t126 * t45 - t80 * t54 + t106;
t75 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t119;
t76 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t118;
t104 = -m(6) * t19 - t95 * t12 - t91 * t13;
t53 = -t80 * mrSges(5,2) - t67 * mrSges(5,3);
t56 = -t67 * mrSges(6,2) + t80 * mrSges(6,3);
t8 = m(5) * t108 + (t53 + t56) * t80 + t122 * t68 + (mrSges(5,1) + mrSges(6,1)) * t66 + t126 * t46 + t104;
t132 = m(4) * t41 - t74 * mrSges(4,1) + t73 * mrSges(4,2) + (t75 * t93 - t76 * t96) * qJD(1) + t129 * t8 + t92 * t7;
t110 = m(7) * (-t65 * pkin(9) + (-pkin(4) - pkin(5)) * t45 + (-pkin(4) * t80 + t130 + t57) * t68 - t133) + t25 * mrSges(7,2) - t24 * mrSges(7,1) + t48 * t36 - t47 * t35;
t103 = m(6) * (-0.2e1 * qJD(5) * t68 + (t68 * t80 + t45) * pkin(4) + t133) + t45 * mrSges(6,1) + t67 * t56 - t110;
t131 = -m(5) * t102 + t45 * mrSges(5,1) + (t54 - t55) * t68 + (mrSges(5,2) - mrSges(6,3)) * t46 + t67 * t53 + t103;
t70 = (-mrSges(4,1) * t96 + mrSges(4,2) * t93) * qJD(1);
t10 = m(4) * t124 + qJDD(3) * mrSges(4,1) - t73 * mrSges(4,3) + qJD(3) * t76 - t119 * t70 - t131;
t6 = m(4) * t123 - qJDD(3) * mrSges(4,2) + t74 * mrSges(4,3) - qJD(3) * t75 + t118 * t70 + t129 * t7 - t92 * t8;
t116 = m(3) * t88 + t96 * t10 + t93 * t6;
t4 = m(3) * t111 + qJDD(1) * mrSges(3,1) - t99 * mrSges(3,2) - t132;
t3 = m(3) * t120 - t99 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t93 * t10 + t96 * t6;
t2 = m(2) * t109 - t99 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t90 * t3 - t89 * t4;
t1 = m(2) * t115 + qJDD(1) * mrSges(2,1) - t99 * mrSges(2,2) + t89 * t3 + t90 * t4;
t5 = [-m(1) * g(1) - t94 * t1 + t97 * t2, t2, t3, t6, t7, -t45 * mrSges(6,2) - t67 * t51 + t106, t13; -m(1) * g(2) + t97 * t1 + t94 * t2, t1, t4, t10, t8, -t46 * mrSges(6,3) - t68 * t55 + t103, t12; (-m(1) - m(2)) * g(3) + t116, -m(2) * g(3) + t116, t116, t132, t131, -t66 * mrSges(6,1) + t46 * mrSges(6,2) + t68 * t51 - t80 * t56 - t104, t110;];
f_new  = t5;
