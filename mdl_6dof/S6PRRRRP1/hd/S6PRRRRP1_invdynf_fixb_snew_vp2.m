% Calculate vector of cutting forces with Newton-Euler
% S6PRRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-05-05 09:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRRRP1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:27:45
% EndTime: 2019-05-05 09:27:52
% DurationCPUTime: 2.62s
% Computational Cost: add. (35255->175), mult. (68901->229), div. (0->0), fcn. (49152->12), ass. (0->91)
t100 = cos(qJ(3));
t95 = sin(qJ(4));
t96 = sin(qJ(3));
t99 = cos(qJ(4));
t70 = (t100 * t99 - t95 * t96) * qJD(2);
t101 = cos(qJ(2));
t90 = sin(pkin(11));
t92 = cos(pkin(11));
t78 = t90 * g(1) - t92 * g(2);
t93 = cos(pkin(6));
t130 = t78 * t93;
t79 = -t92 * g(1) - t90 * g(2);
t89 = -g(3) + qJDD(1);
t91 = sin(pkin(6));
t97 = sin(qJ(2));
t133 = (t89 * t91 + t130) * t101 - t97 * t79;
t102 = qJD(2) ^ 2;
t105 = -qJDD(2) * pkin(2) - t133;
t122 = qJD(2) * t96;
t120 = qJD(2) * qJD(3);
t77 = t100 * qJDD(2) - t96 * t120;
t83 = qJD(3) * pkin(3) - pkin(9) * t122;
t88 = t100 ^ 2;
t103 = -t77 * pkin(3) + t83 * t122 + (-pkin(9) * t88 - pkin(8)) * t102 + t105;
t111 = t100 * t120;
t129 = t91 * t97;
t116 = t101 * t79 + t89 * t129 + t97 * t130;
t56 = -t102 * pkin(2) + qJDD(2) * pkin(8) + t116;
t67 = -t91 * t78 + t93 * t89;
t112 = t100 * t67 - t96 * t56;
t76 = t96 * qJDD(2) + t111;
t30 = (-t76 + t111) * pkin(9) + (t100 * t102 * t96 + qJDD(3)) * pkin(3) + t112;
t125 = t100 * t56 + t96 * t67;
t31 = -t88 * t102 * pkin(3) + t77 * pkin(9) - qJD(3) * t83 + t125;
t126 = t95 * t30 + t99 * t31;
t71 = (t100 * t95 + t96 * t99) * qJD(2);
t59 = -t70 * pkin(4) - t71 * pkin(10);
t87 = qJD(3) + qJD(4);
t85 = t87 ^ 2;
t86 = qJDD(3) + qJDD(4);
t23 = -t85 * pkin(4) + t86 * pkin(10) + t70 * t59 + t126;
t48 = -t71 * qJD(4) - t95 * t76 + t99 * t77;
t49 = t70 * qJD(4) + t99 * t76 + t95 * t77;
t26 = (-t70 * t87 - t49) * pkin(10) + (t71 * t87 - t48) * pkin(4) + t103;
t94 = sin(qJ(5));
t98 = cos(qJ(5));
t114 = -t94 * t23 + t98 * t26;
t61 = -t94 * t71 + t98 * t87;
t35 = t61 * qJD(5) + t98 * t49 + t94 * t86;
t46 = qJDD(5) - t48;
t68 = qJD(5) - t70;
t50 = -t68 * mrSges(7,2) + t61 * mrSges(7,3);
t62 = t98 * t71 + t94 * t87;
t119 = m(7) * (-0.2e1 * qJD(6) * t62 + (t61 * t68 - t35) * qJ(6) + (t61 * t62 + t46) * pkin(5) + t114) + t68 * t50 + t46 * mrSges(7,1);
t41 = -t61 * mrSges(7,1) + t62 * mrSges(7,2);
t42 = -t61 * mrSges(6,1) + t62 * mrSges(6,2);
t51 = -t68 * mrSges(6,2) + t61 * mrSges(6,3);
t13 = m(6) * t114 + t46 * mrSges(6,1) + t68 * t51 + (-t42 - t41) * t62 + (-mrSges(6,3) - mrSges(7,3)) * t35 + t119;
t127 = t98 * t23 + t94 * t26;
t34 = -t62 * qJD(5) - t94 * t49 + t98 * t86;
t52 = t68 * pkin(5) - t62 * qJ(6);
t60 = t61 ^ 2;
t118 = m(7) * (-t60 * pkin(5) + t34 * qJ(6) + 0.2e1 * qJD(6) * t61 - t68 * t52 + t127) + t34 * mrSges(7,3) + t61 * t41;
t53 = t68 * mrSges(7,1) - t62 * mrSges(7,3);
t54 = t68 * mrSges(6,1) - t62 * mrSges(6,3);
t15 = m(6) * t127 + t34 * mrSges(6,3) + t61 * t42 + (-t54 - t53) * t68 + (-mrSges(6,2) - mrSges(7,2)) * t46 + t118;
t65 = -t87 * mrSges(5,2) + t70 * mrSges(5,3);
t66 = t87 * mrSges(5,1) - t71 * mrSges(5,3);
t107 = -m(5) * t103 + t48 * mrSges(5,1) - t49 * mrSges(5,2) - t98 * t13 - t94 * t15 + t70 * t65 - t71 * t66;
t80 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t122;
t121 = qJD(2) * t100;
t81 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t121;
t132 = -(t100 * t81 - t96 * t80) * qJD(2) + m(4) * (-t102 * pkin(8) + t105) - t77 * mrSges(4,1) + t76 * mrSges(4,2) - t107;
t113 = t99 * t30 - t95 * t31;
t22 = -t86 * pkin(4) - t85 * pkin(10) + t71 * t59 - t113;
t117 = m(7) * (-t34 * pkin(5) - t60 * qJ(6) + t62 * t52 + qJDD(6) + t22) + t35 * mrSges(7,2) + t62 * t53;
t131 = m(6) * t22 + t35 * mrSges(6,2) - (t51 + t50) * t61 - (mrSges(6,1) + mrSges(7,1)) * t34 + t62 * t54 + t117;
t10 = m(3) * t133 + qJDD(2) * mrSges(3,1) - t102 * mrSges(3,2) - t132;
t123 = t10 * t101;
t58 = -t70 * mrSges(5,1) + t71 * mrSges(5,2);
t11 = m(5) * t126 - t86 * mrSges(5,2) + t48 * mrSges(5,3) - t94 * t13 + t98 * t15 + t70 * t58 - t87 * t66;
t16 = m(5) * t113 + t86 * mrSges(5,1) - t49 * mrSges(5,3) - t71 * t58 + t87 * t65 - t131;
t75 = (-mrSges(4,1) * t100 + mrSges(4,2) * t96) * qJD(2);
t7 = m(4) * t112 + qJDD(3) * mrSges(4,1) - t76 * mrSges(4,3) + qJD(3) * t81 + t95 * t11 - t75 * t122 + t99 * t16;
t8 = m(4) * t125 - qJDD(3) * mrSges(4,2) + t77 * mrSges(4,3) - qJD(3) * t80 + t99 * t11 + t75 * t121 - t95 * t16;
t4 = m(3) * t116 - t102 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t100 * t8 - t96 * t7;
t6 = m(3) * t67 + t100 * t7 + t96 * t8;
t115 = m(2) * t89 + t91 * t123 + t4 * t129 + t93 * t6;
t2 = m(2) * t79 - t97 * t10 + t101 * t4;
t1 = m(2) * t78 - t91 * t6 + (t4 * t97 + t123) * t93;
t3 = [-m(1) * g(1) - t90 * t1 + t92 * t2, t2, t4, t8, t11, t15, -t46 * mrSges(7,2) - t68 * t53 + t118; -m(1) * g(2) + t92 * t1 + t90 * t2, t1, t10, t7, t16, t13, -t35 * mrSges(7,3) - t62 * t41 + t119; -m(1) * g(3) + t115, t115, t6, t132, -t107, t131, -t34 * mrSges(7,1) - t61 * t50 + t117;];
f_new  = t3;
