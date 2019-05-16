% Calculate vector of cutting forces with Newton-Euler
% S6RPRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-05-06 02:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRRR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 02:43:49
% EndTime: 2019-05-06 02:43:56
% DurationCPUTime: 3.59s
% Computational Cost: add. (53362->180), mult. (104301->237), div. (0->0), fcn. (72020->12), ass. (0->96)
t104 = cos(qJ(3));
t106 = qJD(1) ^ 2;
t102 = cos(qJ(5));
t100 = sin(qJ(1));
t105 = cos(qJ(1));
t120 = t100 * g(1) - t105 * g(2);
t76 = qJDD(1) * pkin(1) + t120;
t114 = -t105 * g(1) - t100 * g(2);
t78 = -t106 * pkin(1) + t114;
t94 = sin(pkin(11));
t95 = cos(pkin(11));
t117 = t95 * t76 - t94 * t78;
t112 = -qJDD(1) * pkin(2) - t117;
t99 = sin(qJ(3));
t124 = qJD(1) * t99;
t122 = qJD(1) * qJD(3);
t80 = t104 * qJDD(1) - t99 * t122;
t83 = qJD(3) * pkin(3) - pkin(8) * t124;
t92 = t104 ^ 2;
t109 = -t80 * pkin(3) + t83 * t124 + (-pkin(8) * t92 - pkin(7)) * t106 + t112;
t101 = cos(qJ(6));
t103 = cos(qJ(4));
t116 = t104 * t122;
t125 = t94 * t76 + t95 * t78;
t58 = -t106 * pkin(2) + qJDD(1) * pkin(7) + t125;
t93 = -g(3) + qJDD(2);
t118 = t104 * t93 - t99 * t58;
t79 = t99 * qJDD(1) + t116;
t39 = (-t79 + t116) * pkin(8) + (t104 * t106 * t99 + qJDD(3)) * pkin(3) + t118;
t126 = t104 * t58 + t99 * t93;
t40 = -t92 * t106 * pkin(3) + t80 * pkin(8) - qJD(3) * t83 + t126;
t98 = sin(qJ(4));
t127 = t103 * t40 + t98 * t39;
t123 = qJD(1) * t104;
t72 = t103 * t123 - t98 * t124;
t73 = (t103 * t99 + t104 * t98) * qJD(1);
t60 = -t72 * pkin(4) - t73 * pkin(9);
t91 = qJD(3) + qJD(4);
t89 = t91 ^ 2;
t90 = qJDD(3) + qJDD(4);
t26 = -t89 * pkin(4) + t90 * pkin(9) + t72 * t60 + t127;
t50 = -t73 * qJD(4) + t103 * t80 - t98 * t79;
t51 = t72 * qJD(4) + t103 * t79 + t98 * t80;
t29 = (-t72 * t91 - t51) * pkin(9) + (t73 * t91 - t50) * pkin(4) + t109;
t97 = sin(qJ(5));
t119 = t102 * t29 - t97 * t26;
t62 = t102 * t91 - t97 * t73;
t33 = t62 * qJD(5) + t102 * t51 + t97 * t90;
t49 = qJDD(5) - t50;
t63 = t102 * t73 + t97 * t91;
t68 = qJD(5) - t72;
t17 = (t62 * t68 - t33) * pkin(10) + (t62 * t63 + t49) * pkin(5) + t119;
t128 = t102 * t26 + t97 * t29;
t32 = -t63 * qJD(5) + t102 * t90 - t97 * t51;
t54 = t68 * pkin(5) - t63 * pkin(10);
t61 = t62 ^ 2;
t18 = -t61 * pkin(5) + t32 * pkin(10) - t68 * t54 + t128;
t96 = sin(qJ(6));
t43 = t101 * t62 - t96 * t63;
t23 = t43 * qJD(6) + t101 * t33 + t96 * t32;
t44 = t101 * t63 + t96 * t62;
t31 = -t43 * mrSges(7,1) + t44 * mrSges(7,2);
t66 = qJD(6) + t68;
t37 = -t66 * mrSges(7,2) + t43 * mrSges(7,3);
t46 = qJDD(6) + t49;
t15 = m(7) * (t101 * t17 - t96 * t18) - t23 * mrSges(7,3) + t46 * mrSges(7,1) - t44 * t31 + t66 * t37;
t22 = -t44 * qJD(6) + t101 * t32 - t96 * t33;
t38 = t66 * mrSges(7,1) - t44 * mrSges(7,3);
t16 = m(7) * (t101 * t18 + t96 * t17) + t22 * mrSges(7,3) - t46 * mrSges(7,2) + t43 * t31 - t66 * t38;
t45 = -t62 * mrSges(6,1) + t63 * mrSges(6,2);
t52 = -t68 * mrSges(6,2) + t62 * mrSges(6,3);
t12 = m(6) * t119 + t49 * mrSges(6,1) - t33 * mrSges(6,3) + t101 * t15 + t96 * t16 - t63 * t45 + t68 * t52;
t53 = t68 * mrSges(6,1) - t63 * mrSges(6,3);
t13 = m(6) * t128 - t49 * mrSges(6,2) + t32 * mrSges(6,3) + t101 * t16 - t96 * t15 + t62 * t45 - t68 * t53;
t64 = -t91 * mrSges(5,2) + t72 * mrSges(5,3);
t65 = t91 * mrSges(5,1) - t73 * mrSges(5,3);
t110 = -m(5) * t109 + t50 * mrSges(5,1) - t51 * mrSges(5,2) - t102 * t12 - t97 * t13 + t72 * t64 - t73 * t65;
t81 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t124;
t82 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t123;
t129 = -(t104 * t82 - t99 * t81) * qJD(1) + m(4) * (-t106 * pkin(7) + t112) - t80 * mrSges(4,1) + t79 * mrSges(4,2) - t110;
t115 = t103 * t39 - t98 * t40;
t25 = -t90 * pkin(4) - t89 * pkin(9) + t73 * t60 - t115;
t111 = t22 * mrSges(7,1) + t43 * t37 - m(7) * (-t32 * pkin(5) - t61 * pkin(10) + t63 * t54 + t25) - t23 * mrSges(7,2) - t44 * t38;
t107 = m(6) * t25 - t32 * mrSges(6,1) + t33 * mrSges(6,2) - t62 * t52 + t63 * t53 - t111;
t59 = -t72 * mrSges(5,1) + t73 * mrSges(5,2);
t14 = m(5) * t115 + t90 * mrSges(5,1) - t51 * mrSges(5,3) - t73 * t59 + t91 * t64 - t107;
t77 = (-mrSges(4,1) * t104 + mrSges(4,2) * t99) * qJD(1);
t9 = m(5) * t127 - t90 * mrSges(5,2) + t50 * mrSges(5,3) + t102 * t13 - t97 * t12 + t72 * t59 - t91 * t65;
t6 = m(4) * t118 + qJDD(3) * mrSges(4,1) - t79 * mrSges(4,3) + qJD(3) * t82 + t103 * t14 - t77 * t124 + t98 * t9;
t7 = m(4) * t126 - qJDD(3) * mrSges(4,2) + t80 * mrSges(4,3) - qJD(3) * t81 + t103 * t9 + t77 * t123 - t98 * t14;
t121 = m(3) * t93 + t104 * t6 + t99 * t7;
t8 = m(3) * t117 + qJDD(1) * mrSges(3,1) - t106 * mrSges(3,2) - t129;
t3 = m(3) * t125 - t106 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t104 * t7 - t99 * t6;
t2 = m(2) * t114 - t106 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t95 * t3 - t94 * t8;
t1 = m(2) * t120 + qJDD(1) * mrSges(2,1) - t106 * mrSges(2,2) + t94 * t3 + t95 * t8;
t4 = [-m(1) * g(1) - t100 * t1 + t105 * t2, t2, t3, t7, t9, t13, t16; -m(1) * g(2) + t105 * t1 + t100 * t2, t1, t8, t6, t14, t12, t15; (-m(1) - m(2)) * g(3) + t121, -m(2) * g(3) + t121, t121, t129, -t110, t107, -t111;];
f_new  = t4;
