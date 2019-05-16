% Calculate vector of cutting forces with Newton-Euler
% S6PRRRPR4
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-05-05 07:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRRPR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 07:41:31
% EndTime: 2019-05-05 07:41:43
% DurationCPUTime: 5.83s
% Computational Cost: add. (86658->177), mult. (174364->240), div. (0->0), fcn. (126873->14), ass. (0->100)
t101 = sin(qJ(2));
t105 = cos(qJ(2));
t93 = sin(pkin(11));
t96 = cos(pkin(11));
t82 = t93 * g(1) - t96 * g(2);
t97 = cos(pkin(6));
t128 = t82 * t97;
t83 = -t96 * g(1) - t93 * g(2);
t91 = -g(3) + qJDD(1);
t94 = sin(pkin(6));
t131 = -t101 * t83 + t105 * (t91 * t94 + t128);
t103 = cos(qJ(4));
t106 = qJD(3) ^ 2;
t104 = cos(qJ(3));
t123 = t104 * qJD(2);
t100 = sin(qJ(3));
t107 = qJD(2) ^ 2;
t125 = t101 * t94;
t120 = t101 * t128 + t105 * t83 + t91 * t125;
t51 = -t107 * pkin(2) + qJDD(2) * pkin(8) + t120;
t68 = -t94 * t82 + t97 * t91;
t126 = t100 * t68 + t104 * t51;
t79 = (-pkin(3) * t104 - pkin(9) * t100) * qJD(2);
t35 = -t106 * pkin(3) + qJDD(3) * pkin(9) + t123 * t79 + t126;
t122 = qJD(2) * qJD(3);
t117 = t104 * t122;
t50 = -qJDD(2) * pkin(2) - t107 * pkin(8) - t131;
t80 = t100 * qJDD(2) + t117;
t89 = t100 * t122;
t81 = t104 * qJDD(2) - t89;
t40 = (-t80 - t117) * pkin(9) + (-t81 + t89) * pkin(3) + t50;
t99 = sin(qJ(4));
t118 = t103 * t40 - t99 * t35;
t102 = cos(qJ(6));
t124 = qJD(2) * t100;
t76 = t103 * qJD(3) - t124 * t99;
t59 = t76 * qJD(4) + t99 * qJDD(3) + t103 * t80;
t73 = qJDD(4) - t81;
t77 = t99 * qJD(3) + t103 * t124;
t88 = qJD(4) - t123;
t23 = (t76 * t88 - t59) * qJ(5) + (t76 * t77 + t73) * pkin(4) + t118;
t127 = t103 * t35 + t99 * t40;
t58 = -t77 * qJD(4) + t103 * qJDD(3) - t99 * t80;
t66 = t88 * pkin(4) - t77 * qJ(5);
t72 = t76 ^ 2;
t25 = -t72 * pkin(4) + t58 * qJ(5) - t88 * t66 + t127;
t92 = sin(pkin(12));
t95 = cos(pkin(12));
t62 = t92 * t76 + t95 * t77;
t115 = -0.2e1 * qJD(5) * t62 + t95 * t23 - t92 * t25;
t43 = t92 * t58 + t95 * t59;
t61 = t95 * t76 - t92 * t77;
t17 = (t61 * t88 - t43) * pkin(10) + (t61 * t62 + t73) * pkin(5) + t115;
t121 = 0.2e1 * qJD(5) * t61 + t92 * t23 + t95 * t25;
t42 = t95 * t58 - t92 * t59;
t54 = t88 * pkin(5) - t62 * pkin(10);
t60 = t61 ^ 2;
t18 = -t60 * pkin(5) + t42 * pkin(10) - t88 * t54 + t121;
t98 = sin(qJ(6));
t45 = t102 * t61 - t98 * t62;
t28 = t45 * qJD(6) + t102 * t43 + t98 * t42;
t46 = t102 * t62 + t98 * t61;
t32 = -t45 * mrSges(7,1) + t46 * mrSges(7,2);
t87 = qJD(6) + t88;
t38 = -t87 * mrSges(7,2) + t45 * mrSges(7,3);
t71 = qJDD(6) + t73;
t14 = m(7) * (t102 * t17 - t98 * t18) - t28 * mrSges(7,3) + t71 * mrSges(7,1) - t46 * t32 + t87 * t38;
t27 = -t46 * qJD(6) + t102 * t42 - t98 * t43;
t39 = t87 * mrSges(7,1) - t46 * mrSges(7,3);
t15 = m(7) * (t102 * t18 + t98 * t17) + t27 * mrSges(7,3) - t71 * mrSges(7,2) + t45 * t32 - t87 * t39;
t47 = -t61 * mrSges(6,1) + t62 * mrSges(6,2);
t52 = -t88 * mrSges(6,2) + t61 * mrSges(6,3);
t12 = m(6) * t115 + t73 * mrSges(6,1) - t43 * mrSges(6,3) + t102 * t14 + t98 * t15 - t62 * t47 + t88 * t52;
t53 = t88 * mrSges(6,1) - t62 * mrSges(6,3);
t13 = m(6) * t121 - t73 * mrSges(6,2) + t42 * mrSges(6,3) + t102 * t15 - t98 * t14 + t61 * t47 - t88 * t53;
t63 = -t76 * mrSges(5,1) + t77 * mrSges(5,2);
t65 = -t88 * mrSges(5,2) + t76 * mrSges(5,3);
t10 = m(5) * t118 + t73 * mrSges(5,1) - t59 * mrSges(5,3) + t95 * t12 + t92 * t13 - t77 * t63 + t88 * t65;
t67 = t88 * mrSges(5,1) - t77 * mrSges(5,3);
t11 = m(5) * t127 - t73 * mrSges(5,2) + t58 * mrSges(5,3) - t92 * t12 + t95 * t13 + t76 * t63 - t88 * t67;
t84 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t124;
t85 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t123;
t130 = m(4) * t50 - t81 * mrSges(4,1) + t80 * mrSges(4,2) + t103 * t10 + t99 * t11 + (t100 * t84 - t104 * t85) * qJD(2);
t8 = m(3) * t131 + qJDD(2) * mrSges(3,1) - t107 * mrSges(3,2) - t130;
t129 = t105 * t8;
t116 = -t100 * t51 + t104 * t68;
t34 = -qJDD(3) * pkin(3) - t106 * pkin(9) + t79 * t124 - t116;
t109 = -t58 * pkin(4) - t72 * qJ(5) + t77 * t66 + qJDD(5) + t34;
t112 = t27 * mrSges(7,1) + t45 * t38 - m(7) * (-t42 * pkin(5) - t60 * pkin(10) + t62 * t54 + t109) - t28 * mrSges(7,2) - t46 * t39;
t110 = -m(6) * t109 + t42 * mrSges(6,1) - t43 * mrSges(6,2) + t61 * t52 - t62 * t53 + t112;
t108 = m(5) * t34 - t58 * mrSges(5,1) + t59 * mrSges(5,2) - t76 * t65 + t77 * t67 - t110;
t78 = (-mrSges(4,1) * t104 + mrSges(4,2) * t100) * qJD(2);
t16 = m(4) * t116 + qJDD(3) * mrSges(4,1) - t80 * mrSges(4,3) + qJD(3) * t85 - t124 * t78 - t108;
t9 = m(4) * t126 - qJDD(3) * mrSges(4,2) + t81 * mrSges(4,3) - qJD(3) * t84 - t99 * t10 + t103 * t11 + t123 * t78;
t4 = m(3) * t120 - t107 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t100 * t16 + t104 * t9;
t6 = m(3) * t68 + t100 * t9 + t104 * t16;
t119 = m(2) * t91 + t4 * t125 + t129 * t94 + t97 * t6;
t2 = m(2) * t83 - t101 * t8 + t105 * t4;
t1 = m(2) * t82 - t94 * t6 + (t101 * t4 + t129) * t97;
t3 = [-m(1) * g(1) - t93 * t1 + t96 * t2, t2, t4, t9, t11, t13, t15; -m(1) * g(2) + t96 * t1 + t93 * t2, t1, t8, t16, t10, t12, t14; -m(1) * g(3) + t119, t119, t6, t130, t108, -t110, -t112;];
f_new  = t3;
