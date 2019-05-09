% Calculate vector of cutting forces with Newton-Euler
% S6RPRRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-05-06 01:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRRP8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP8_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP8_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP8_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP8_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP8_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:46:18
% EndTime: 2019-05-06 01:46:21
% DurationCPUTime: 1.39s
% Computational Cost: add. (16644->176), mult. (32549->217), div. (0->0), fcn. (21207->8), ass. (0->85)
t87 = sin(qJ(1));
t90 = cos(qJ(1));
t103 = -t90 * g(1) - t87 * g(2);
t99 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t103;
t85 = sin(qJ(4));
t86 = sin(qJ(3));
t88 = cos(qJ(4));
t89 = cos(qJ(3));
t62 = (t85 * t89 + t86 * t88) * qJD(1);
t112 = qJD(1) * qJD(3);
t105 = t86 * t112;
t126 = -pkin(1) - pkin(7);
t107 = t87 * g(1) - t90 * g(2);
t91 = qJD(1) ^ 2;
t97 = -t91 * qJ(2) + qJDD(2) - t107;
t57 = t126 * qJDD(1) + t97;
t115 = t86 * g(3) + t89 * t57;
t70 = t89 * qJDD(1) - t105;
t29 = (-t70 - t105) * pkin(8) + (-t86 * t89 * t91 + qJDD(3)) * pkin(3) + t115;
t106 = -t89 * g(3) + t86 * t57;
t69 = -t86 * qJDD(1) - t89 * t112;
t113 = qJD(1) * t89;
t73 = qJD(3) * pkin(3) - pkin(8) * t113;
t83 = t86 ^ 2;
t30 = -t83 * t91 * pkin(3) + t69 * pkin(8) - qJD(3) * t73 + t106;
t104 = t88 * t29 - t85 * t30;
t63 = (-t85 * t86 + t88 * t89) * qJD(1);
t48 = t62 * pkin(4) - t63 * pkin(9);
t81 = qJD(3) + qJD(4);
t79 = t81 ^ 2;
t80 = qJDD(3) + qJDD(4);
t19 = -t80 * pkin(4) - t79 * pkin(9) + t63 * t48 - t104;
t124 = cos(qJ(5));
t41 = -t62 * qJD(4) + t85 * t69 + t88 * t70;
t84 = sin(qJ(5));
t50 = t124 * t63 + t84 * t81;
t22 = t50 * qJD(5) - t124 * t80 + t84 * t41;
t49 = -t124 * t81 + t84 * t63;
t23 = -t49 * qJD(5) + t124 * t41 + t84 * t80;
t61 = qJD(5) + t62;
t42 = -t49 * mrSges(7,2) + t61 * mrSges(7,3);
t109 = m(7) * (-0.2e1 * qJD(6) * t50 + (t49 * t61 - t23) * qJ(6) + (t50 * t61 + t22) * pkin(5) + t19) + t22 * mrSges(7,1) + t49 * t42;
t43 = -t61 * mrSges(6,2) - t49 * mrSges(6,3);
t44 = t61 * mrSges(6,1) - t50 * mrSges(6,3);
t45 = -t61 * mrSges(7,1) + t50 * mrSges(7,2);
t128 = m(6) * t19 + t22 * mrSges(6,1) + (t44 - t45) * t50 + (mrSges(6,2) - mrSges(7,3)) * t23 + t49 * t43 + t109;
t127 = -m(2) - m(3);
t40 = -t63 * qJD(4) + t88 * t69 - t85 * t70;
t93 = -t69 * pkin(3) + t73 * t113 + (-pkin(8) * t83 + t126) * t91 + t99;
t17 = (t62 * t81 - t41) * pkin(9) + (t63 * t81 - t40) * pkin(4) + t93;
t118 = t85 * t29 + t88 * t30;
t20 = -t79 * pkin(4) + t80 * pkin(9) - t62 * t48 + t118;
t100 = t124 * t17 - t84 * t20;
t31 = t49 * pkin(5) - t50 * qJ(6);
t39 = qJDD(5) - t40;
t59 = t61 ^ 2;
t125 = m(7) * (-t39 * pkin(5) - t59 * qJ(6) + t50 * t31 + qJDD(6) - t100);
t123 = (mrSges(2,1) - mrSges(3,2));
t122 = -mrSges(2,2) + mrSges(3,3);
t120 = -mrSges(6,3) - mrSges(7,2);
t119 = t124 * t20 + t84 * t17;
t32 = t49 * mrSges(7,1) - t50 * mrSges(7,3);
t117 = -t49 * mrSges(6,1) - t50 * mrSges(6,2) - t32;
t114 = qJD(1) * t86;
t110 = m(7) * (-t59 * pkin(5) + t39 * qJ(6) + 0.2e1 * qJD(6) * t61 - t49 * t31 + t119) + t61 * t45 + t39 * mrSges(7,3);
t11 = m(6) * t100 - t125 + (t43 + t42) * t61 + t117 * t50 + (mrSges(6,1) + mrSges(7,1)) * t39 + t120 * t23;
t47 = t62 * mrSges(5,1) + t63 * mrSges(5,2);
t55 = t81 * mrSges(5,1) - t63 * mrSges(5,3);
t9 = m(6) * t119 - t39 * mrSges(6,2) + t117 * t49 + t120 * t22 - t61 * t44 + t110;
t6 = m(5) * t118 - t80 * mrSges(5,2) + t40 * mrSges(5,3) - t84 * t11 + t124 * t9 - t62 * t47 - t81 * t55;
t68 = (mrSges(4,1) * t86 + mrSges(4,2) * t89) * qJD(1);
t54 = -t81 * mrSges(5,2) - t62 * mrSges(5,3);
t7 = m(5) * t104 + t80 * mrSges(5,1) - t41 * mrSges(5,3) - t63 * t47 + t81 * t54 - t128;
t71 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t114;
t3 = m(4) * t115 + qJDD(3) * mrSges(4,1) - t70 * mrSges(4,3) + qJD(3) * t71 - t68 * t113 + t85 * t6 + t88 * t7;
t72 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t113;
t4 = m(4) * t106 - qJDD(3) * mrSges(4,2) + t69 * mrSges(4,3) - qJD(3) * t72 - t68 * t114 + t88 * t6 - t85 * t7;
t108 = -t86 * t3 + t89 * t4;
t98 = -m(3) * (-qJDD(1) * pkin(1) + t97) - t89 * t3 - t86 * t4;
t96 = m(5) * t93 - t40 * mrSges(5,1) + t41 * mrSges(5,2) + t124 * t11 + t62 * t54 + t63 * t55 + t84 * t9;
t94 = -t69 * mrSges(4,1) + m(4) * (t126 * t91 + t99) + t71 * t114 + t72 * t113 + t70 * mrSges(4,2) + t96;
t92 = -m(3) * (t91 * pkin(1) - t99) + t94;
t5 = m(2) * t103 + t122 * qJDD(1) - (t123 * t91) + t92;
t1 = m(2) * t107 + t123 * qJDD(1) + t122 * t91 + t98;
t2 = [-m(1) * g(1) - t87 * t1 + t90 * t5, t5, -m(3) * g(3) + t108, t4, t6, t9, -t22 * mrSges(7,2) - t49 * t32 + t110; -m(1) * g(2) + t90 * t1 + t87 * t5, t1, -(t91 * mrSges(3,2)) - qJDD(1) * mrSges(3,3) - t92, t3, t7, t11, -t23 * mrSges(7,3) - t50 * t45 + t109; (-m(1) + t127) * g(3) + t108, t127 * g(3) + t108, qJDD(1) * mrSges(3,2) - t91 * mrSges(3,3) - t98, t94, t96, t128, -t39 * mrSges(7,1) + t23 * mrSges(7,2) + t50 * t32 - t61 * t42 + t125;];
f_new  = t2;
