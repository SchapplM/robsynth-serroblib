% Calculate vector of cutting forces with Newton-Euler
% S6RPRRRR9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-05-06 04:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRRR9_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR9_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR9_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR9_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR9_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR9_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR9_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 04:26:21
% EndTime: 2019-05-06 04:26:28
% DurationCPUTime: 3.00s
% Computational Cost: add. (43807->180), mult. (86767->226), div. (0->0), fcn. (58324->10), ass. (0->96)
t93 = sin(qJ(1));
t98 = cos(qJ(1));
t112 = -t98 * g(1) - t93 * g(2);
t128 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t112;
t127 = -m(2) - m(3);
t126 = (-pkin(1) - pkin(7));
t125 = (mrSges(2,1) - mrSges(3,2));
t124 = -mrSges(2,2) + mrSges(3,3);
t100 = qJD(1) ^ 2;
t104 = (t126 * t100) - t128;
t120 = qJD(1) * qJD(3);
t92 = sin(qJ(3));
t115 = t92 * t120;
t97 = cos(qJ(3));
t83 = t97 * t120;
t76 = -t92 * qJDD(1) - t83;
t77 = t97 * qJDD(1) - t115;
t41 = (-t77 + t115) * pkin(8) + (-t76 + t83) * pkin(3) + t104;
t118 = t93 * g(1) - t98 * g(2);
t108 = -t100 * qJ(2) + qJDD(2) - t118;
t61 = t126 * qJDD(1) + t108;
t117 = -t97 * g(3) + t92 * t61;
t75 = (pkin(3) * t92 - pkin(8) * t97) * qJD(1);
t85 = t92 * qJD(1);
t99 = qJD(3) ^ 2;
t44 = -t99 * pkin(3) + qJDD(3) * pkin(8) - t75 * t85 + t117;
t91 = sin(qJ(4));
t96 = cos(qJ(4));
t113 = t96 * t41 - t91 * t44;
t121 = qJD(1) * t97;
t72 = t96 * qJD(3) - t91 * t121;
t50 = t72 * qJD(4) + t91 * qJDD(3) + t96 * t77;
t71 = qJDD(4) - t76;
t73 = t91 * qJD(3) + t96 * t121;
t82 = t85 + qJD(4);
t23 = (t72 * t82 - t50) * pkin(9) + (t72 * t73 + t71) * pkin(4) + t113;
t122 = t91 * t41 + t96 * t44;
t49 = -t73 * qJD(4) + t96 * qJDD(3) - t91 * t77;
t59 = t82 * pkin(4) - t73 * pkin(9);
t70 = t72 ^ 2;
t25 = -t70 * pkin(4) + t49 * pkin(9) - t82 * t59 + t122;
t90 = sin(qJ(5));
t95 = cos(qJ(5));
t123 = t90 * t23 + t95 * t25;
t111 = t92 * g(3) + t97 * t61;
t43 = -qJDD(3) * pkin(3) - t99 * pkin(8) + t75 * t121 - t111;
t102 = -t49 * pkin(4) - t70 * pkin(9) + t73 * t59 + t43;
t53 = t90 * t72 + t95 * t73;
t31 = -t53 * qJD(5) + t95 * t49 - t90 * t50;
t52 = t95 * t72 - t90 * t73;
t32 = t52 * qJD(5) + t90 * t49 + t95 * t50;
t89 = sin(qJ(6));
t94 = cos(qJ(6));
t37 = t89 * t52 + t94 * t53;
t19 = -t37 * qJD(6) + t94 * t31 - t89 * t32;
t36 = t94 * t52 - t89 * t53;
t20 = t36 * qJD(6) + t89 * t31 + t94 * t32;
t81 = qJD(5) + t82;
t78 = qJD(6) + t81;
t33 = -t78 * mrSges(7,2) + t36 * mrSges(7,3);
t34 = t78 * mrSges(7,1) - t37 * mrSges(7,3);
t47 = t81 * pkin(5) - t53 * pkin(10);
t51 = t52 ^ 2;
t107 = t19 * mrSges(7,1) + t36 * t33 - m(7) * (-t31 * pkin(5) - t51 * pkin(10) + t53 * t47 + t102) - t20 * mrSges(7,2) - t37 * t34;
t45 = -t81 * mrSges(6,2) + t52 * mrSges(6,3);
t46 = t81 * mrSges(6,1) - t53 * mrSges(6,3);
t103 = -m(6) * t102 + t31 * mrSges(6,1) - t32 * mrSges(6,2) + t52 * t45 - t53 * t46 + t107;
t55 = -t82 * mrSges(5,2) + t72 * mrSges(5,3);
t56 = t82 * mrSges(5,1) - t73 * mrSges(5,3);
t101 = m(5) * t43 - t49 * mrSges(5,1) + t50 * mrSges(5,2) - t72 * t55 + t73 * t56 - t103;
t74 = (mrSges(4,1) * t92 + mrSges(4,2) * t97) * qJD(1);
t79 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t85;
t11 = m(4) * t111 + qJDD(3) * mrSges(4,1) - t77 * mrSges(4,3) + qJD(3) * t79 - t74 * t121 - t101;
t114 = t95 * t23 - t90 * t25;
t66 = qJDD(5) + t71;
t14 = (t52 * t81 - t32) * pkin(10) + (t52 * t53 + t66) * pkin(5) + t114;
t15 = -t51 * pkin(5) + t31 * pkin(10) - t81 * t47 + t123;
t27 = -t36 * mrSges(7,1) + t37 * mrSges(7,2);
t65 = qJDD(6) + t66;
t12 = m(7) * (t94 * t14 - t89 * t15) - t20 * mrSges(7,3) + t65 * mrSges(7,1) - t37 * t27 + t78 * t33;
t13 = m(7) * (t89 * t14 + t94 * t15) + t19 * mrSges(7,3) - t65 * mrSges(7,2) + t36 * t27 - t78 * t34;
t38 = -t52 * mrSges(6,1) + t53 * mrSges(6,2);
t10 = m(6) * t123 - t66 * mrSges(6,2) + t31 * mrSges(6,3) - t89 * t12 + t94 * t13 + t52 * t38 - t81 * t46;
t54 = -t72 * mrSges(5,1) + t73 * mrSges(5,2);
t9 = m(6) * t114 + t66 * mrSges(6,1) - t32 * mrSges(6,3) + t94 * t12 + t89 * t13 - t53 * t38 + t81 * t45;
t7 = m(5) * t113 + t71 * mrSges(5,1) - t50 * mrSges(5,3) + t90 * t10 - t73 * t54 + t82 * t55 + t95 * t9;
t8 = m(5) * t122 - t71 * mrSges(5,2) + t49 * mrSges(5,3) + t95 * t10 + t72 * t54 - t82 * t56 - t90 * t9;
t80 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t121;
t4 = m(4) * t117 - qJDD(3) * mrSges(4,2) + t76 * mrSges(4,3) - qJD(3) * t80 - t91 * t7 - t74 * t85 + t96 * t8;
t116 = -t92 * t11 + t97 * t4;
t109 = -m(3) * (-qJDD(1) * pkin(1) + t108) - t97 * t11 - t92 * t4;
t106 = m(4) * t104 - t76 * mrSges(4,1) + t77 * mrSges(4,2) + t80 * t121 + t96 * t7 + t79 * t85 + t91 * t8;
t105 = -m(3) * (t100 * pkin(1) + t128) + t106;
t2 = m(2) * t112 + t124 * qJDD(1) - (t125 * t100) + t105;
t1 = m(2) * t118 + t125 * qJDD(1) + t124 * t100 + t109;
t3 = [-m(1) * g(1) - t93 * t1 + t98 * t2, t2, -m(3) * g(3) + t116, t4, t8, t10, t13; -m(1) * g(2) + t98 * t1 + t93 * t2, t1, -(t100 * mrSges(3,2)) - qJDD(1) * mrSges(3,3) - t105, t11, t7, t9, t12; (-m(1) + t127) * g(3) + t116, t127 * g(3) + t116, qJDD(1) * mrSges(3,2) - t100 * mrSges(3,3) - t109, t106, t101, -t103, -t107;];
f_new  = t3;
