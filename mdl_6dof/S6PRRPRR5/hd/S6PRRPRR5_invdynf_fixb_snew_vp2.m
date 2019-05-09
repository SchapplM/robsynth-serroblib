% Calculate vector of cutting forces with Newton-Euler
% S6PRRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-05-05 05:29
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRPRR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 05:22:03
% EndTime: 2019-05-05 05:22:13
% DurationCPUTime: 5.68s
% Computational Cost: add. (81868->177), mult. (170831->240), div. (0->0), fcn. (124486->14), ass. (0->100)
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
t131 = -t101 * t83 + (t91 * t94 + t128) * t105;
t103 = cos(qJ(5));
t106 = qJD(3) ^ 2;
t104 = cos(qJ(3));
t123 = t104 * qJD(2);
t100 = sin(qJ(3));
t107 = qJD(2) ^ 2;
t125 = t101 * t94;
t120 = t101 * t128 + t105 * t83 + t91 * t125;
t51 = -t107 * pkin(2) + qJDD(2) * pkin(8) + t120;
t66 = -t94 * t82 + t97 * t91;
t126 = t100 * t66 + t104 * t51;
t78 = (-pkin(3) * t104 - qJ(4) * t100) * qJD(2);
t35 = -pkin(3) * t106 + qJDD(3) * qJ(4) + t78 * t123 + t126;
t122 = qJD(2) * qJD(3);
t117 = t104 * t122;
t50 = -qJDD(2) * pkin(2) - t107 * pkin(8) - t131;
t80 = t100 * qJDD(2) + t117;
t89 = t100 * t122;
t81 = t104 * qJDD(2) - t89;
t41 = (-t80 - t117) * qJ(4) + (-t81 + t89) * pkin(3) + t50;
t124 = qJD(2) * t100;
t92 = sin(pkin(12));
t95 = cos(pkin(12));
t74 = t92 * qJD(3) + t95 * t124;
t115 = -0.2e1 * qJD(4) * t74 - t92 * t35 + t95 * t41;
t102 = cos(qJ(6));
t64 = t92 * qJDD(3) + t95 * t80;
t73 = t95 * qJD(3) - t92 * t124;
t23 = (-t73 * t123 - t64) * pkin(9) + (t73 * t74 - t81) * pkin(4) + t115;
t121 = 0.2e1 * qJD(4) * t73 + t95 * t35 + t92 * t41;
t63 = t95 * qJDD(3) - t92 * t80;
t65 = -pkin(4) * t123 - t74 * pkin(9);
t72 = t73 ^ 2;
t25 = -pkin(4) * t72 + pkin(9) * t63 + t65 * t123 + t121;
t99 = sin(qJ(5));
t118 = t103 * t23 - t99 * t25;
t57 = t103 * t73 - t99 * t74;
t43 = t57 * qJD(5) + t103 * t64 + t99 * t63;
t58 = t103 * t74 + t99 * t73;
t75 = qJDD(5) - t81;
t88 = qJD(5) - t123;
t17 = (t57 * t88 - t43) * pkin(10) + (t57 * t58 + t75) * pkin(5) + t118;
t127 = t103 * t25 + t99 * t23;
t42 = -t58 * qJD(5) + t103 * t63 - t99 * t64;
t54 = t88 * pkin(5) - t58 * pkin(10);
t56 = t57 ^ 2;
t18 = -pkin(5) * t56 + pkin(10) * t42 - t54 * t88 + t127;
t98 = sin(qJ(6));
t45 = t102 * t57 - t98 * t58;
t28 = qJD(6) * t45 + t102 * t43 + t42 * t98;
t46 = t102 * t58 + t98 * t57;
t32 = -mrSges(7,1) * t45 + mrSges(7,2) * t46;
t87 = qJD(6) + t88;
t38 = -mrSges(7,2) * t87 + mrSges(7,3) * t45;
t71 = qJDD(6) + t75;
t14 = m(7) * (t102 * t17 - t18 * t98) - t28 * mrSges(7,3) + t71 * mrSges(7,1) - t46 * t32 + t87 * t38;
t27 = -qJD(6) * t46 + t102 * t42 - t43 * t98;
t39 = mrSges(7,1) * t87 - mrSges(7,3) * t46;
t15 = m(7) * (t102 * t18 + t17 * t98) + t27 * mrSges(7,3) - t71 * mrSges(7,2) + t45 * t32 - t87 * t39;
t47 = -t57 * mrSges(6,1) + t58 * mrSges(6,2);
t52 = -t88 * mrSges(6,2) + t57 * mrSges(6,3);
t12 = m(6) * t118 + t75 * mrSges(6,1) - t43 * mrSges(6,3) + t102 * t14 + t98 * t15 - t58 * t47 + t88 * t52;
t53 = t88 * mrSges(6,1) - t58 * mrSges(6,3);
t13 = m(6) * t127 - t75 * mrSges(6,2) + t42 * mrSges(6,3) + t102 * t15 - t98 * t14 + t57 * t47 - t88 * t53;
t59 = -t73 * mrSges(5,1) + t74 * mrSges(5,2);
t61 = mrSges(5,2) * t123 + t73 * mrSges(5,3);
t10 = m(5) * t115 - t81 * mrSges(5,1) - t64 * mrSges(5,3) + t103 * t12 - t61 * t123 + t99 * t13 - t74 * t59;
t62 = -mrSges(5,1) * t123 - t74 * mrSges(5,3);
t11 = m(5) * t121 + t81 * mrSges(5,2) + t63 * mrSges(5,3) + t103 * t13 - t99 * t12 + t62 * t123 + t73 * t59;
t84 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t124;
t85 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t123;
t130 = m(4) * t50 - t81 * mrSges(4,1) + t80 * mrSges(4,2) + t95 * t10 + t92 * t11 + (t100 * t84 - t104 * t85) * qJD(2);
t8 = m(3) * t131 + qJDD(2) * mrSges(3,1) - t107 * mrSges(3,2) - t130;
t129 = t105 * t8;
t116 = -t100 * t51 + t104 * t66;
t34 = -qJDD(3) * pkin(3) - qJ(4) * t106 + t78 * t124 + qJDD(4) - t116;
t109 = -pkin(4) * t63 - pkin(9) * t72 + t74 * t65 + t34;
t112 = t27 * mrSges(7,1) + t45 * t38 - m(7) * (-pkin(5) * t42 - pkin(10) * t56 + t54 * t58 + t109) - t28 * mrSges(7,2) - t46 * t39;
t110 = -m(6) * t109 + t42 * mrSges(6,1) - t43 * mrSges(6,2) + t57 * t52 - t58 * t53 + t112;
t108 = m(5) * t34 - t63 * mrSges(5,1) + t64 * mrSges(5,2) - t73 * t61 + t74 * t62 - t110;
t79 = (-mrSges(4,1) * t104 + mrSges(4,2) * t100) * qJD(2);
t16 = m(4) * t116 + qJDD(3) * mrSges(4,1) - t80 * mrSges(4,3) + qJD(3) * t85 - t79 * t124 - t108;
t9 = m(4) * t126 - qJDD(3) * mrSges(4,2) + t81 * mrSges(4,3) - qJD(3) * t84 - t92 * t10 + t95 * t11 + t79 * t123;
t4 = m(3) * t120 - t107 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t100 * t16 + t104 * t9;
t6 = m(3) * t66 + t100 * t9 + t104 * t16;
t119 = m(2) * t91 + t4 * t125 + t94 * t129 + t97 * t6;
t2 = m(2) * t83 - t101 * t8 + t105 * t4;
t1 = m(2) * t82 - t94 * t6 + (t101 * t4 + t129) * t97;
t3 = [-m(1) * g(1) - t1 * t93 + t2 * t96, t2, t4, t9, t11, t13, t15; -m(1) * g(2) + t1 * t96 + t2 * t93, t1, t8, t16, t10, t12, t14; -m(1) * g(3) + t119, t119, t6, t130, t108, -t110, -t112;];
f_new  = t3;
