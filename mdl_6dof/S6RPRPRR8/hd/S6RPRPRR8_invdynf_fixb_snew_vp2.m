% Calculate vector of cutting forces with Newton-Euler
% S6RPRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-05-05 19:23
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRPRR8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR8_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR8_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR8_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR8_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR8_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 19:19:44
% EndTime: 2019-05-05 19:19:50
% DurationCPUTime: 2.38s
% Computational Cost: add. (30519->179), mult. (66932->231), div. (0->0), fcn. (45045->10), ass. (0->92)
t101 = cos(qJ(1));
t97 = sin(qJ(1));
t114 = -t101 * g(1) - t97 * g(2);
t112 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t114;
t100 = cos(qJ(3));
t103 = qJD(1) ^ 2;
t123 = qJD(1) * qJD(3);
t96 = sin(qJ(3));
t118 = t96 * t123;
t119 = t97 * g(1) - t101 * g(2);
t110 = -t103 * qJ(2) + qJDD(2) - t119;
t131 = -pkin(1) - pkin(7);
t64 = t131 * qJDD(1) + t110;
t127 = t96 * g(3) + t100 * t64;
t80 = t100 * qJDD(1) - t118;
t36 = (-t80 - t118) * qJ(4) + (-t100 * t103 * t96 + qJDD(3)) * pkin(3) + t127;
t117 = -t100 * g(3) + t96 * t64;
t79 = -t96 * qJDD(1) - t100 * t123;
t124 = qJD(1) * t100;
t82 = qJD(3) * pkin(3) - qJ(4) * t124;
t91 = t96 ^ 2;
t37 = -t91 * t103 * pkin(3) + t79 * qJ(4) - qJD(3) * t82 + t117;
t126 = qJD(1) * t96;
t92 = sin(pkin(10));
t93 = cos(pkin(10));
t72 = t93 * t124 - t92 * t126;
t133 = -0.2e1 * qJD(4) * t72 + t93 * t36 - t92 * t37;
t132 = -m(2) - m(3);
t130 = (mrSges(2,1) - mrSges(3,2));
t129 = -mrSges(2,2) + mrSges(3,3);
t102 = qJD(3) ^ 2;
t71 = -t92 * t124 - t93 * t126;
t121 = 0.2e1 * qJD(4) * t71 + t92 * t36 + t93 * t37;
t49 = -t71 * pkin(4) - t72 * pkin(8);
t20 = -t102 * pkin(4) + qJDD(3) * pkin(8) + t71 * t49 + t121;
t106 = -t79 * pkin(3) + qJDD(4) + t82 * t124 + (-qJ(4) * t91 + t131) * t103 + t112;
t53 = t93 * t79 - t92 * t80;
t54 = t92 * t79 + t93 * t80;
t26 = (-qJD(3) * t71 - t54) * pkin(8) + (qJD(3) * t72 - t53) * pkin(4) + t106;
t95 = sin(qJ(5));
t99 = cos(qJ(5));
t128 = t99 * t20 + t95 * t26;
t69 = qJD(5) - t71;
t19 = -qJDD(3) * pkin(4) - t102 * pkin(8) + t72 * t49 - t133;
t57 = t95 * qJD(3) + t99 * t72;
t31 = -t57 * qJD(5) + t99 * qJDD(3) - t95 * t54;
t56 = t99 * qJD(3) - t95 * t72;
t32 = t56 * qJD(5) + t95 * qJDD(3) + t99 * t54;
t94 = sin(qJ(6));
t98 = cos(qJ(6));
t41 = t94 * t56 + t98 * t57;
t22 = -t41 * qJD(6) + t98 * t31 - t94 * t32;
t40 = t98 * t56 - t94 * t57;
t23 = t40 * qJD(6) + t94 * t31 + t98 * t32;
t66 = qJD(6) + t69;
t29 = -t66 * mrSges(7,2) + t40 * mrSges(7,3);
t30 = t66 * mrSges(7,1) - t41 * mrSges(7,3);
t46 = t69 * pkin(5) - t57 * pkin(9);
t55 = t56 ^ 2;
t109 = t22 * mrSges(7,1) + t40 * t29 - m(7) * (-t31 * pkin(5) - t55 * pkin(9) + t57 * t46 + t19) - t23 * mrSges(7,2) - t41 * t30;
t44 = -t69 * mrSges(6,2) + t56 * mrSges(6,3);
t45 = t69 * mrSges(6,1) - t57 * mrSges(6,3);
t104 = m(6) * t19 - t31 * mrSges(6,1) + t32 * mrSges(6,2) - t56 * t44 + t57 * t45 - t109;
t48 = -t71 * mrSges(5,1) + t72 * mrSges(5,2);
t62 = -qJD(3) * mrSges(5,2) + t71 * mrSges(5,3);
t11 = m(5) * t133 + qJDD(3) * mrSges(5,1) - t54 * mrSges(5,3) + qJD(3) * t62 - t72 * t48 - t104;
t116 = -t95 * t20 + t99 * t26;
t52 = qJDD(5) - t53;
t14 = (t56 * t69 - t32) * pkin(9) + (t56 * t57 + t52) * pkin(5) + t116;
t15 = -t55 * pkin(5) + t31 * pkin(9) - t69 * t46 + t128;
t28 = -t40 * mrSges(7,1) + t41 * mrSges(7,2);
t50 = qJDD(6) + t52;
t12 = m(7) * (t98 * t14 - t94 * t15) - t23 * mrSges(7,3) + t50 * mrSges(7,1) - t41 * t28 + t66 * t29;
t13 = m(7) * (t94 * t14 + t98 * t15) + t22 * mrSges(7,3) - t50 * mrSges(7,2) + t40 * t28 - t66 * t30;
t42 = -t56 * mrSges(6,1) + t57 * mrSges(6,2);
t10 = m(6) * t128 - t52 * mrSges(6,2) + t31 * mrSges(6,3) - t94 * t12 + t98 * t13 + t56 * t42 - t69 * t45;
t63 = qJD(3) * mrSges(5,1) - t72 * mrSges(5,3);
t9 = m(6) * t116 + t52 * mrSges(6,1) - t32 * mrSges(6,3) + t98 * t12 + t94 * t13 - t57 * t42 + t69 * t44;
t6 = m(5) * t121 - qJDD(3) * mrSges(5,2) + t53 * mrSges(5,3) - qJD(3) * t63 + t99 * t10 + t71 * t48 - t95 * t9;
t78 = (mrSges(4,1) * t96 + mrSges(4,2) * t100) * qJD(1);
t81 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t126;
t3 = m(4) * t127 + qJDD(3) * mrSges(4,1) - t80 * mrSges(4,3) + qJD(3) * t81 + t93 * t11 - t78 * t124 + t92 * t6;
t83 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t124;
t4 = m(4) * t117 - qJDD(3) * mrSges(4,2) + t79 * mrSges(4,3) - qJD(3) * t83 - t92 * t11 - t78 * t126 + t93 * t6;
t120 = t100 * t4 - t96 * t3;
t111 = -m(3) * (-qJDD(1) * pkin(1) + t110) - t100 * t3 - t96 * t4;
t108 = m(5) * t106 - t53 * mrSges(5,1) + t54 * mrSges(5,2) + t95 * t10 - t71 * t62 + t72 * t63 + t99 * t9;
t107 = -t79 * mrSges(4,1) + t108 + m(4) * (t131 * t103 + t112) + t81 * t126 + t83 * t124 + t80 * mrSges(4,2);
t105 = -m(3) * (t103 * pkin(1) - t112) + t107;
t5 = m(2) * t114 + t129 * qJDD(1) - (t130 * t103) + t105;
t1 = m(2) * t119 + t130 * qJDD(1) + t129 * t103 + t111;
t2 = [-m(1) * g(1) - t97 * t1 + t101 * t5, t5, -m(3) * g(3) + t120, t4, t6, t10, t13; -m(1) * g(2) + t101 * t1 + t97 * t5, t1, -(t103 * mrSges(3,2)) - qJDD(1) * mrSges(3,3) - t105, t3, t11, t9, t12; (-m(1) + t132) * g(3) + t120, t132 * g(3) + t120, qJDD(1) * mrSges(3,2) - t103 * mrSges(3,3) - t111, t107, t108, t104, -t109;];
f_new  = t2;
