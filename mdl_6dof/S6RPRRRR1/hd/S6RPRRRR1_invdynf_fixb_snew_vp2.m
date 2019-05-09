% Calculate vector of cutting forces with Newton-Euler
% S6RPRRRR1
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
% Datum: 2019-05-06 02:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRRR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 02:31:32
% EndTime: 2019-05-06 02:31:40
% DurationCPUTime: 3.88s
% Computational Cost: add. (54339->180), mult. (113078->238), div. (0->0), fcn. (79909->12), ass. (0->96)
t100 = cos(qJ(3));
t102 = qJD(1) ^ 2;
t101 = cos(qJ(1));
t96 = sin(qJ(1));
t117 = t96 * g(1) - t101 * g(2);
t71 = qJDD(1) * pkin(1) + t117;
t112 = -t101 * g(1) - t96 * g(2);
t73 = -t102 * pkin(1) + t112;
t90 = sin(pkin(11));
t91 = cos(pkin(11));
t114 = t91 * t71 - t90 * t73;
t109 = -qJDD(1) * pkin(2) - t114;
t95 = sin(qJ(3));
t121 = qJD(1) * t95;
t119 = qJD(1) * qJD(3);
t75 = t100 * qJDD(1) - t95 * t119;
t78 = qJD(3) * pkin(3) - pkin(8) * t121;
t88 = t100 ^ 2;
t106 = -t75 * pkin(3) + t78 * t121 + (-pkin(8) * t88 - pkin(7)) * t102 + t109;
t94 = sin(qJ(4));
t99 = cos(qJ(4));
t69 = (t100 * t94 + t95 * t99) * qJD(1);
t113 = t100 * t119;
t74 = t95 * qJDD(1) + t113;
t48 = -t69 * qJD(4) - t94 * t74 + t99 * t75;
t87 = qJD(3) + qJD(4);
t63 = t87 * pkin(4) - t69 * pkin(9);
t68 = (t100 * t99 - t94 * t95) * qJD(1);
t64 = t68 ^ 2;
t104 = -t48 * pkin(4) - t64 * pkin(9) + t69 * t63 + t106;
t122 = t90 * t71 + t91 * t73;
t57 = -t102 * pkin(2) + qJDD(1) * pkin(7) + t122;
t89 = -g(3) + qJDD(2);
t115 = t100 * t89 - t95 * t57;
t39 = (-t74 + t113) * pkin(8) + (t100 * t102 * t95 + qJDD(3)) * pkin(3) + t115;
t123 = t100 * t57 + t95 * t89;
t42 = -t88 * t102 * pkin(3) + t75 * pkin(8) - qJD(3) * t78 + t123;
t116 = t99 * t39 - t94 * t42;
t49 = t68 * qJD(4) + t99 * t74 + t94 * t75;
t86 = qJDD(3) + qJDD(4);
t21 = (t68 * t87 - t49) * pkin(9) + (t68 * t69 + t86) * pkin(4) + t116;
t124 = t94 * t39 + t99 * t42;
t23 = -t64 * pkin(4) + t48 * pkin(9) - t87 * t63 + t124;
t93 = sin(qJ(5));
t98 = cos(qJ(5));
t125 = t93 * t21 + t98 * t23;
t58 = t98 * t68 - t93 * t69;
t59 = t93 * t68 + t98 * t69;
t41 = -t58 * pkin(5) - t59 * pkin(10);
t83 = qJD(5) + t87;
t81 = t83 ^ 2;
t82 = qJDD(5) + t86;
t18 = -t81 * pkin(5) + t82 * pkin(10) + t58 * t41 + t125;
t30 = -t59 * qJD(5) + t98 * t48 - t93 * t49;
t31 = t58 * qJD(5) + t93 * t48 + t98 * t49;
t19 = (-t58 * t83 - t31) * pkin(10) + (t59 * t83 - t30) * pkin(5) + t104;
t92 = sin(qJ(6));
t97 = cos(qJ(6));
t45 = -t92 * t59 + t97 * t83;
t25 = t45 * qJD(6) + t97 * t31 + t92 * t82;
t29 = qJDD(6) - t30;
t46 = t97 * t59 + t92 * t83;
t32 = -t45 * mrSges(7,1) + t46 * mrSges(7,2);
t53 = qJD(6) - t58;
t33 = -t53 * mrSges(7,2) + t45 * mrSges(7,3);
t15 = m(7) * (-t92 * t18 + t97 * t19) - t25 * mrSges(7,3) + t29 * mrSges(7,1) - t46 * t32 + t53 * t33;
t24 = -t46 * qJD(6) - t92 * t31 + t97 * t82;
t34 = t53 * mrSges(7,1) - t46 * mrSges(7,3);
t16 = m(7) * (t97 * t18 + t92 * t19) + t24 * mrSges(7,3) - t29 * mrSges(7,2) + t45 * t32 - t53 * t34;
t50 = -t83 * mrSges(6,2) + t58 * mrSges(6,3);
t51 = t83 * mrSges(6,1) - t59 * mrSges(6,3);
t108 = -m(6) * t104 + t30 * mrSges(6,1) - t31 * mrSges(6,2) - t97 * t15 - t92 * t16 + t58 * t50 - t59 * t51;
t61 = -t87 * mrSges(5,2) + t68 * mrSges(5,3);
t62 = t87 * mrSges(5,1) - t69 * mrSges(5,3);
t105 = -m(5) * t106 + t48 * mrSges(5,1) - t49 * mrSges(5,2) + t68 * t61 - t69 * t62 + t108;
t76 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t121;
t120 = qJD(1) * t100;
t77 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t120;
t126 = -(t100 * t77 - t95 * t76) * qJD(1) + m(4) * (-t102 * pkin(7) + t109) - t75 * mrSges(4,1) + t74 * mrSges(4,2) - t105;
t72 = (-mrSges(4,1) * t100 + mrSges(4,2) * t95) * qJD(1);
t40 = -t58 * mrSges(6,1) + t59 * mrSges(6,2);
t11 = m(6) * t125 - t82 * mrSges(6,2) + t30 * mrSges(6,3) - t92 * t15 + t97 * t16 + t58 * t40 - t83 * t51;
t111 = t98 * t21 - t93 * t23;
t107 = m(7) * (-t82 * pkin(5) - t81 * pkin(10) + t59 * t41 - t111) - t24 * mrSges(7,1) + t25 * mrSges(7,2) - t45 * t33 + t46 * t34;
t12 = m(6) * t111 + t82 * mrSges(6,1) - t31 * mrSges(6,3) - t59 * t40 + t83 * t50 - t107;
t60 = -t68 * mrSges(5,1) + t69 * mrSges(5,2);
t8 = m(5) * t116 + t86 * mrSges(5,1) - t49 * mrSges(5,3) + t93 * t11 + t98 * t12 - t69 * t60 + t87 * t61;
t9 = m(5) * t124 - t86 * mrSges(5,2) + t48 * mrSges(5,3) + t98 * t11 - t93 * t12 + t68 * t60 - t87 * t62;
t6 = m(4) * t115 + qJDD(3) * mrSges(4,1) - t74 * mrSges(4,3) + qJD(3) * t77 - t72 * t121 + t99 * t8 + t94 * t9;
t7 = m(4) * t123 - qJDD(3) * mrSges(4,2) + t75 * mrSges(4,3) - qJD(3) * t76 + t72 * t120 - t94 * t8 + t99 * t9;
t118 = m(3) * t89 + t100 * t6 + t95 * t7;
t10 = m(3) * t114 + qJDD(1) * mrSges(3,1) - t102 * mrSges(3,2) - t126;
t3 = m(3) * t122 - t102 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t100 * t7 - t95 * t6;
t2 = m(2) * t112 - t102 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t90 * t10 + t91 * t3;
t1 = m(2) * t117 + qJDD(1) * mrSges(2,1) - t102 * mrSges(2,2) + t91 * t10 + t90 * t3;
t4 = [-m(1) * g(1) - t96 * t1 + t101 * t2, t2, t3, t7, t9, t11, t16; -m(1) * g(2) + t101 * t1 + t96 * t2, t1, t10, t6, t8, t12, t15; (-m(1) - m(2)) * g(3) + t118, -m(2) * g(3) + t118, t118, t126, -t105, -t108, t107;];
f_new  = t4;
