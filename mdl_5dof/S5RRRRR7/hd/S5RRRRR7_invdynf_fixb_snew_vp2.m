% Calculate vector of cutting forces with Newton-Euler
% S5RRRRR7
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x6]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRRR7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR7_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR7_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR7_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR7_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR7_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:21:28
% EndTime: 2019-12-31 22:21:33
% DurationCPUTime: 2.33s
% Computational Cost: add. (27902->167), mult. (60486->222), div. (0->0), fcn. (43507->10), ass. (0->87)
t103 = qJD(1) * qJD(2);
t82 = sin(qJ(2));
t87 = cos(qJ(2));
t67 = t82 * qJDD(1) + t87 * t103;
t68 = t87 * qJDD(1) - t82 * t103;
t105 = qJD(1) * t82;
t69 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t105;
t104 = qJD(1) * t87;
t70 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t104;
t89 = qJD(1) ^ 2;
t81 = sin(qJ(3));
t86 = cos(qJ(3));
t62 = (t81 * t87 + t82 * t86) * qJD(1);
t45 = -t62 * qJD(3) - t81 * t67 + t86 * t68;
t61 = (-t81 * t82 + t86 * t87) * qJD(1);
t46 = t61 * qJD(3) + t86 * t67 + t81 * t68;
t77 = qJD(2) + qJD(3);
t56 = -t77 * mrSges(4,2) + t61 * mrSges(4,3);
t57 = t77 * mrSges(4,1) - t62 * mrSges(4,3);
t71 = qJD(2) * pkin(2) - pkin(7) * t105;
t78 = t87 ^ 2;
t83 = sin(qJ(1));
t88 = cos(qJ(1));
t102 = t83 * g(1) - t88 * g(2);
t96 = -qJDD(1) * pkin(1) - t102;
t93 = -t68 * pkin(2) + t71 * t105 + (-pkin(7) * t78 - pkin(6)) * t89 + t96;
t99 = -t88 * g(1) - t83 * g(2);
t64 = -t89 * pkin(1) + qJDD(1) * pkin(6) + t99;
t108 = t82 * t64;
t109 = pkin(2) * t89;
t39 = qJDD(2) * pkin(2) - t67 * pkin(7) - t108 + (pkin(7) * t103 + t82 * t109 - g(3)) * t87;
t101 = -t82 * g(3) + t87 * t64;
t40 = t68 * pkin(7) - qJD(2) * t71 - t78 * t109 + t101;
t100 = t86 * t39 - t81 * t40;
t76 = qJDD(2) + qJDD(3);
t19 = (t61 * t77 - t46) * pkin(8) + (t61 * t62 + t76) * pkin(3) + t100;
t106 = t81 * t39 + t86 * t40;
t58 = t77 * pkin(3) - t62 * pkin(8);
t60 = t61 ^ 2;
t21 = -t60 * pkin(3) + t45 * pkin(8) - t77 * t58 + t106;
t80 = sin(qJ(4));
t85 = cos(qJ(4));
t107 = t80 * t19 + t85 * t21;
t53 = t85 * t61 - t80 * t62;
t54 = t80 * t61 + t85 * t62;
t35 = -t53 * pkin(4) - t54 * pkin(9);
t74 = qJD(4) + t77;
t72 = t74 ^ 2;
t73 = qJDD(4) + t76;
t16 = -t72 * pkin(4) + t73 * pkin(9) + t53 * t35 + t107;
t28 = -t54 * qJD(4) + t85 * t45 - t80 * t46;
t29 = t53 * qJD(4) + t80 * t45 + t85 * t46;
t91 = -t45 * pkin(3) - t60 * pkin(8) + t62 * t58 + t93;
t17 = (-t53 * t74 - t29) * pkin(9) + (t54 * t74 - t28) * pkin(4) + t91;
t79 = sin(qJ(5));
t84 = cos(qJ(5));
t41 = -t79 * t54 + t84 * t74;
t23 = t41 * qJD(5) + t84 * t29 + t79 * t73;
t27 = qJDD(5) - t28;
t42 = t84 * t54 + t79 * t74;
t30 = -t41 * mrSges(6,1) + t42 * mrSges(6,2);
t50 = qJD(5) - t53;
t31 = -t50 * mrSges(6,2) + t41 * mrSges(6,3);
t13 = m(6) * (-t79 * t16 + t84 * t17) - t23 * mrSges(6,3) + t27 * mrSges(6,1) - t42 * t30 + t50 * t31;
t22 = -t42 * qJD(5) - t79 * t29 + t84 * t73;
t32 = t50 * mrSges(6,1) - t42 * mrSges(6,3);
t14 = m(6) * (t84 * t16 + t79 * t17) + t22 * mrSges(6,3) - t27 * mrSges(6,2) + t41 * t30 - t50 * t32;
t48 = -t74 * mrSges(5,2) + t53 * mrSges(5,3);
t49 = t74 * mrSges(5,1) - t54 * mrSges(5,3);
t95 = -m(5) * t91 + t28 * mrSges(5,1) - t29 * mrSges(5,2) - t84 * t13 - t79 * t14 + t53 * t48 - t54 * t49;
t92 = -m(4) * t93 + t45 * mrSges(4,1) - t46 * mrSges(4,2) + t61 * t56 - t62 * t57 + t95;
t111 = (t82 * t69 - t87 * t70) * qJD(1) + m(3) * (-t89 * pkin(6) + t96) - t68 * mrSges(3,1) + t67 * mrSges(3,2) - t92;
t34 = -t53 * mrSges(5,1) + t54 * mrSges(5,2);
t98 = t85 * t19 - t80 * t21;
t94 = m(6) * (-t73 * pkin(4) - t72 * pkin(9) + t54 * t35 - t98) - t22 * mrSges(6,1) + t23 * mrSges(6,2) - t41 * t31 + t42 * t32;
t10 = m(5) * t98 + t73 * mrSges(5,1) - t29 * mrSges(5,3) - t54 * t34 + t74 * t48 - t94;
t55 = -t61 * mrSges(4,1) + t62 * mrSges(4,2);
t9 = m(5) * t107 - t73 * mrSges(5,2) + t28 * mrSges(5,3) - t79 * t13 + t84 * t14 + t53 * t34 - t74 * t49;
t6 = m(4) * t100 + t76 * mrSges(4,1) - t46 * mrSges(4,3) + t85 * t10 - t62 * t55 + t77 * t56 + t80 * t9;
t66 = (-mrSges(3,1) * t87 + mrSges(3,2) * t82) * qJD(1);
t7 = m(4) * t106 - t76 * mrSges(4,2) + t45 * mrSges(4,3) - t80 * t10 + t61 * t55 - t77 * t57 + t85 * t9;
t4 = m(3) * (-t87 * g(3) - t108) - t67 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t66 * t105 + qJD(2) * t70 + t81 * t7 + t86 * t6;
t5 = m(3) * t101 - qJDD(2) * mrSges(3,2) + t68 * mrSges(3,3) - qJD(2) * t69 + t66 * t104 - t81 * t6 + t86 * t7;
t110 = t87 * t4 + t82 * t5;
t8 = m(2) * t102 + qJDD(1) * mrSges(2,1) - t89 * mrSges(2,2) - t111;
t1 = m(2) * t99 - t89 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t82 * t4 + t87 * t5;
t2 = [-m(1) * g(1) + t88 * t1 - t83 * t8, t1, t5, t7, t9, t14; -m(1) * g(2) + t83 * t1 + t88 * t8, t8, t4, t6, t10, t13; (-m(1) - m(2)) * g(3) + t110, -m(2) * g(3) + t110, t111, -t92, -t95, t94;];
f_new = t2;
