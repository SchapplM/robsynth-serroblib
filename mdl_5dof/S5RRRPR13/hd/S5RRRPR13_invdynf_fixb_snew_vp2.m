% Calculate vector of cutting forces with Newton-Euler
% S5RRRPR13
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRPR13_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR13_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR13_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR13_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR13_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR13_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR13_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:43:23
% EndTime: 2019-12-31 21:43:27
% DurationCPUTime: 1.59s
% Computational Cost: add. (17302->172), mult. (37479->227), div. (0->0), fcn. (27825->10), ass. (0->91)
t104 = qJD(1) * qJD(2);
t80 = sin(pkin(5));
t84 = sin(qJ(2));
t87 = cos(qJ(2));
t69 = (-qJDD(1) * t87 + t84 * t104) * t80;
t107 = qJD(1) * t80;
t103 = t84 * t107;
t118 = cos(qJ(3));
t81 = cos(pkin(5));
t77 = t81 * qJD(1) + qJD(2);
t83 = sin(qJ(3));
t57 = t83 * t103 - t118 * t77;
t106 = qJD(1) * t87;
t102 = t80 * t106;
t73 = -qJD(3) + t102;
t48 = t57 * mrSges(5,1) + t73 * mrSges(5,3);
t109 = t73 * mrSges(4,2) - t57 * mrSges(4,3) - t48;
t113 = mrSges(4,1) - mrSges(5,2);
t67 = (-pkin(2) * t87 - pkin(8) * t84) * t107;
t75 = t77 ^ 2;
t76 = t81 * qJDD(1) + qJDD(2);
t114 = t80 * t87;
t85 = sin(qJ(1));
t88 = cos(qJ(1));
t101 = t85 * g(1) - t88 * g(2);
t120 = pkin(7) * t80;
t89 = qJD(1) ^ 2;
t64 = qJDD(1) * pkin(1) + t89 * t120 + t101;
t116 = t64 * t81;
t100 = -t88 * g(1) - t85 * g(2);
t65 = -t89 * pkin(1) + qJDD(1) * t120 + t100;
t98 = -g(3) * t114 + t87 * t116 - t84 * t65;
t29 = -t76 * pkin(2) - t75 * pkin(8) + t67 * t103 - t98;
t58 = t118 * t103 + t83 * t77;
t68 = (qJDD(1) * t84 + t87 * t104) * t80;
t39 = t58 * qJD(3) - t118 * t76 + t83 * t68;
t40 = -t57 * qJD(3) + t118 * t68 + t83 * t76;
t47 = -t73 * mrSges(4,1) - t58 * mrSges(4,3);
t117 = t57 * t73;
t41 = t57 * pkin(3) - t58 * qJ(4);
t61 = qJDD(3) + t69;
t72 = t73 ^ 2;
t108 = t84 * t116 + t87 * t65;
t30 = -t75 * pkin(2) + t76 * pkin(8) + (-g(3) * t84 + t67 * t106) * t80 + t108;
t119 = t81 * g(3);
t31 = t69 * pkin(2) - t68 * pkin(8) - t119 + (-t64 + (pkin(2) * t84 - pkin(8) * t87) * t77 * qJD(1)) * t80;
t99 = t118 * t31 - t83 * t30;
t19 = -t61 * pkin(3) - t72 * qJ(4) + t58 * t41 + qJDD(4) - t99;
t14 = (t57 * t58 - t61) * pkin(9) + (t40 - t117) * pkin(4) + t19;
t50 = t58 * pkin(4) + t73 * pkin(9);
t56 = t57 ^ 2;
t121 = -2 * qJD(4);
t90 = (-t40 - t117) * qJ(4) + t29 + (-t73 * pkin(3) + t121) * t58;
t17 = -t56 * pkin(4) - t58 * t50 + (pkin(3) + pkin(9)) * t39 + t90;
t82 = sin(qJ(5));
t86 = cos(qJ(5));
t44 = t86 * t57 + t82 * t73;
t24 = t44 * qJD(5) + t82 * t39 + t86 * t61;
t45 = t82 * t57 - t86 * t73;
t32 = -t44 * mrSges(6,1) + t45 * mrSges(6,2);
t55 = qJD(5) + t58;
t34 = -t55 * mrSges(6,2) + t44 * mrSges(6,3);
t37 = qJDD(5) + t40;
t12 = m(6) * (t86 * t14 - t82 * t17) - t24 * mrSges(6,3) + t37 * mrSges(6,1) - t45 * t32 + t55 * t34;
t23 = -t45 * qJD(5) + t86 * t39 - t82 * t61;
t35 = t55 * mrSges(6,1) - t45 * mrSges(6,3);
t13 = m(6) * (t82 * t14 + t86 * t17) + t23 * mrSges(6,3) - t37 * mrSges(6,2) + t44 * t32 - t55 * t35;
t49 = t58 * mrSges(5,1) - t73 * mrSges(5,2);
t97 = t82 * t12 - t86 * t13 - m(5) * (t39 * pkin(3) + t90) + t40 * mrSges(5,3) + t58 * t49;
t122 = m(4) * t29 + t40 * mrSges(4,2) + t109 * t57 + t113 * t39 + t58 * t47 - t97;
t115 = t80 * t84;
t112 = -mrSges(4,3) - mrSges(5,1);
t111 = t118 * t30 + t83 * t31;
t43 = -t57 * mrSges(5,2) - t58 * mrSges(5,3);
t110 = -t57 * mrSges(4,1) - t58 * mrSges(4,2) - t43;
t93 = -t72 * pkin(3) + t61 * qJ(4) - t57 * t41 + t111;
t94 = -t23 * mrSges(6,1) - t44 * t34 + m(6) * (-t39 * pkin(4) - t56 * pkin(9) + (t121 - t50) * t73 + t93) + t24 * mrSges(6,2) + t45 * t35;
t92 = -m(5) * (0.2e1 * qJD(4) * t73 - t93) + t94;
t10 = m(4) * t111 + (t47 - t49) * t73 + (-mrSges(4,2) + mrSges(5,3)) * t61 + t110 * t57 + t112 * t39 + t92;
t62 = t77 * mrSges(3,1) - mrSges(3,3) * t103;
t66 = (-mrSges(3,1) * t87 + mrSges(3,2) * t84) * t107;
t95 = -m(5) * t19 - t86 * t12 - t82 * t13;
t9 = m(4) * t99 - t109 * t73 + t110 * t58 + t112 * t40 + t113 * t61 + t95;
t4 = m(3) * (-g(3) * t115 + t108) - t69 * mrSges(3,3) - t76 * mrSges(3,2) + t66 * t102 - t77 * t62 + t118 * t10 - t83 * t9;
t63 = -t77 * mrSges(3,2) + mrSges(3,3) * t102;
t6 = m(3) * (-t80 * t64 - t119) + t68 * mrSges(3,2) + t69 * mrSges(3,1) + t83 * t10 + t118 * t9 + (t62 * t84 - t63 * t87) * t107;
t8 = m(3) * t98 + t76 * mrSges(3,1) - t68 * mrSges(3,3) - t66 * t103 + t77 * t63 - t122;
t105 = t8 * t114 + t4 * t115 + t81 * t6;
t2 = m(2) * t100 - t89 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t87 * t4 - t84 * t8;
t1 = m(2) * t101 + qJDD(1) * mrSges(2,1) - t89 * mrSges(2,2) - t80 * t6 + (t4 * t84 + t8 * t87) * t81;
t3 = [-m(1) * g(1) - t1 * t85 + t2 * t88, t2, t4, t10, -t39 * mrSges(5,2) - t57 * t48 - t97, t13; -m(1) * g(2) + t1 * t88 + t2 * t85, t1, t8, t9, t39 * mrSges(5,1) - t61 * mrSges(5,3) + t57 * t43 + t73 * t49 - t92, t12; (-m(1) - m(2)) * g(3) + t105, -m(2) * g(3) + t105, t6, t122, t40 * mrSges(5,1) + t61 * mrSges(5,2) + t58 * t43 - t73 * t48 - t95, t94;];
f_new = t3;
