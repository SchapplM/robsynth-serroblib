% Calculate vector of cutting forces with Newton-Euler
% S5RRRPR12
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRPR12_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR12_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR12_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR12_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR12_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR12_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR12_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:37:17
% EndTime: 2019-12-31 21:37:24
% DurationCPUTime: 3.71s
% Computational Cost: add. (47725->173), mult. (104744->244), div. (0->0), fcn. (81360->12), ass. (0->95)
t108 = qJD(1) * qJD(2);
t85 = sin(pkin(5));
t90 = sin(qJ(2));
t93 = cos(qJ(2));
t74 = (-qJDD(1) * t93 + t90 * t108) * t85;
t119 = pkin(7) * t85;
t87 = cos(pkin(5));
t118 = t87 * g(3);
t117 = cos(qJ(3));
t91 = sin(qJ(1));
t94 = cos(qJ(1));
t104 = t91 * g(1) - t94 * g(2);
t95 = qJD(1) ^ 2;
t69 = qJDD(1) * pkin(1) + t95 * t119 + t104;
t116 = t69 * t87;
t115 = t85 * t90;
t114 = t85 * t93;
t110 = qJD(1) * t93;
t102 = -t94 * g(1) - t91 * g(2);
t70 = -t95 * pkin(1) + qJDD(1) * t119 + t102;
t112 = t90 * t116 + t93 * t70;
t111 = qJD(1) * t85;
t72 = (-pkin(2) * t93 - pkin(8) * t90) * t111;
t81 = t87 * qJD(1) + qJD(2);
t79 = t81 ^ 2;
t80 = t87 * qJDD(1) + qJDD(2);
t36 = -t79 * pkin(2) + t80 * pkin(8) + (-g(3) * t90 + t72 * t110) * t85 + t112;
t73 = (qJDD(1) * t90 + t93 * t108) * t85;
t37 = t74 * pkin(2) - t73 * pkin(8) - t118 + (-t69 + (pkin(2) * t90 - pkin(8) * t93) * t81 * qJD(1)) * t85;
t89 = sin(qJ(3));
t113 = t117 * t36 + t89 * t37;
t105 = t85 * t110;
t101 = t117 * t37 - t89 * t36;
t106 = t90 * t111;
t62 = t89 * t106 - t117 * t81;
t49 = -t62 * qJD(3) + t117 * t73 + t89 * t80;
t63 = t117 * t106 + t89 * t81;
t51 = t62 * mrSges(4,1) + t63 * mrSges(4,2);
t78 = qJD(3) - t105;
t57 = -t78 * mrSges(4,2) - t62 * mrSges(4,3);
t66 = qJDD(3) + t74;
t50 = t62 * pkin(3) - t63 * qJ(4);
t77 = t78 ^ 2;
t20 = -t66 * pkin(3) - t77 * qJ(4) + t63 * t50 + qJDD(4) - t101;
t84 = sin(pkin(10));
t86 = cos(pkin(10));
t40 = -t84 * t49 + t86 * t66;
t41 = t86 * t49 + t84 * t66;
t55 = -t84 * t63 + t86 * t78;
t44 = -t62 * mrSges(5,2) + t55 * mrSges(5,3);
t56 = t86 * t63 + t84 * t78;
t45 = t62 * mrSges(5,1) - t56 * mrSges(5,3);
t88 = sin(qJ(5));
t92 = cos(qJ(5));
t39 = t88 * t55 + t92 * t56;
t26 = -t39 * qJD(5) + t92 * t40 - t88 * t41;
t38 = t92 * t55 - t88 * t56;
t27 = t38 * qJD(5) + t88 * t40 + t92 * t41;
t61 = qJD(5) + t62;
t30 = -t61 * mrSges(6,2) + t38 * mrSges(6,3);
t31 = t61 * mrSges(6,1) - t39 * mrSges(6,3);
t46 = t62 * pkin(4) - t56 * pkin(9);
t54 = t55 ^ 2;
t98 = t26 * mrSges(6,1) + t38 * t30 - m(6) * (-t40 * pkin(4) - t54 * pkin(9) + t56 * t46 + t20) - t27 * mrSges(6,2) - t39 * t31;
t96 = m(5) * t20 - t40 * mrSges(5,1) + t41 * mrSges(5,2) - t55 * t44 + t56 * t45 - t98;
t12 = m(4) * t101 + t66 * mrSges(4,1) - t49 * mrSges(4,3) - t63 * t51 + t78 * t57 - t96;
t67 = t81 * mrSges(3,1) - mrSges(3,3) * t106;
t71 = (-mrSges(3,1) * t93 + mrSges(3,2) * t90) * t111;
t21 = -t77 * pkin(3) + t66 * qJ(4) - t62 * t50 + t113;
t100 = -g(3) * t114 + t93 * t116 - t90 * t70;
t35 = -t80 * pkin(2) - t79 * pkin(8) + t72 * t106 - t100;
t48 = t63 * qJD(3) - t117 * t80 + t89 * t73;
t24 = (t62 * t78 - t49) * qJ(4) + (t63 * t78 + t48) * pkin(3) + t35;
t103 = -0.2e1 * qJD(4) * t56 - t84 * t21 + t86 * t24;
t15 = (t55 * t62 - t41) * pkin(9) + (t55 * t56 + t48) * pkin(4) + t103;
t107 = 0.2e1 * qJD(4) * t55 + t86 * t21 + t84 * t24;
t16 = -t54 * pkin(4) + t40 * pkin(9) - t62 * t46 + t107;
t29 = -t38 * mrSges(6,1) + t39 * mrSges(6,2);
t47 = qJDD(5) + t48;
t13 = m(6) * (t92 * t15 - t88 * t16) - t27 * mrSges(6,3) + t47 * mrSges(6,1) - t39 * t29 + t61 * t30;
t14 = m(6) * (t88 * t15 + t92 * t16) + t26 * mrSges(6,3) - t47 * mrSges(6,2) + t38 * t29 - t61 * t31;
t42 = -t55 * mrSges(5,1) + t56 * mrSges(5,2);
t10 = m(5) * t103 + t48 * mrSges(5,1) - t41 * mrSges(5,3) + t92 * t13 + t88 * t14 - t56 * t42 + t62 * t44;
t11 = m(5) * t107 - t48 * mrSges(5,2) + t40 * mrSges(5,3) - t88 * t13 + t92 * t14 + t55 * t42 - t62 * t45;
t58 = t78 * mrSges(4,1) - t63 * mrSges(4,3);
t9 = m(4) * t113 - t66 * mrSges(4,2) - t48 * mrSges(4,3) - t84 * t10 + t86 * t11 - t62 * t51 - t78 * t58;
t4 = m(3) * (-g(3) * t115 + t112) - t74 * mrSges(3,3) - t80 * mrSges(3,2) + t71 * t105 - t81 * t67 + t117 * t9 - t89 * t12;
t68 = -t81 * mrSges(3,2) + mrSges(3,3) * t105;
t6 = m(3) * (-t85 * t69 - t118) + t73 * mrSges(3,2) + t74 * mrSges(3,1) + t89 * t9 + t117 * t12 + (t67 * t90 - t68 * t93) * t111;
t97 = m(4) * t35 + t48 * mrSges(4,1) + t49 * mrSges(4,2) + t86 * t10 + t84 * t11 + t62 * t57 + t63 * t58;
t8 = m(3) * t100 + t80 * mrSges(3,1) - t73 * mrSges(3,3) - t71 * t106 + t81 * t68 - t97;
t109 = t8 * t114 + t4 * t115 + t87 * t6;
t2 = m(2) * t102 - t95 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t93 * t4 - t90 * t8;
t1 = m(2) * t104 + qJDD(1) * mrSges(2,1) - t95 * mrSges(2,2) - t85 * t6 + (t90 * t4 + t93 * t8) * t87;
t3 = [-m(1) * g(1) - t91 * t1 + t94 * t2, t2, t4, t9, t11, t14; -m(1) * g(2) + t94 * t1 + t91 * t2, t1, t8, t12, t10, t13; (-m(1) - m(2)) * g(3) + t109, -m(2) * g(3) + t109, t6, t97, t96, -t98;];
f_new = t3;
