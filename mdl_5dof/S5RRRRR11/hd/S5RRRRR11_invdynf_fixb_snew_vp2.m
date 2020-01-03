% Calculate vector of cutting forces with Newton-Euler
% S5RRRRR11
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRRR11_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR11_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR11_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR11_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR11_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR11_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR11_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:39:00
% EndTime: 2019-12-31 22:39:09
% DurationCPUTime: 3.80s
% Computational Cost: add. (50430->174), mult. (107910->242), div. (0->0), fcn. (84456->12), ass. (0->97)
t109 = qJD(1) * qJD(2);
t85 = sin(pkin(5));
t90 = sin(qJ(2));
t95 = cos(qJ(2));
t76 = (-qJDD(1) * t95 + t90 * t109) * t85;
t120 = pkin(7) * t85;
t86 = cos(pkin(5));
t119 = t86 * g(3);
t91 = sin(qJ(1));
t96 = cos(qJ(1));
t106 = t91 * g(1) - t96 * g(2);
t97 = qJD(1) ^ 2;
t71 = qJDD(1) * pkin(1) + t97 * t120 + t106;
t118 = t71 * t86;
t117 = t85 * t90;
t116 = t85 * t95;
t111 = qJD(1) * t95;
t103 = -t96 * g(1) - t91 * g(2);
t72 = -t97 * pkin(1) + qJDD(1) * t120 + t103;
t113 = t90 * t118 + t95 * t72;
t112 = qJD(1) * t85;
t74 = (-pkin(2) * t95 - pkin(8) * t90) * t112;
t82 = t86 * qJD(1) + qJD(2);
t80 = t82 ^ 2;
t81 = t86 * qJDD(1) + qJDD(2);
t38 = -t80 * pkin(2) + t81 * pkin(8) + (-g(3) * t90 + t74 * t111) * t85 + t113;
t75 = (qJDD(1) * t90 + t95 * t109) * t85;
t39 = t76 * pkin(2) - t75 * pkin(8) - t119 + (-t71 + (pkin(2) * t90 - pkin(8) * t95) * t82 * qJD(1)) * t85;
t89 = sin(qJ(3));
t94 = cos(qJ(3));
t114 = t94 * t38 + t89 * t39;
t108 = t90 * t112;
t63 = -t89 * t108 + t94 * t82;
t64 = t94 * t108 + t89 * t82;
t52 = -t63 * pkin(3) - t64 * pkin(9);
t68 = qJDD(3) + t76;
t107 = t85 * t111;
t79 = qJD(3) - t107;
t77 = t79 ^ 2;
t24 = -t77 * pkin(3) + t68 * pkin(9) + t63 * t52 + t114;
t102 = -g(3) * t116 + t95 * t118 - t90 * t72;
t37 = -t81 * pkin(2) - t80 * pkin(8) + t74 * t108 - t102;
t49 = -t64 * qJD(3) - t89 * t75 + t94 * t81;
t50 = t63 * qJD(3) + t94 * t75 + t89 * t81;
t27 = (-t63 * t79 - t50) * pkin(9) + (t64 * t79 - t49) * pkin(3) + t37;
t88 = sin(qJ(4));
t93 = cos(qJ(4));
t115 = t93 * t24 + t88 * t27;
t104 = -t89 * t38 + t94 * t39;
t51 = -t63 * mrSges(4,1) + t64 * mrSges(4,2);
t56 = -t79 * mrSges(4,2) + t63 * mrSges(4,3);
t55 = t93 * t64 + t88 * t79;
t30 = -t55 * qJD(4) - t88 * t50 + t93 * t68;
t54 = -t88 * t64 + t93 * t79;
t31 = t54 * qJD(4) + t93 * t50 + t88 * t68;
t87 = sin(qJ(5));
t92 = cos(qJ(5));
t41 = t87 * t54 + t92 * t55;
t20 = -t41 * qJD(5) + t92 * t30 - t87 * t31;
t40 = t92 * t54 - t87 * t55;
t21 = t40 * qJD(5) + t87 * t30 + t92 * t31;
t23 = -t68 * pkin(3) - t77 * pkin(9) + t64 * t52 - t104;
t62 = qJD(4) - t63;
t60 = qJD(5) + t62;
t32 = -t60 * mrSges(6,2) + t40 * mrSges(6,3);
t33 = t60 * mrSges(6,1) - t41 * mrSges(6,3);
t46 = t62 * pkin(4) - t55 * pkin(10);
t53 = t54 ^ 2;
t100 = t20 * mrSges(6,1) + t40 * t32 - m(6) * (-t30 * pkin(4) - t53 * pkin(10) + t55 * t46 + t23) - t21 * mrSges(6,2) - t41 * t33;
t44 = -t62 * mrSges(5,2) + t54 * mrSges(5,3);
t45 = t62 * mrSges(5,1) - t55 * mrSges(5,3);
t98 = m(5) * t23 - t30 * mrSges(5,1) + t31 * mrSges(5,2) - t54 * t44 + t55 * t45 - t100;
t12 = m(4) * t104 + t68 * mrSges(4,1) - t50 * mrSges(4,3) - t64 * t51 + t79 * t56 - t98;
t69 = t82 * mrSges(3,1) - mrSges(3,3) * t108;
t73 = (-mrSges(3,1) * t95 + mrSges(3,2) * t90) * t112;
t105 = -t88 * t24 + t93 * t27;
t48 = qJDD(4) - t49;
t15 = (t54 * t62 - t31) * pkin(10) + (t54 * t55 + t48) * pkin(4) + t105;
t16 = -t53 * pkin(4) + t30 * pkin(10) - t62 * t46 + t115;
t29 = -t40 * mrSges(6,1) + t41 * mrSges(6,2);
t47 = qJDD(5) + t48;
t13 = m(6) * (t92 * t15 - t87 * t16) - t21 * mrSges(6,3) + t47 * mrSges(6,1) - t41 * t29 + t60 * t32;
t14 = m(6) * (t87 * t15 + t92 * t16) + t20 * mrSges(6,3) - t47 * mrSges(6,2) + t40 * t29 - t60 * t33;
t42 = -t54 * mrSges(5,1) + t55 * mrSges(5,2);
t10 = m(5) * t105 + t48 * mrSges(5,1) - t31 * mrSges(5,3) + t92 * t13 + t87 * t14 - t55 * t42 + t62 * t44;
t11 = m(5) * t115 - t48 * mrSges(5,2) + t30 * mrSges(5,3) - t87 * t13 + t92 * t14 + t54 * t42 - t62 * t45;
t57 = t79 * mrSges(4,1) - t64 * mrSges(4,3);
t9 = m(4) * t114 - t68 * mrSges(4,2) + t49 * mrSges(4,3) - t88 * t10 + t93 * t11 + t63 * t51 - t79 * t57;
t4 = m(3) * (-g(3) * t117 + t113) - t76 * mrSges(3,3) - t81 * mrSges(3,2) + t73 * t107 - t82 * t69 + t94 * t9 - t89 * t12;
t70 = -t82 * mrSges(3,2) + mrSges(3,3) * t107;
t6 = m(3) * (-t85 * t71 - t119) + t75 * mrSges(3,2) + t76 * mrSges(3,1) + t89 * t9 + t94 * t12 + (t69 * t90 - t70 * t95) * t112;
t99 = m(4) * t37 - t49 * mrSges(4,1) + t50 * mrSges(4,2) + t93 * t10 + t88 * t11 - t63 * t56 + t64 * t57;
t8 = m(3) * t102 + t81 * mrSges(3,1) - t75 * mrSges(3,3) - t73 * t108 + t82 * t70 - t99;
t110 = t8 * t116 + t4 * t117 + t86 * t6;
t2 = m(2) * t103 - t97 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t95 * t4 - t90 * t8;
t1 = m(2) * t106 + qJDD(1) * mrSges(2,1) - t97 * mrSges(2,2) - t85 * t6 + (t90 * t4 + t95 * t8) * t86;
t3 = [-m(1) * g(1) - t91 * t1 + t96 * t2, t2, t4, t9, t11, t14; -m(1) * g(2) + t96 * t1 + t91 * t2, t1, t8, t12, t10, t13; (-m(1) - m(2)) * g(3) + t110, -m(2) * g(3) + t110, t6, t99, t98, -t100;];
f_new = t3;
