% Calculate vector of cutting forces with Newton-Euler
% S5RRRRR10
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
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRRR10_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR10_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR10_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR10_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR10_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR10_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR10_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:32:12
% EndTime: 2019-12-31 22:32:19
% DurationCPUTime: 3.65s
% Computational Cost: add. (48311->174), mult. (103815->240), div. (0->0), fcn. (81875->12), ass. (0->96)
t85 = cos(pkin(5));
t118 = t85 * g(3);
t90 = sin(qJ(1));
t95 = cos(qJ(1));
t105 = t90 * g(1) - t95 * g(2);
t84 = sin(pkin(5));
t96 = qJD(1) ^ 2;
t69 = t96 * t84 * pkin(7) + qJDD(1) * pkin(1) + t105;
t117 = t69 * t85;
t89 = sin(qJ(2));
t116 = t84 * t89;
t103 = -t95 * g(1) - t90 * g(2);
t109 = qJDD(1) * t84;
t70 = -t96 * pkin(1) + pkin(7) * t109 + t103;
t94 = cos(qJ(2));
t101 = -t89 * t70 + (-g(3) * t84 + t117) * t94;
t111 = qJD(1) * t84;
t107 = t89 * t111;
t110 = qJD(1) * t94;
t106 = t84 * t110;
t81 = t85 * qJD(1) + qJD(2);
t68 = -t81 * mrSges(3,2) + mrSges(3,3) * t106;
t71 = (-mrSges(3,1) * t94 + mrSges(3,2) * t89) * t111;
t73 = (qJD(2) * t110 + qJDD(1) * t89) * t84;
t80 = t85 * qJDD(1) + qJDD(2);
t112 = t89 * t117 + t94 * t70;
t72 = (-pkin(2) * t94 - pkin(8) * t89) * t111;
t79 = t81 ^ 2;
t42 = -t79 * pkin(2) + t80 * pkin(8) + (-g(3) * t89 + t110 * t72) * t84 + t112;
t74 = -qJD(2) * t107 + t109 * t94;
t43 = -t74 * pkin(2) - t73 * pkin(8) - t118 + (-t69 + (pkin(2) * t89 - pkin(8) * t94) * t81 * qJD(1)) * t84;
t88 = sin(qJ(3));
t93 = cos(qJ(3));
t104 = -t88 * t42 + t93 * t43;
t61 = -t107 * t88 + t93 * t81;
t50 = t61 * qJD(3) + t93 * t73 + t88 * t80;
t62 = t107 * t93 + t88 * t81;
t66 = qJDD(3) - t74;
t77 = qJD(3) - t106;
t21 = (t61 * t77 - t50) * pkin(9) + (t61 * t62 + t66) * pkin(3) + t104;
t113 = t93 * t42 + t88 * t43;
t49 = -t62 * qJD(3) - t88 * t73 + t93 * t80;
t57 = t77 * pkin(3) - t62 * pkin(9);
t60 = t61 ^ 2;
t23 = -t60 * pkin(3) + t49 * pkin(9) - t77 * t57 + t113;
t87 = sin(qJ(4));
t92 = cos(qJ(4));
t114 = t87 * t21 + t92 * t23;
t52 = t92 * t61 - t87 * t62;
t53 = t87 * t61 + t92 * t62;
t37 = -t52 * pkin(4) - t53 * pkin(10);
t65 = qJDD(4) + t66;
t76 = qJD(4) + t77;
t75 = t76 ^ 2;
t18 = -t75 * pkin(4) + t65 * pkin(10) + t52 * t37 + t114;
t30 = -t53 * qJD(4) + t92 * t49 - t87 * t50;
t31 = t52 * qJD(4) + t87 * t49 + t92 * t50;
t41 = -t80 * pkin(2) - t79 * pkin(8) + t72 * t107 - t101;
t98 = -t49 * pkin(3) - t60 * pkin(9) + t62 * t57 + t41;
t19 = (-t52 * t76 - t31) * pkin(10) + (t53 * t76 - t30) * pkin(4) + t98;
t86 = sin(qJ(5));
t91 = cos(qJ(5));
t44 = -t86 * t53 + t91 * t76;
t25 = t44 * qJD(5) + t91 * t31 + t86 * t65;
t29 = qJDD(5) - t30;
t45 = t91 * t53 + t86 * t76;
t32 = -t44 * mrSges(6,1) + t45 * mrSges(6,2);
t51 = qJD(5) - t52;
t33 = -t51 * mrSges(6,2) + t44 * mrSges(6,3);
t15 = m(6) * (-t86 * t18 + t91 * t19) - t25 * mrSges(6,3) + t29 * mrSges(6,1) - t45 * t32 + t51 * t33;
t24 = -t45 * qJD(5) - t86 * t31 + t91 * t65;
t34 = t51 * mrSges(6,1) - t45 * mrSges(6,3);
t16 = m(6) * (t91 * t18 + t86 * t19) + t24 * mrSges(6,3) - t29 * mrSges(6,2) + t44 * t32 - t51 * t34;
t46 = -t76 * mrSges(5,2) + t52 * mrSges(5,3);
t47 = t76 * mrSges(5,1) - t53 * mrSges(5,3);
t100 = -m(5) * t98 + t30 * mrSges(5,1) - t31 * mrSges(5,2) - t91 * t15 - t86 * t16 + t52 * t46 - t53 * t47;
t55 = -t77 * mrSges(4,2) + t61 * mrSges(4,3);
t56 = t77 * mrSges(4,1) - t62 * mrSges(4,3);
t97 = m(4) * t41 - t49 * mrSges(4,1) + t50 * mrSges(4,2) - t61 * t55 + t62 * t56 - t100;
t10 = m(3) * t101 + t80 * mrSges(3,1) - t73 * mrSges(3,3) - t107 * t71 + t81 * t68 - t97;
t115 = t94 * t10;
t67 = t81 * mrSges(3,1) - mrSges(3,3) * t107;
t36 = -t52 * mrSges(5,1) + t53 * mrSges(5,2);
t11 = m(5) * t114 - t65 * mrSges(5,2) + t30 * mrSges(5,3) - t86 * t15 + t91 * t16 + t52 * t36 - t76 * t47;
t102 = t92 * t21 - t87 * t23;
t99 = m(6) * (-t65 * pkin(4) - t75 * pkin(10) + t53 * t37 - t102) - t24 * mrSges(6,1) + t25 * mrSges(6,2) - t44 * t33 + t45 * t34;
t12 = m(5) * t102 + t65 * mrSges(5,1) - t31 * mrSges(5,3) - t53 * t36 + t76 * t46 - t99;
t54 = -t61 * mrSges(4,1) + t62 * mrSges(4,2);
t7 = m(4) * t104 + t66 * mrSges(4,1) - t50 * mrSges(4,3) + t87 * t11 + t92 * t12 - t62 * t54 + t77 * t55;
t8 = m(4) * t113 - t66 * mrSges(4,2) + t49 * mrSges(4,3) + t92 * t11 - t87 * t12 + t61 * t54 - t77 * t56;
t4 = m(3) * (-g(3) * t116 + t112) + t74 * mrSges(3,3) - t80 * mrSges(3,2) + t71 * t106 - t81 * t67 + t93 * t8 - t88 * t7;
t6 = m(3) * (-t84 * t69 - t118) + t73 * mrSges(3,2) - t74 * mrSges(3,1) + t88 * t8 + t93 * t7 + (t67 * t89 - t68 * t94) * t111;
t108 = t84 * t115 + t4 * t116 + t85 * t6;
t2 = m(2) * t103 - t96 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t89 * t10 + t94 * t4;
t1 = m(2) * t105 + qJDD(1) * mrSges(2,1) - t96 * mrSges(2,2) - t84 * t6 + (t89 * t4 + t115) * t85;
t3 = [-m(1) * g(1) - t90 * t1 + t95 * t2, t2, t4, t8, t11, t16; -m(1) * g(2) + t95 * t1 + t90 * t2, t1, t10, t7, t12, t15; (-m(1) - m(2)) * g(3) + t108, -m(2) * g(3) + t108, t6, t97, -t100, t99;];
f_new = t3;
