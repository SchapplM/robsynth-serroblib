% Calculate vector of cutting forces with Newton-Euler
% S5RRPRR14
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPRR14_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR14_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR14_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR14_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR14_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR14_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR14_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:35:58
% EndTime: 2019-12-31 20:36:05
% DurationCPUTime: 3.45s
% Computational Cost: add. (42526->173), mult. (97309->242), div. (0->0), fcn. (76110->12), ass. (0->94)
t87 = cos(pkin(5));
t118 = t87 * g(3);
t91 = sin(qJ(1));
t95 = cos(qJ(1));
t105 = t91 * g(1) - t95 * g(2);
t85 = sin(pkin(5));
t96 = qJD(1) ^ 2;
t70 = t96 * t85 * pkin(7) + qJDD(1) * pkin(1) + t105;
t117 = t70 * t87;
t90 = sin(qJ(2));
t116 = t85 * t90;
t103 = -t95 * g(1) - t91 * g(2);
t110 = qJDD(1) * t85;
t71 = -t96 * pkin(1) + pkin(7) * t110 + t103;
t94 = cos(qJ(2));
t101 = -t90 * t71 + (-g(3) * t85 + t117) * t94;
t112 = qJD(1) * t85;
t107 = t90 * t112;
t111 = qJD(1) * t94;
t106 = t85 * t111;
t81 = t87 * qJD(1) + qJD(2);
t69 = -t81 * mrSges(3,2) + mrSges(3,3) * t106;
t73 = (-mrSges(3,1) * t94 + mrSges(3,2) * t90) * t112;
t74 = (qJD(2) * t111 + qJDD(1) * t90) * t85;
t80 = t87 * qJDD(1) + qJDD(2);
t113 = t90 * t117 + t94 * t71;
t72 = (-pkin(2) * t94 - qJ(3) * t90) * t112;
t79 = t81 ^ 2;
t42 = -t79 * pkin(2) + t80 * qJ(3) + (-g(3) * t90 + t72 * t111) * t85 + t113;
t75 = -qJD(2) * t107 + t94 * t110;
t43 = -t75 * pkin(2) - t118 - t74 * qJ(3) + (-t70 + (pkin(2) * t90 - qJ(3) * t94) * t81 * qJD(1)) * t85;
t84 = sin(pkin(10));
t86 = cos(pkin(10));
t64 = t86 * t107 + t84 * t81;
t104 = -0.2e1 * qJD(3) * t64 - t84 * t42 + t86 * t43;
t56 = t86 * t74 + t84 * t80;
t63 = -t84 * t107 + t86 * t81;
t21 = (-t63 * t106 - t56) * pkin(8) + (t63 * t64 - t75) * pkin(3) + t104;
t108 = 0.2e1 * qJD(3) * t63 + t86 * t42 + t84 * t43;
t55 = -t84 * t74 + t86 * t80;
t57 = -pkin(3) * t106 - t64 * pkin(8);
t62 = t63 ^ 2;
t23 = -t62 * pkin(3) + t55 * pkin(8) + t57 * t106 + t108;
t89 = sin(qJ(4));
t93 = cos(qJ(4));
t114 = t89 * t21 + t93 * t23;
t50 = t93 * t63 - t89 * t64;
t51 = t89 * t63 + t93 * t64;
t37 = -t50 * pkin(4) - t51 * pkin(9);
t67 = qJDD(4) - t75;
t77 = qJD(4) - t106;
t76 = t77 ^ 2;
t18 = -t76 * pkin(4) + t67 * pkin(9) + t50 * t37 + t114;
t31 = -t51 * qJD(4) + t93 * t55 - t89 * t56;
t32 = t50 * qJD(4) + t89 * t55 + t93 * t56;
t41 = -t80 * pkin(2) - t79 * qJ(3) + t72 * t107 + qJDD(3) - t101;
t98 = -t55 * pkin(3) - t62 * pkin(8) + t64 * t57 + t41;
t19 = (-t50 * t77 - t32) * pkin(9) + (t51 * t77 - t31) * pkin(4) + t98;
t88 = sin(qJ(5));
t92 = cos(qJ(5));
t44 = -t88 * t51 + t92 * t77;
t25 = t44 * qJD(5) + t92 * t32 + t88 * t67;
t45 = t92 * t51 + t88 * t77;
t28 = -t44 * mrSges(6,1) + t45 * mrSges(6,2);
t30 = qJDD(5) - t31;
t49 = qJD(5) - t50;
t33 = -t49 * mrSges(6,2) + t44 * mrSges(6,3);
t15 = m(6) * (-t88 * t18 + t92 * t19) - t25 * mrSges(6,3) + t30 * mrSges(6,1) - t45 * t28 + t49 * t33;
t24 = -t45 * qJD(5) - t88 * t32 + t92 * t67;
t34 = t49 * mrSges(6,1) - t45 * mrSges(6,3);
t16 = m(6) * (t92 * t18 + t88 * t19) + t24 * mrSges(6,3) - t30 * mrSges(6,2) + t44 * t28 - t49 * t34;
t46 = -t77 * mrSges(5,2) + t50 * mrSges(5,3);
t47 = t77 * mrSges(5,1) - t51 * mrSges(5,3);
t100 = -m(5) * t98 + t31 * mrSges(5,1) - t32 * mrSges(5,2) - t92 * t15 - t88 * t16 + t50 * t46 - t51 * t47;
t53 = mrSges(4,2) * t106 + t63 * mrSges(4,3);
t54 = -mrSges(4,1) * t106 - t64 * mrSges(4,3);
t97 = m(4) * t41 - t55 * mrSges(4,1) + t56 * mrSges(4,2) - t63 * t53 + t64 * t54 - t100;
t10 = m(3) * t101 + t80 * mrSges(3,1) - t74 * mrSges(3,3) - t73 * t107 + t81 * t69 - t97;
t115 = t94 * t10;
t68 = t81 * mrSges(3,1) - mrSges(3,3) * t107;
t36 = -t50 * mrSges(5,1) + t51 * mrSges(5,2);
t11 = m(5) * t114 - t67 * mrSges(5,2) + t31 * mrSges(5,3) - t88 * t15 + t92 * t16 + t50 * t36 - t77 * t47;
t102 = t93 * t21 - t89 * t23;
t99 = m(6) * (-t67 * pkin(4) - t76 * pkin(9) + t51 * t37 - t102) - t24 * mrSges(6,1) + t25 * mrSges(6,2) - t44 * t33 + t45 * t34;
t12 = m(5) * t102 + t67 * mrSges(5,1) - t32 * mrSges(5,3) - t51 * t36 + t77 * t46 - t99;
t52 = -t63 * mrSges(4,1) + t64 * mrSges(4,2);
t7 = m(4) * t104 - t75 * mrSges(4,1) - t56 * mrSges(4,3) - t53 * t106 + t89 * t11 + t93 * t12 - t64 * t52;
t8 = m(4) * t108 + t75 * mrSges(4,2) + t55 * mrSges(4,3) + t54 * t106 + t93 * t11 - t89 * t12 + t63 * t52;
t4 = m(3) * (-g(3) * t116 + t113) + t75 * mrSges(3,3) - t80 * mrSges(3,2) + t73 * t106 - t81 * t68 + t86 * t8 - t84 * t7;
t6 = m(3) * (-t85 * t70 - t118) + t74 * mrSges(3,2) - t75 * mrSges(3,1) + t84 * t8 + t86 * t7 + (t68 * t90 - t69 * t94) * t112;
t109 = t85 * t115 + t4 * t116 + t87 * t6;
t2 = m(2) * t103 - t96 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t90 * t10 + t94 * t4;
t1 = m(2) * t105 + qJDD(1) * mrSges(2,1) - t96 * mrSges(2,2) - t85 * t6 + (t90 * t4 + t115) * t87;
t3 = [-m(1) * g(1) - t91 * t1 + t95 * t2, t2, t4, t8, t11, t16; -m(1) * g(2) + t95 * t1 + t91 * t2, t1, t10, t7, t12, t15; (-m(1) - m(2)) * g(3) + t109, -m(2) * g(3) + t109, t6, t97, -t100, t99;];
f_new = t3;
