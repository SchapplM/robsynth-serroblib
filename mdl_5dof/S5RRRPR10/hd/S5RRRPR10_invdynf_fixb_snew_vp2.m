% Calculate vector of cutting forces with Newton-Euler
% S5RRRPR10
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
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRPR10_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR10_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR10_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR10_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR10_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR10_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR10_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:26:50
% EndTime: 2019-12-31 21:26:57
% DurationCPUTime: 3.50s
% Computational Cost: add. (44996->174), mult. (99767->243), div. (0->0), fcn. (77963->12), ass. (0->96)
t108 = qJD(1) * qJD(2);
t83 = sin(pkin(5));
t88 = sin(qJ(2));
t92 = cos(qJ(2));
t74 = (-qJDD(1) * t92 + t88 * t108) * t83;
t119 = 2 * qJD(4);
t118 = pkin(7) * t83;
t85 = cos(pkin(5));
t117 = t85 * g(3);
t89 = sin(qJ(1));
t93 = cos(qJ(1));
t104 = t89 * g(1) - t93 * g(2);
t94 = qJD(1) ^ 2;
t69 = qJDD(1) * pkin(1) + t94 * t118 + t104;
t116 = t69 * t85;
t115 = t83 * t88;
t102 = -t93 * g(1) - t89 * g(2);
t70 = -t94 * pkin(1) + qJDD(1) * t118 + t102;
t100 = -t88 * t70 + (-g(3) * t83 + t116) * t92;
t111 = qJD(1) * t83;
t106 = t88 * t111;
t110 = qJD(1) * t92;
t105 = t83 * t110;
t79 = t85 * qJD(1) + qJD(2);
t68 = -t79 * mrSges(3,2) + mrSges(3,3) * t105;
t71 = (-mrSges(3,1) * t92 + mrSges(3,2) * t88) * t111;
t73 = (qJDD(1) * t88 + t92 * t108) * t83;
t78 = t85 * qJDD(1) + qJDD(2);
t72 = (-pkin(2) * t92 - pkin(8) * t88) * t111;
t77 = t79 ^ 2;
t41 = -t78 * pkin(2) - t77 * pkin(8) + t72 * t106 - t100;
t87 = sin(qJ(3));
t91 = cos(qJ(3));
t63 = t91 * t106 + t87 * t79;
t49 = -t63 * qJD(3) - t87 * t73 + t91 * t78;
t62 = -t87 * t106 + t91 * t79;
t50 = t62 * qJD(3) + t91 * t73 + t87 * t78;
t76 = qJD(3) - t105;
t56 = -t76 * mrSges(4,2) + t62 * mrSges(4,3);
t58 = t76 * mrSges(4,1) - t63 * mrSges(4,3);
t112 = t88 * t116 + t92 * t70;
t42 = -t77 * pkin(2) + t78 * pkin(8) + (-g(3) * t88 + t72 * t110) * t83 + t112;
t43 = t74 * pkin(2) - t73 * pkin(8) - t117 + (-t69 + (pkin(2) * t88 - pkin(8) * t92) * t79 * qJD(1)) * t83;
t103 = -t87 * t42 + t91 * t43;
t66 = qJDD(3) + t74;
t21 = (t62 * t76 - t50) * qJ(4) + (t62 * t63 + t66) * pkin(3) + t103;
t113 = t91 * t42 + t87 * t43;
t57 = t76 * pkin(3) - t63 * qJ(4);
t61 = t62 ^ 2;
t23 = -t61 * pkin(3) + t49 * qJ(4) - t76 * t57 + t113;
t82 = sin(pkin(10));
t84 = cos(pkin(10));
t53 = t84 * t62 - t82 * t63;
t107 = t53 * t119 + t82 * t21 + t84 * t23;
t54 = t82 * t62 + t84 * t63;
t37 = -t53 * pkin(4) - t54 * pkin(9);
t75 = t76 ^ 2;
t18 = -t75 * pkin(4) + t66 * pkin(9) + t53 * t37 + t107;
t33 = t84 * t49 - t82 * t50;
t34 = t82 * t49 + t84 * t50;
t96 = -t49 * pkin(3) - t61 * qJ(4) + t63 * t57 + qJDD(4) + t41;
t19 = (-t53 * t76 - t34) * pkin(9) + (t54 * t76 - t33) * pkin(4) + t96;
t86 = sin(qJ(5));
t90 = cos(qJ(5));
t44 = -t86 * t54 + t90 * t76;
t27 = t44 * qJD(5) + t90 * t34 + t86 * t66;
t45 = t90 * t54 + t86 * t76;
t28 = -t44 * mrSges(6,1) + t45 * mrSges(6,2);
t52 = qJD(5) - t53;
t29 = -t52 * mrSges(6,2) + t44 * mrSges(6,3);
t32 = qJDD(5) - t33;
t15 = m(6) * (-t86 * t18 + t90 * t19) - t27 * mrSges(6,3) + t32 * mrSges(6,1) - t45 * t28 + t52 * t29;
t26 = -t45 * qJD(5) - t86 * t34 + t90 * t66;
t30 = t52 * mrSges(6,1) - t45 * mrSges(6,3);
t16 = m(6) * (t90 * t18 + t86 * t19) + t26 * mrSges(6,3) - t32 * mrSges(6,2) + t44 * t28 - t52 * t30;
t46 = -t76 * mrSges(5,2) + t53 * mrSges(5,3);
t47 = t76 * mrSges(5,1) - t54 * mrSges(5,3);
t98 = -m(5) * t96 + t33 * mrSges(5,1) - t34 * mrSges(5,2) - t90 * t15 - t86 * t16 + t53 * t46 - t54 * t47;
t95 = m(4) * t41 - t49 * mrSges(4,1) + t50 * mrSges(4,2) - t62 * t56 + t63 * t58 - t98;
t10 = m(3) * t100 + t78 * mrSges(3,1) - t73 * mrSges(3,3) - t71 * t106 + t79 * t68 - t95;
t114 = t92 * t10;
t67 = t79 * mrSges(3,1) - mrSges(3,3) * t106;
t36 = -t53 * mrSges(5,1) + t54 * mrSges(5,2);
t11 = m(5) * t107 - t66 * mrSges(5,2) + t33 * mrSges(5,3) - t86 * t15 + t90 * t16 + t53 * t36 - t76 * t47;
t101 = -t84 * t21 + t82 * t23;
t97 = m(6) * (-t66 * pkin(4) - t75 * pkin(9) + (t119 + t37) * t54 + t101) - t26 * mrSges(6,1) + t27 * mrSges(6,2) - t44 * t29 + t45 * t30;
t12 = m(5) * (-0.2e1 * qJD(4) * t54 - t101) - t34 * mrSges(5,3) + t66 * mrSges(5,1) - t54 * t36 + t76 * t46 - t97;
t55 = -t62 * mrSges(4,1) + t63 * mrSges(4,2);
t7 = m(4) * t103 + t66 * mrSges(4,1) - t50 * mrSges(4,3) + t82 * t11 + t84 * t12 - t63 * t55 + t76 * t56;
t8 = m(4) * t113 - t66 * mrSges(4,2) + t49 * mrSges(4,3) + t84 * t11 - t82 * t12 + t62 * t55 - t76 * t58;
t4 = m(3) * (-g(3) * t115 + t112) - t74 * mrSges(3,3) - t78 * mrSges(3,2) + t71 * t105 - t79 * t67 + t91 * t8 - t87 * t7;
t6 = m(3) * (-t83 * t69 - t117) + t73 * mrSges(3,2) + t74 * mrSges(3,1) + t87 * t8 + t91 * t7 + (t67 * t88 - t68 * t92) * t111;
t109 = t83 * t114 + t4 * t115 + t85 * t6;
t2 = m(2) * t102 - t94 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t88 * t10 + t92 * t4;
t1 = m(2) * t104 + qJDD(1) * mrSges(2,1) - t94 * mrSges(2,2) - t83 * t6 + (t88 * t4 + t114) * t85;
t3 = [-m(1) * g(1) - t89 * t1 + t93 * t2, t2, t4, t8, t11, t16; -m(1) * g(2) + t93 * t1 + t89 * t2, t1, t10, t7, t12, t15; (-m(1) - m(2)) * g(3) + t109, -m(2) * g(3) + t109, t6, t95, -t98, t97;];
f_new = t3;
