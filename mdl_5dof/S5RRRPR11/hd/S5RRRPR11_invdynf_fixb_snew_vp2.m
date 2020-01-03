% Calculate vector of cutting forces with Newton-Euler
% S5RRRPR11
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRPR11_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR11_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR11_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR11_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR11_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR11_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR11_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:32:28
% EndTime: 2019-12-31 21:32:31
% DurationCPUTime: 0.99s
% Computational Cost: add. (9793->163), mult. (19211->202), div. (0->0), fcn. (12048->8), ass. (0->81)
t81 = sin(qJ(2));
t106 = qJD(1) * t81;
t114 = cos(qJ(3));
t80 = sin(qJ(3));
t62 = -t114 * qJD(2) + t80 * t106;
t84 = cos(qJ(2));
t105 = t84 * qJD(1);
t73 = qJD(3) - t105;
t113 = t62 * t73;
t104 = qJD(1) * qJD(2);
t100 = t84 * t104;
t66 = t81 * qJDD(1) + t100;
t40 = -t62 * qJD(3) + t80 * qJDD(2) + t114 * t66;
t87 = qJD(1) ^ 2;
t82 = sin(qJ(1));
t85 = cos(qJ(1));
t97 = -t85 * g(1) - t82 * g(2);
t56 = -t87 * pkin(1) + qJDD(1) * pkin(6) + t97;
t107 = -t84 * g(3) - t81 * t56;
t65 = (-pkin(2) * t84 - pkin(7) * t81) * qJD(1);
t86 = qJD(2) ^ 2;
t90 = qJDD(2) * pkin(2) + t86 * pkin(7) - t65 * t106 + t107;
t119 = (-t40 + t113) * qJ(4) - t90;
t63 = t80 * qJD(2) + t114 * t106;
t45 = t62 * mrSges(5,1) - t63 * mrSges(5,3);
t109 = -t62 * mrSges(4,1) - t63 * mrSges(4,2) - t45;
t101 = t81 * t104;
t103 = t82 * g(1) - t85 * g(2);
t55 = -qJDD(1) * pkin(1) - t87 * pkin(6) - t103;
t67 = t84 * qJDD(1) - t101;
t28 = (-t66 - t100) * pkin(7) + (-t67 + t101) * pkin(2) + t55;
t102 = -t81 * g(3) + t84 * t56;
t32 = -t86 * pkin(2) + qJDD(2) * pkin(7) + t65 * t105 + t102;
t110 = t114 * t32 + t80 * t28;
t111 = -mrSges(4,3) - mrSges(5,2);
t39 = t63 * qJD(3) - t114 * qJDD(2) + t80 * t66;
t48 = t73 * mrSges(4,1) - t63 * mrSges(4,3);
t61 = qJDD(3) - t67;
t44 = t62 * pkin(3) - t63 * qJ(4);
t72 = t73 ^ 2;
t96 = t114 * t28 - t80 * t32;
t19 = -t61 * pkin(3) - t72 * qJ(4) + t63 * t44 + qJDD(4) - t96;
t12 = (-t40 - t113) * pkin(8) + (t62 * t63 - t61) * pkin(4) + t19;
t51 = -t73 * pkin(4) - t63 * pkin(8);
t60 = t62 ^ 2;
t116 = 2 * qJD(4);
t93 = -t72 * pkin(3) + t61 * qJ(4) + t73 * t116 - t62 * t44 + t110;
t14 = -t60 * pkin(4) + t39 * pkin(8) + t73 * t51 + t93;
t79 = sin(qJ(5));
t83 = cos(qJ(5));
t41 = t83 * t62 - t79 * t63;
t23 = t41 * qJD(5) + t79 * t39 + t83 * t40;
t42 = t79 * t62 + t83 * t63;
t26 = -t41 * mrSges(6,1) + t42 * mrSges(6,2);
t71 = qJD(5) - t73;
t33 = -t71 * mrSges(6,2) + t41 * mrSges(6,3);
t59 = qJDD(5) - t61;
t10 = m(6) * (t83 * t12 - t79 * t14) - t23 * mrSges(6,3) + t59 * mrSges(6,1) - t42 * t26 + t71 * t33;
t22 = -t42 * qJD(5) + t83 * t39 - t79 * t40;
t34 = t71 * mrSges(6,1) - t42 * mrSges(6,3);
t11 = m(6) * (t79 * t12 + t83 * t14) + t22 * mrSges(6,3) - t59 * mrSges(6,2) + t41 * t26 - t71 * t34;
t49 = -t73 * mrSges(5,1) + t63 * mrSges(5,2);
t94 = m(5) * t93 + t61 * mrSges(5,3) - t79 * t10 + t83 * t11 + t73 * t49;
t5 = m(4) * t110 - t61 * mrSges(4,2) + t109 * t62 + t111 * t39 - t73 * t48 + t94;
t47 = -t73 * mrSges(4,2) - t62 * mrSges(4,3);
t50 = -t62 * mrSges(5,2) + t73 * mrSges(5,3);
t92 = -m(5) * t19 - t83 * t10 - t79 * t11;
t6 = m(4) * t96 + (t47 + t50) * t73 + t109 * t63 + (mrSges(4,1) + mrSges(5,1)) * t61 + t111 * t40 + t92;
t68 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t106;
t69 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t105;
t118 = m(3) * t55 - t67 * mrSges(3,1) + t66 * mrSges(3,2) + (t68 * t81 - t69 * t84) * qJD(1) + t114 * t6 + t80 * t5;
t98 = m(6) * (-t60 * pkin(8) + (-pkin(3) - pkin(4)) * t39 + (-pkin(3) * t73 + t116 + t51) * t63 - t119) + t23 * mrSges(6,2) - t22 * mrSges(6,1) + t42 * t34 - t41 * t33;
t91 = m(5) * (-0.2e1 * qJD(4) * t63 + (t63 * t73 + t39) * pkin(3) + t119) + t39 * mrSges(5,1) + t62 * t50 - t98;
t117 = -m(4) * t90 + t39 * mrSges(4,1) + (t48 - t49) * t63 + (mrSges(4,2) - mrSges(5,3)) * t40 + t62 * t47 + t91;
t64 = (-mrSges(3,1) * t84 + mrSges(3,2) * t81) * qJD(1);
t4 = m(3) * t102 - qJDD(2) * mrSges(3,2) + t67 * mrSges(3,3) - qJD(2) * t68 + t64 * t105 + t114 * t5 - t80 * t6;
t8 = m(3) * t107 + qJDD(2) * mrSges(3,1) - t66 * mrSges(3,3) + qJD(2) * t69 - t64 * t106 - t117;
t115 = t81 * t4 + t84 * t8;
t2 = m(2) * t103 + qJDD(1) * mrSges(2,1) - t87 * mrSges(2,2) - t118;
t1 = m(2) * t97 - t87 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t84 * t4 - t81 * t8;
t3 = [-m(1) * g(1) + t85 * t1 - t82 * t2, t1, t4, t5, -t39 * mrSges(5,2) - t62 * t45 + t94, t11; -m(1) * g(2) + t82 * t1 + t85 * t2, t2, t8, t6, -t40 * mrSges(5,3) - t63 * t49 + t91, t10; (-m(1) - m(2)) * g(3) + t115, -m(2) * g(3) + t115, t118, t117, -t61 * mrSges(5,1) + t40 * mrSges(5,2) + t63 * t45 - t73 * t50 - t92, t98;];
f_new = t3;
