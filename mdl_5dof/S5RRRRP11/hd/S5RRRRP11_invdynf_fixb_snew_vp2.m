% Calculate vector of cutting forces with Newton-Euler
% S5RRRRP11
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRRP11_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP11_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP11_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP11_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP11_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP11_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP11_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:14:32
% EndTime: 2019-12-31 22:14:37
% DurationCPUTime: 1.89s
% Computational Cost: add. (22416->170), mult. (47970->228), div. (0->0), fcn. (36392->10), ass. (0->88)
t77 = sin(pkin(5));
t81 = sin(qJ(2));
t84 = cos(qJ(2));
t99 = qJD(1) * qJD(2);
t68 = (-qJDD(1) * t84 + t81 * t99) * t77;
t78 = cos(pkin(5));
t74 = t78 * qJD(1) + qJD(2);
t80 = sin(qJ(3));
t83 = cos(qJ(3));
t102 = qJD(1) * t77;
t96 = t81 * t102;
t56 = t83 * t74 - t80 * t96;
t57 = t80 * t74 + t83 * t96;
t46 = -t56 * pkin(3) - t57 * pkin(9);
t60 = qJDD(3) + t68;
t101 = qJD(1) * t84;
t95 = t77 * t101;
t71 = qJD(3) - t95;
t70 = t71 ^ 2;
t115 = pkin(7) * t77;
t86 = qJD(1) ^ 2;
t82 = sin(qJ(1));
t85 = cos(qJ(1));
t94 = t82 * g(1) - t85 * g(2);
t63 = qJDD(1) * pkin(1) + t86 * t115 + t94;
t112 = t63 * t78;
t92 = -t85 * g(1) - t82 * g(2);
t64 = -t86 * pkin(1) + qJDD(1) * t115 + t92;
t103 = t81 * t112 + t84 * t64;
t66 = (-pkin(2) * t84 - pkin(8) * t81) * t102;
t72 = t74 ^ 2;
t73 = t78 * qJDD(1) + qJDD(2);
t30 = -t72 * pkin(2) + t73 * pkin(8) + (-g(3) * t81 + t66 * t101) * t77 + t103;
t114 = t78 * g(3);
t67 = (qJDD(1) * t81 + t84 * t99) * t77;
t31 = t68 * pkin(2) - t67 * pkin(8) - t114 + (-t63 + (pkin(2) * t81 - pkin(8) * t84) * t74 * qJD(1)) * t77;
t93 = -t80 * t30 + t83 * t31;
t18 = -t60 * pkin(3) - t70 * pkin(9) + t57 * t46 - t93;
t113 = cos(qJ(4));
t44 = t56 * qJD(3) + t83 * t67 + t80 * t73;
t79 = sin(qJ(4));
t48 = t113 * t57 + t79 * t71;
t23 = t48 * qJD(4) - t113 * t60 + t79 * t44;
t47 = -t113 * t71 + t79 * t57;
t24 = -t47 * qJD(4) + t113 * t44 + t79 * t60;
t54 = qJD(4) - t56;
t38 = -t54 * mrSges(5,2) - t47 * mrSges(5,3);
t39 = t54 * mrSges(5,1) - t48 * mrSges(5,3);
t40 = -t54 * mrSges(6,1) + t48 * mrSges(6,2);
t37 = -t47 * mrSges(6,2) + t54 * mrSges(6,3);
t97 = m(6) * (-0.2e1 * qJD(5) * t48 + (t47 * t54 - t24) * qJ(5) + (t48 * t54 + t23) * pkin(4) + t18) + t23 * mrSges(6,1) + t47 * t37;
t117 = m(5) * t18 + t23 * mrSges(5,1) + (t39 - t40) * t48 + (mrSges(5,2) - mrSges(6,3)) * t24 + t47 * t38 + t97;
t33 = t47 * pkin(4) - t48 * qJ(5);
t43 = -t57 * qJD(3) - t80 * t67 + t83 * t73;
t42 = qJDD(4) - t43;
t53 = t54 ^ 2;
t106 = t83 * t30 + t80 * t31;
t19 = -t70 * pkin(3) + t60 * pkin(9) + t56 * t46 + t106;
t110 = t77 * t84;
t91 = -g(3) * t110 + t84 * t112 - t81 * t64;
t29 = -t73 * pkin(2) - t72 * pkin(8) + t66 * t96 - t91;
t21 = (-t56 * t71 - t44) * pkin(9) + (t57 * t71 - t43) * pkin(3) + t29;
t90 = t113 * t21 - t79 * t19;
t116 = m(6) * (-t42 * pkin(4) - t53 * qJ(5) + t48 * t33 + qJDD(5) - t90);
t111 = t77 * t81;
t108 = -mrSges(5,3) - mrSges(6,2);
t107 = t113 * t19 + t79 * t21;
t34 = t47 * mrSges(6,1) - t48 * mrSges(6,3);
t105 = -t47 * mrSges(5,1) - t48 * mrSges(5,2) - t34;
t45 = -t56 * mrSges(4,1) + t57 * mrSges(4,2);
t49 = -t71 * mrSges(4,2) + t56 * mrSges(4,3);
t10 = m(4) * t93 + t60 * mrSges(4,1) - t44 * mrSges(4,3) - t57 * t45 + t71 * t49 - t117;
t61 = t74 * mrSges(3,1) - mrSges(3,3) * t96;
t65 = (-mrSges(3,1) * t84 + mrSges(3,2) * t81) * t102;
t98 = m(6) * (-t53 * pkin(4) + t42 * qJ(5) + 0.2e1 * qJD(5) * t54 - t47 * t33 + t107) + t54 * t40 + t42 * mrSges(6,3);
t11 = m(5) * t107 - t42 * mrSges(5,2) + t105 * t47 + t108 * t23 - t54 * t39 + t98;
t12 = m(5) * t90 - t116 + (t38 + t37) * t54 + t105 * t48 + (mrSges(5,1) + mrSges(6,1)) * t42 + t108 * t24;
t50 = t71 * mrSges(4,1) - t57 * mrSges(4,3);
t9 = m(4) * t106 - t60 * mrSges(4,2) + t43 * mrSges(4,3) + t113 * t11 - t79 * t12 + t56 * t45 - t71 * t50;
t4 = m(3) * (-g(3) * t111 + t103) - t68 * mrSges(3,3) - t73 * mrSges(3,2) + t65 * t95 - t74 * t61 + t83 * t9 - t80 * t10;
t62 = -t74 * mrSges(3,2) + mrSges(3,3) * t95;
t6 = m(3) * (-t77 * t63 - t114) + t67 * mrSges(3,2) + t68 * mrSges(3,1) + t80 * t9 + t83 * t10 + (t61 * t81 - t62 * t84) * t102;
t87 = m(4) * t29 - t43 * mrSges(4,1) + t44 * mrSges(4,2) + t79 * t11 + t113 * t12 - t56 * t49 + t57 * t50;
t8 = m(3) * t91 + t73 * mrSges(3,1) - t67 * mrSges(3,3) + t74 * t62 - t65 * t96 - t87;
t100 = t8 * t110 + t4 * t111 + t78 * t6;
t2 = m(2) * t92 - t86 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t84 * t4 - t81 * t8;
t1 = m(2) * t94 + qJDD(1) * mrSges(2,1) - t86 * mrSges(2,2) - t77 * t6 + (t81 * t4 + t84 * t8) * t78;
t3 = [-m(1) * g(1) - t82 * t1 + t85 * t2, t2, t4, t9, t11, -t23 * mrSges(6,2) - t47 * t34 + t98; -m(1) * g(2) + t85 * t1 + t82 * t2, t1, t8, t10, t12, -t24 * mrSges(6,3) - t48 * t40 + t97; (-m(1) - m(2)) * g(3) + t100, -m(2) * g(3) + t100, t6, t87, t117, -t42 * mrSges(6,1) + t24 * mrSges(6,2) + t48 * t34 - t54 * t37 + t116;];
f_new = t3;
