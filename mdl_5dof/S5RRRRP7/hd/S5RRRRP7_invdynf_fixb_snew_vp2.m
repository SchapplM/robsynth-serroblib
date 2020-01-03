% Calculate vector of cutting forces with Newton-Euler
% S5RRRRP7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRRP7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP7_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP7_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP7_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP7_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP7_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:56:11
% EndTime: 2019-12-31 21:56:14
% DurationCPUTime: 1.18s
% Computational Cost: add. (12393->163), mult. (25069->208), div. (0->0), fcn. (16776->8), ass. (0->78)
t76 = sin(qJ(3));
t77 = sin(qJ(2));
t79 = cos(qJ(3));
t80 = cos(qJ(2));
t57 = (t76 * t77 - t79 * t80) * qJD(1);
t97 = qJD(1) * qJD(2);
t63 = t77 * qJDD(1) + t80 * t97;
t64 = t80 * qJDD(1) - t77 * t97;
t99 = qJD(1) * t77;
t65 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t99;
t98 = qJD(1) * t80;
t66 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t98;
t82 = qJD(1) ^ 2;
t107 = cos(qJ(4));
t58 = (t76 * t80 + t77 * t79) * qJD(1);
t73 = qJD(2) + qJD(3);
t75 = sin(qJ(4));
t50 = -t107 * t73 + t75 * t58;
t51 = t107 * t58 + t75 * t73;
t31 = t50 * mrSges(6,1) - t51 * mrSges(6,3);
t101 = -t50 * mrSges(5,1) - t51 * mrSges(5,2) - t31;
t40 = -t58 * qJD(3) - t76 * t63 + t79 * t64;
t41 = -t57 * qJD(3) + t79 * t63 + t76 * t64;
t67 = qJD(2) * pkin(2) - pkin(7) * t99;
t74 = t80 ^ 2;
t78 = sin(qJ(1));
t81 = cos(qJ(1));
t94 = t78 * g(1) - t81 * g(2);
t88 = -qJDD(1) * pkin(1) - t94;
t84 = -t64 * pkin(2) + t67 * t99 + (-pkin(7) * t74 - pkin(6)) * t82 + t88;
t18 = (t57 * t73 - t41) * pkin(8) + (t58 * t73 - t40) * pkin(3) + t84;
t91 = -t81 * g(1) - t78 * g(2);
t60 = -t82 * pkin(1) + qJDD(1) * pkin(6) + t91;
t106 = t77 * t60;
t108 = pkin(2) * t82;
t33 = qJDD(2) * pkin(2) - t63 * pkin(7) - t106 + (pkin(7) * t97 + t77 * t108 - g(3)) * t80;
t93 = -t77 * g(3) + t80 * t60;
t34 = t64 * pkin(7) - qJD(2) * t67 - t74 * t108 + t93;
t102 = t76 * t33 + t79 * t34;
t49 = t57 * pkin(3) - t58 * pkin(8);
t71 = t73 ^ 2;
t72 = qJDD(2) + qJDD(3);
t21 = -t71 * pkin(3) + t72 * pkin(8) - t57 * t49 + t102;
t103 = t107 * t21 + t75 * t18;
t104 = -mrSges(5,3) - mrSges(6,2);
t23 = t51 * qJD(4) - t107 * t72 + t75 * t41;
t38 = qJDD(4) - t40;
t56 = qJD(4) + t57;
t45 = t56 * mrSges(5,1) - t51 * mrSges(5,3);
t30 = t50 * pkin(4) - t51 * qJ(5);
t46 = -t56 * mrSges(6,1) + t51 * mrSges(6,2);
t55 = t56 ^ 2;
t96 = m(6) * (-t55 * pkin(4) + t38 * qJ(5) + 0.2e1 * qJD(5) * t56 - t50 * t30 + t103) + t56 * t46 + t38 * mrSges(6,3);
t10 = m(5) * t103 - t38 * mrSges(5,2) + t101 * t50 + t104 * t23 - t56 * t45 + t96;
t87 = t107 * t18 - t75 * t21;
t109 = m(6) * (-t38 * pkin(4) - t55 * qJ(5) + t51 * t30 + qJDD(5) - t87);
t24 = -t50 * qJD(4) + t107 * t41 + t75 * t72;
t43 = -t50 * mrSges(6,2) + t56 * mrSges(6,3);
t44 = -t56 * mrSges(5,2) - t50 * mrSges(5,3);
t12 = m(5) * t87 - t109 + (t44 + t43) * t56 + t101 * t51 + (mrSges(5,1) + mrSges(6,1)) * t38 + t104 * t24;
t52 = -t73 * mrSges(4,2) - t57 * mrSges(4,3);
t53 = t73 * mrSges(4,1) - t58 * mrSges(4,3);
t86 = -m(4) * t84 + t40 * mrSges(4,1) - t41 * mrSges(4,2) - t75 * t10 - t107 * t12 - t57 * t52 - t58 * t53;
t112 = (t77 * t65 - t80 * t66) * qJD(1) + m(3) * (-t82 * pkin(6) + t88) - t64 * mrSges(3,1) + t63 * mrSges(3,2) - t86;
t92 = t79 * t33 - t76 * t34;
t20 = -t72 * pkin(3) - t71 * pkin(8) + t58 * t49 - t92;
t95 = m(6) * (-0.2e1 * qJD(5) * t51 + (t50 * t56 - t24) * qJ(5) + (t51 * t56 + t23) * pkin(4) + t20) + t23 * mrSges(6,1) + t50 * t43;
t111 = m(5) * t20 + t23 * mrSges(5,1) + (t45 - t46) * t51 + (mrSges(5,2) - mrSges(6,3)) * t24 + t50 * t44 + t95;
t62 = (-mrSges(3,1) * t80 + mrSges(3,2) * t77) * qJD(1);
t48 = t57 * mrSges(4,1) + t58 * mrSges(4,2);
t7 = m(4) * t102 - t72 * mrSges(4,2) + t40 * mrSges(4,3) + t107 * t10 - t75 * t12 - t57 * t48 - t73 * t53;
t8 = m(4) * t92 + t72 * mrSges(4,1) - t41 * mrSges(4,3) - t58 * t48 + t73 * t52 - t111;
t4 = m(3) * (-t80 * g(3) - t106) - t63 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t62 * t99 + qJD(2) * t66 + t76 * t7 + t79 * t8;
t5 = m(3) * t93 - qJDD(2) * mrSges(3,2) + t64 * mrSges(3,3) - qJD(2) * t65 + t62 * t98 + t79 * t7 - t76 * t8;
t110 = t80 * t4 + t77 * t5;
t6 = m(2) * t94 + qJDD(1) * mrSges(2,1) - t82 * mrSges(2,2) - t112;
t1 = m(2) * t91 - t82 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t77 * t4 + t80 * t5;
t2 = [-m(1) * g(1) + t81 * t1 - t78 * t6, t1, t5, t7, t10, -t23 * mrSges(6,2) - t50 * t31 + t96; -m(1) * g(2) + t78 * t1 + t81 * t6, t6, t4, t8, t12, -t24 * mrSges(6,3) - t51 * t46 + t95; (-m(1) - m(2)) * g(3) + t110, -m(2) * g(3) + t110, t112, -t86, t111, -t38 * mrSges(6,1) + t24 * mrSges(6,2) + t51 * t31 - t56 * t43 + t109;];
f_new = t2;
