% Calculate vector of cutting forces with Newton-Euler
% S5RRRRP6
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
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRRP6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP6_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP6_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP6_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP6_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP6_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:52:52
% EndTime: 2019-12-31 21:52:56
% DurationCPUTime: 1.19s
% Computational Cost: add. (12578->163), mult. (25429->208), div. (0->0), fcn. (17086->8), ass. (0->77)
t77 = sin(qJ(3));
t78 = sin(qJ(2));
t81 = cos(qJ(3));
t82 = cos(qJ(2));
t60 = (t77 * t78 - t81 * t82) * qJD(1);
t100 = qJD(1) * qJD(2);
t66 = t78 * qJDD(1) + t82 * t100;
t67 = t82 * qJDD(1) - t78 * t100;
t102 = qJD(1) * t78;
t68 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t102;
t101 = qJD(1) * t82;
t69 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t101;
t84 = qJD(1) ^ 2;
t61 = (t77 * t82 + t78 * t81) * qJD(1);
t42 = -t61 * qJD(3) - t77 * t66 + t81 * t67;
t43 = -t60 * qJD(3) + t81 * t66 + t77 * t67;
t74 = qJD(2) + qJD(3);
t70 = qJD(2) * pkin(2) - pkin(7) * t102;
t75 = t82 ^ 2;
t79 = sin(qJ(1));
t83 = cos(qJ(1));
t96 = t79 * g(1) - t83 * g(2);
t89 = -qJDD(1) * pkin(1) - t96;
t86 = -t67 * pkin(2) + t70 * t102 + (-pkin(7) * t75 - pkin(6)) * t84 + t89;
t19 = (t60 * t74 - t43) * pkin(8) + (t61 * t74 - t42) * pkin(3) + t86;
t92 = -t83 * g(1) - t79 * g(2);
t63 = -t84 * pkin(1) + qJDD(1) * pkin(6) + t92;
t107 = t78 * t63;
t108 = pkin(2) * t84;
t35 = qJDD(2) * pkin(2) - t66 * pkin(7) - t107 + (pkin(7) * t100 + t78 * t108 - g(3)) * t82;
t95 = -t78 * g(3) + t82 * t63;
t36 = t67 * pkin(7) - qJD(2) * t70 - t75 * t108 + t95;
t104 = t77 * t35 + t81 * t36;
t52 = t60 * pkin(3) - t61 * pkin(8);
t72 = t74 ^ 2;
t73 = qJDD(2) + qJDD(3);
t22 = -t72 * pkin(3) + t73 * pkin(8) - t60 * t52 + t104;
t76 = sin(qJ(4));
t80 = cos(qJ(4));
t105 = t76 * t19 + t80 * t22;
t55 = t80 * t61 + t76 * t74;
t25 = -t55 * qJD(4) - t76 * t43 + t80 * t73;
t54 = -t76 * t61 + t80 * t74;
t34 = -t54 * mrSges(5,1) + t55 * mrSges(5,2);
t40 = qJDD(4) - t42;
t59 = qJD(4) + t60;
t48 = t59 * mrSges(6,1) - t55 * mrSges(6,3);
t49 = t59 * mrSges(5,1) - t55 * mrSges(5,3);
t33 = -t54 * mrSges(6,1) + t55 * mrSges(6,2);
t47 = t59 * pkin(4) - t55 * qJ(5);
t53 = t54 ^ 2;
t98 = m(6) * (-t53 * pkin(4) + t25 * qJ(5) + 0.2e1 * qJD(5) * t54 - t59 * t47 + t105) + t25 * mrSges(6,3) + t54 * t33;
t11 = m(5) * t105 + t25 * mrSges(5,3) + t54 * t34 + (-t49 - t48) * t59 + (-mrSges(5,2) - mrSges(6,2)) * t40 + t98;
t56 = -t74 * mrSges(4,2) - t60 * mrSges(4,3);
t57 = t74 * mrSges(4,1) - t61 * mrSges(4,3);
t26 = t54 * qJD(4) + t80 * t43 + t76 * t73;
t46 = -t59 * mrSges(5,2) + t54 * mrSges(5,3);
t94 = t80 * t19 - t76 * t22;
t45 = -t59 * mrSges(6,2) + t54 * mrSges(6,3);
t99 = m(6) * (-0.2e1 * qJD(5) * t55 + (t54 * t59 - t26) * qJ(5) + (t54 * t55 + t40) * pkin(4) + t94) + t59 * t45 + t40 * mrSges(6,1);
t9 = m(5) * t94 + t40 * mrSges(5,1) + t59 * t46 + (-t34 - t33) * t55 + (-mrSges(5,3) - mrSges(6,3)) * t26 + t99;
t88 = -m(4) * t86 + t42 * mrSges(4,1) - t43 * mrSges(4,2) - t76 * t11 - t60 * t56 - t61 * t57 - t80 * t9;
t111 = (t78 * t68 - t82 * t69) * qJD(1) + m(3) * (-t84 * pkin(6) + t89) - t67 * mrSges(3,1) + t66 * mrSges(3,2) - t88;
t93 = t81 * t35 - t77 * t36;
t21 = -t73 * pkin(3) - t72 * pkin(8) + t61 * t52 - t93;
t97 = m(6) * (-t25 * pkin(4) - t53 * qJ(5) + t55 * t47 + qJDD(5) + t21) + t26 * mrSges(6,2) + t55 * t48;
t110 = m(5) * t21 + t26 * mrSges(5,2) - (t46 + t45) * t54 - (mrSges(5,1) + mrSges(6,1)) * t25 + t55 * t49 + t97;
t51 = t60 * mrSges(4,1) + t61 * mrSges(4,2);
t12 = m(4) * t93 + t73 * mrSges(4,1) - t43 * mrSges(4,3) - t61 * t51 + t74 * t56 - t110;
t65 = (-mrSges(3,1) * t82 + mrSges(3,2) * t78) * qJD(1);
t7 = m(4) * t104 - t73 * mrSges(4,2) + t42 * mrSges(4,3) + t80 * t11 - t60 * t51 - t74 * t57 - t76 * t9;
t4 = m(3) * (-t82 * g(3) - t107) - t66 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t65 * t102 + qJD(2) * t69 + t77 * t7 + t81 * t12;
t5 = m(3) * t95 - qJDD(2) * mrSges(3,2) + t67 * mrSges(3,3) - qJD(2) * t68 + t65 * t101 - t77 * t12 + t81 * t7;
t109 = t82 * t4 + t78 * t5;
t6 = m(2) * t96 + qJDD(1) * mrSges(2,1) - t84 * mrSges(2,2) - t111;
t1 = m(2) * t92 - t84 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t78 * t4 + t82 * t5;
t2 = [-m(1) * g(1) + t83 * t1 - t79 * t6, t1, t5, t7, t11, -t40 * mrSges(6,2) - t59 * t48 + t98; -m(1) * g(2) + t79 * t1 + t83 * t6, t6, t4, t12, t9, -t26 * mrSges(6,3) - t55 * t33 + t99; (-m(1) - m(2)) * g(3) + t109, -m(2) * g(3) + t109, t111, -t88, t110, -t25 * mrSges(6,1) - t54 * t45 + t97;];
f_new = t2;
