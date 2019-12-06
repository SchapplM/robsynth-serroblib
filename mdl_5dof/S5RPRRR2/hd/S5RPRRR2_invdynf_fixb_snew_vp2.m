% Calculate vector of cutting forces with Newton-Euler
% S5RPRRR2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRRR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR2_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR2_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR2_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR2_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR2_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:11:24
% EndTime: 2019-12-05 18:11:32
% DurationCPUTime: 2.74s
% Computational Cost: add. (31292->155), mult. (77475->201), div. (0->0), fcn. (58396->10), ass. (0->85)
t86 = qJD(1) ^ 2;
t77 = cos(pkin(9));
t74 = t77 ^ 2;
t76 = sin(pkin(9));
t107 = t76 ^ 2 + t74;
t112 = t107 * mrSges(3,3);
t104 = qJD(1) * qJD(2);
t102 = -t77 * g(3) - 0.2e1 * t76 * t104;
t106 = pkin(6) * qJDD(1);
t110 = pkin(2) * t86;
t81 = sin(qJ(1));
t85 = cos(qJ(1));
t98 = -t85 * g(1) - t81 * g(2);
t65 = -t86 * pkin(1) + qJDD(1) * qJ(2) + t98;
t44 = (t77 * t110 - t106 - t65) * t76 + t102;
t99 = -t76 * g(3) + (0.2e1 * t104 + t65) * t77;
t45 = t77 * t106 - t74 * t110 + t99;
t80 = sin(qJ(3));
t84 = cos(qJ(3));
t100 = t84 * t44 - t80 * t45;
t94 = -t76 * t80 + t77 * t84;
t63 = t94 * qJD(1);
t95 = t76 * t84 + t77 * t80;
t64 = t95 * qJD(1);
t52 = -t63 * mrSges(4,1) + t64 * mrSges(4,2);
t105 = t63 * qJD(3);
t56 = t95 * qJDD(1) + t105;
t57 = -qJD(3) * mrSges(4,2) + t63 * mrSges(4,3);
t79 = sin(qJ(4));
t22 = (-t56 + t105) * pkin(7) + (t63 * t64 + qJDD(3)) * pkin(3) + t100;
t108 = t80 * t44 + t84 * t45;
t55 = -t64 * qJD(3) + t94 * qJDD(1);
t59 = qJD(3) * pkin(3) - t64 * pkin(7);
t62 = t63 ^ 2;
t26 = -t62 * pkin(3) + t55 * pkin(7) - qJD(3) * t59 + t108;
t83 = cos(qJ(4));
t101 = t83 * t22 - t79 * t26;
t47 = t83 * t63 - t79 * t64;
t31 = t47 * qJD(4) + t79 * t55 + t83 * t56;
t48 = t79 * t63 + t83 * t64;
t72 = qJDD(3) + qJDD(4);
t75 = qJD(3) + qJD(4);
t13 = (t47 * t75 - t31) * pkin(8) + (t47 * t48 + t72) * pkin(4) + t101;
t109 = t79 * t22 + t83 * t26;
t30 = -t48 * qJD(4) + t83 * t55 - t79 * t56;
t42 = t75 * pkin(4) - t48 * pkin(8);
t46 = t47 ^ 2;
t14 = -t46 * pkin(4) + t30 * pkin(8) - t75 * t42 + t109;
t78 = sin(qJ(5));
t82 = cos(qJ(5));
t35 = t82 * t47 - t78 * t48;
t19 = t35 * qJD(5) + t78 * t30 + t82 * t31;
t36 = t78 * t47 + t82 * t48;
t24 = -t35 * mrSges(6,1) + t36 * mrSges(6,2);
t70 = qJD(5) + t75;
t32 = -t70 * mrSges(6,2) + t35 * mrSges(6,3);
t69 = qJDD(5) + t72;
t11 = m(6) * (t82 * t13 - t78 * t14) - t19 * mrSges(6,3) + t69 * mrSges(6,1) - t36 * t24 + t70 * t32;
t18 = -t36 * qJD(5) + t82 * t30 - t78 * t31;
t33 = t70 * mrSges(6,1) - t36 * mrSges(6,3);
t12 = m(6) * (t78 * t13 + t82 * t14) + t18 * mrSges(6,3) - t69 * mrSges(6,2) + t35 * t24 - t70 * t33;
t37 = -t47 * mrSges(5,1) + t48 * mrSges(5,2);
t40 = -t75 * mrSges(5,2) + t47 * mrSges(5,3);
t8 = m(5) * t101 + t72 * mrSges(5,1) - t31 * mrSges(5,3) + t82 * t11 + t78 * t12 - t48 * t37 + t75 * t40;
t41 = t75 * mrSges(5,1) - t48 * mrSges(5,3);
t9 = m(5) * t109 - t72 * mrSges(5,2) + t30 * mrSges(5,3) - t78 * t11 + t82 * t12 + t47 * t37 - t75 * t41;
t6 = m(4) * t100 + qJDD(3) * mrSges(4,1) - t56 * mrSges(4,3) + qJD(3) * t57 - t64 * t52 + t79 * t9 + t83 * t8;
t58 = qJD(3) * mrSges(4,1) - t64 * mrSges(4,3);
t7 = m(4) * t108 - qJDD(3) * mrSges(4,2) + t55 * mrSges(4,3) - qJD(3) * t58 + t63 * t52 - t79 * t8 + t83 * t9;
t96 = -t77 * mrSges(3,1) + t76 * mrSges(3,2);
t93 = qJDD(1) * mrSges(3,3) + t86 * t96;
t4 = m(3) * t102 + t80 * t7 + t84 * t6 + (-m(3) * t65 - t93) * t76;
t5 = m(3) * t99 - t80 * t6 + t84 * t7 + t93 * t77;
t111 = t77 * t4 + t76 * t5;
t103 = t81 * g(1) - t85 * g(2);
t97 = qJDD(2) - t103;
t90 = (-pkin(2) * t77 - pkin(1)) * qJDD(1) + (-t107 * pkin(6) - qJ(2)) * t86 + t97;
t89 = -t55 * pkin(3) - t62 * pkin(7) + t64 * t59 + t90;
t92 = t18 * mrSges(6,1) + t35 * t32 - m(6) * (-t30 * pkin(4) - t46 * pkin(8) + t48 * t42 + t89) - t19 * mrSges(6,2) - t36 * t33;
t91 = -m(5) * t89 + t30 * mrSges(5,1) - t31 * mrSges(5,2) + t47 * t40 - t48 * t41 + t92;
t88 = -m(4) * t90 + t55 * mrSges(4,1) - t56 * mrSges(4,2) + t63 * t57 - t64 * t58 + t91;
t87 = m(3) * (-qJDD(1) * pkin(1) - t86 * qJ(2) + t97) - t88;
t10 = (-mrSges(2,2) + t112) * t86 + (mrSges(2,1) - t96) * qJDD(1) + m(2) * t103 - t87;
t1 = m(2) * t98 - t86 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t76 * t4 + t77 * t5;
t2 = [-m(1) * g(1) + t85 * t1 - t81 * t10, t1, t5, t7, t9, t12; -m(1) * g(2) + t81 * t1 + t85 * t10, t10, t4, t6, t8, t11; (-m(1) - m(2)) * g(3) + t111, -m(2) * g(3) + t111, t96 * qJDD(1) - t86 * t112 + t87, -t88, -t91, -t92;];
f_new = t2;
