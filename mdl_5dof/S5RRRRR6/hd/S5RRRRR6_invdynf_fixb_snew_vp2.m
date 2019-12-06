% Calculate vector of cutting forces with Newton-Euler
% S5RRRRR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-05 19:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRRR6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR6_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR6_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR6_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 19:00:09
% EndTime: 2019-12-05 19:00:12
% DurationCPUTime: 1.43s
% Computational Cost: add. (25010->142), mult. (32116->188), div. (0->0), fcn. (20711->10), ass. (0->78)
t66 = qJDD(1) + qJDD(2);
t72 = sin(qJ(3));
t77 = cos(qJ(3));
t68 = qJD(1) + qJD(2);
t91 = qJD(3) * t68;
t49 = t72 * t66 + t77 * t91;
t50 = t77 * t66 - t72 * t91;
t97 = t68 * t72;
t56 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t97;
t96 = t68 * t77;
t57 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t96;
t64 = t68 ^ 2;
t71 = sin(qJ(4));
t76 = cos(qJ(4));
t46 = (t71 * t77 + t72 * t76) * t68;
t30 = -t46 * qJD(4) - t71 * t49 + t76 * t50;
t45 = (-t71 * t72 + t76 * t77) * t68;
t31 = t45 * qJD(4) + t76 * t49 + t71 * t50;
t67 = qJD(3) + qJD(4);
t41 = -t67 * mrSges(5,2) + t45 * mrSges(5,3);
t42 = t67 * mrSges(5,1) - t46 * mrSges(5,3);
t58 = qJD(3) * pkin(3) - pkin(8) * t97;
t69 = t77 ^ 2;
t74 = sin(qJ(1));
t79 = cos(qJ(1));
t92 = t79 * g(2) + t74 * g(3);
t54 = qJDD(1) * pkin(1) + t92;
t80 = qJD(1) ^ 2;
t89 = t74 * g(2) - t79 * g(3);
t55 = -t80 * pkin(1) + t89;
t73 = sin(qJ(2));
t78 = cos(qJ(2));
t87 = t78 * t54 - t73 * t55;
t85 = -t66 * pkin(2) - t87;
t83 = -t50 * pkin(3) + t58 * t97 + (-pkin(8) * t69 - pkin(7)) * t64 + t85;
t70 = sin(qJ(5));
t75 = cos(qJ(5));
t36 = t70 * t45 + t75 * t46;
t18 = -t36 * qJD(5) + t75 * t30 - t70 * t31;
t35 = t75 * t45 - t70 * t46;
t19 = t35 * qJD(5) + t70 * t30 + t75 * t31;
t60 = qJD(5) + t67;
t32 = -t60 * mrSges(6,2) + t35 * mrSges(6,3);
t33 = t60 * mrSges(6,1) - t36 * mrSges(6,3);
t43 = t67 * pkin(4) - t46 * pkin(9);
t44 = t45 ^ 2;
t84 = t18 * mrSges(6,1) + t35 * t32 - m(6) * (-t30 * pkin(4) - t44 * pkin(9) + t46 * t43 + t83) - t19 * mrSges(6,2) - t36 * t33;
t82 = -m(5) * t83 + t30 * mrSges(5,1) - t31 * mrSges(5,2) + t45 * t41 - t46 * t42 + t84;
t101 = (t72 * t56 - t77 * t57) * t68 + m(4) * (-t64 * pkin(7) + t85) - t50 * mrSges(4,1) + t49 * mrSges(4,2) - t82;
t100 = -m(2) - m(3);
t48 = (-mrSges(4,1) * t77 + mrSges(4,2) * t72) * t68;
t65 = qJDD(3) + qJDD(4);
t93 = t73 * t54 + t78 * t55;
t40 = -t64 * pkin(2) + t66 * pkin(7) + t93;
t95 = t72 * t40;
t98 = pkin(3) * t64;
t25 = qJDD(3) * pkin(3) - t49 * pkin(8) - t95 + (pkin(8) * t91 + t72 * t98 - g(1)) * t77;
t90 = -t72 * g(1) + t77 * t40;
t26 = t50 * pkin(8) - qJD(3) * t58 - t69 * t98 + t90;
t88 = t76 * t25 - t71 * t26;
t13 = (t45 * t67 - t31) * pkin(9) + (t45 * t46 + t65) * pkin(4) + t88;
t94 = t71 * t25 + t76 * t26;
t14 = -t44 * pkin(4) + t30 * pkin(9) - t67 * t43 + t94;
t21 = -t35 * mrSges(6,1) + t36 * mrSges(6,2);
t59 = qJDD(5) + t65;
t11 = m(6) * (t75 * t13 - t70 * t14) - t19 * mrSges(6,3) + t59 * mrSges(6,1) - t36 * t21 + t60 * t32;
t12 = m(6) * (t70 * t13 + t75 * t14) + t18 * mrSges(6,3) - t59 * mrSges(6,2) + t35 * t21 - t60 * t33;
t38 = -t45 * mrSges(5,1) + t46 * mrSges(5,2);
t8 = m(5) * t88 + t65 * mrSges(5,1) - t31 * mrSges(5,3) + t75 * t11 + t70 * t12 - t46 * t38 + t67 * t41;
t9 = m(5) * t94 - t65 * mrSges(5,2) + t30 * mrSges(5,3) - t70 * t11 + t75 * t12 + t45 * t38 - t67 * t42;
t6 = m(4) * (-t77 * g(1) - t95) - t49 * mrSges(4,3) + qJDD(3) * mrSges(4,1) - t48 * t97 + qJD(3) * t57 + t71 * t9 + t76 * t8;
t7 = m(4) * t90 - qJDD(3) * mrSges(4,2) + t50 * mrSges(4,3) - qJD(3) * t56 + t48 * t96 - t71 * t8 + t76 * t9;
t99 = t77 * t6 + t72 * t7;
t10 = m(3) * t87 + t66 * mrSges(3,1) - t64 * mrSges(3,2) - t101;
t3 = m(3) * t93 - t64 * mrSges(3,1) - t66 * mrSges(3,2) - t72 * t6 + t77 * t7;
t2 = m(2) * t92 + qJDD(1) * mrSges(2,1) - t80 * mrSges(2,2) + t78 * t10 + t73 * t3;
t1 = m(2) * t89 - t80 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t73 * t10 + t78 * t3;
t4 = [(-m(1) + t100) * g(1) + t99, t1, t3, t7, t9, t12; -m(1) * g(2) - t74 * t1 - t79 * t2, t2, t10, t6, t8, t11; -m(1) * g(3) + t79 * t1 - t74 * t2, t100 * g(1) + t99, -m(3) * g(1) + t99, t101, -t82, -t84;];
f_new = t4;
