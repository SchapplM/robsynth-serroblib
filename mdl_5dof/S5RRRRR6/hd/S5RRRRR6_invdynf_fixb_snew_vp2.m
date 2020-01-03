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
% Datum: 2020-01-03 12:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:14:52
% EndTime: 2020-01-03 12:14:55
% DurationCPUTime: 1.42s
% Computational Cost: add. (25010->142), mult. (32116->188), div. (0->0), fcn. (20711->10), ass. (0->78)
t64 = qJDD(1) + qJDD(2);
t70 = sin(qJ(3));
t75 = cos(qJ(3));
t66 = qJD(1) + qJD(2);
t90 = qJD(3) * t66;
t49 = t70 * t64 + t75 * t90;
t50 = t75 * t64 - t70 * t90;
t95 = t66 * t70;
t56 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t95;
t94 = t66 * t75;
t57 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t94;
t62 = t66 ^ 2;
t69 = sin(qJ(4));
t74 = cos(qJ(4));
t46 = (t69 * t75 + t70 * t74) * t66;
t30 = -t46 * qJD(4) - t69 * t49 + t74 * t50;
t45 = (-t69 * t70 + t74 * t75) * t66;
t31 = t45 * qJD(4) + t74 * t49 + t69 * t50;
t65 = qJD(3) + qJD(4);
t41 = -t65 * mrSges(5,2) + t45 * mrSges(5,3);
t42 = t65 * mrSges(5,1) - t46 * mrSges(5,3);
t58 = qJD(3) * pkin(3) - pkin(8) * t95;
t67 = t75 ^ 2;
t72 = sin(qJ(1));
t77 = cos(qJ(1));
t85 = -t77 * g(2) - t72 * g(3);
t54 = qJDD(1) * pkin(1) + t85;
t78 = qJD(1) ^ 2;
t88 = -t72 * g(2) + t77 * g(3);
t55 = -t78 * pkin(1) + t88;
t71 = sin(qJ(2));
t76 = cos(qJ(2));
t86 = t76 * t54 - t71 * t55;
t83 = -t64 * pkin(2) - t86;
t81 = -t50 * pkin(3) + t58 * t95 + (-pkin(8) * t67 - pkin(7)) * t62 + t83;
t68 = sin(qJ(5));
t73 = cos(qJ(5));
t36 = t68 * t45 + t73 * t46;
t18 = -t36 * qJD(5) + t73 * t30 - t68 * t31;
t35 = t73 * t45 - t68 * t46;
t19 = t35 * qJD(5) + t68 * t30 + t73 * t31;
t60 = qJD(5) + t65;
t32 = -t60 * mrSges(6,2) + t35 * mrSges(6,3);
t33 = t60 * mrSges(6,1) - t36 * mrSges(6,3);
t43 = t65 * pkin(4) - t46 * pkin(9);
t44 = t45 ^ 2;
t82 = t18 * mrSges(6,1) + t35 * t32 - m(6) * (-t30 * pkin(4) - t44 * pkin(9) + t46 * t43 + t81) - t19 * mrSges(6,2) - t36 * t33;
t80 = -m(5) * t81 + t30 * mrSges(5,1) - t31 * mrSges(5,2) + t45 * t41 - t46 * t42 + t82;
t99 = (t70 * t56 - t75 * t57) * t66 + m(4) * (-t62 * pkin(7) + t83) - t50 * mrSges(4,1) + t49 * mrSges(4,2) - t80;
t98 = -m(2) - m(3);
t48 = (-mrSges(4,1) * t75 + mrSges(4,2) * t70) * t66;
t63 = qJDD(3) + qJDD(4);
t91 = t71 * t54 + t76 * t55;
t40 = -t62 * pkin(2) + t64 * pkin(7) + t91;
t93 = t70 * t40;
t96 = pkin(3) * t62;
t25 = qJDD(3) * pkin(3) - t49 * pkin(8) - t93 + (pkin(8) * t90 + t70 * t96 - g(1)) * t75;
t89 = -t70 * g(1) + t75 * t40;
t26 = t50 * pkin(8) - qJD(3) * t58 - t67 * t96 + t89;
t87 = t74 * t25 - t69 * t26;
t13 = (t45 * t65 - t31) * pkin(9) + (t45 * t46 + t63) * pkin(4) + t87;
t92 = t69 * t25 + t74 * t26;
t14 = -t44 * pkin(4) + t30 * pkin(9) - t65 * t43 + t92;
t21 = -t35 * mrSges(6,1) + t36 * mrSges(6,2);
t59 = qJDD(5) + t63;
t11 = m(6) * (t73 * t13 - t68 * t14) - t19 * mrSges(6,3) + t59 * mrSges(6,1) - t36 * t21 + t60 * t32;
t12 = m(6) * (t68 * t13 + t73 * t14) + t18 * mrSges(6,3) - t59 * mrSges(6,2) + t35 * t21 - t60 * t33;
t38 = -t45 * mrSges(5,1) + t46 * mrSges(5,2);
t8 = m(5) * t87 + t63 * mrSges(5,1) - t31 * mrSges(5,3) + t73 * t11 + t68 * t12 - t46 * t38 + t65 * t41;
t9 = m(5) * t92 - t63 * mrSges(5,2) + t30 * mrSges(5,3) - t68 * t11 + t73 * t12 + t45 * t38 - t65 * t42;
t6 = m(4) * (-t75 * g(1) - t93) - t49 * mrSges(4,3) + qJDD(3) * mrSges(4,1) - t48 * t95 + qJD(3) * t57 + t69 * t9 + t74 * t8;
t7 = m(4) * t89 - qJDD(3) * mrSges(4,2) + t50 * mrSges(4,3) - qJD(3) * t56 + t48 * t94 - t69 * t8 + t74 * t9;
t97 = t75 * t6 + t70 * t7;
t10 = m(3) * t86 + t64 * mrSges(3,1) - t62 * mrSges(3,2) - t99;
t3 = m(3) * t91 - t62 * mrSges(3,1) - t64 * mrSges(3,2) - t70 * t6 + t75 * t7;
t2 = m(2) * t85 + qJDD(1) * mrSges(2,1) - t78 * mrSges(2,2) + t76 * t10 + t71 * t3;
t1 = m(2) * t88 - t78 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t71 * t10 + t76 * t3;
t4 = [(-m(1) + t98) * g(1) + t97, t1, t3, t7, t9, t12; -m(1) * g(2) + t72 * t1 + t77 * t2, t2, t10, t6, t8, t11; -m(1) * g(3) - t77 * t1 + t72 * t2, t98 * g(1) + t97, -m(3) * g(1) + t97, t99, -t80, -t82;];
f_new = t4;
