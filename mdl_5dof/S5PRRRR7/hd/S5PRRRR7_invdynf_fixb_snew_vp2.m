% Calculate vector of cutting forces with Newton-Euler
% S5PRRRR7
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRRRR7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR7_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR7_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR7_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR7_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR7_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:11:53
% EndTime: 2019-12-05 17:11:56
% DurationCPUTime: 1.26s
% Computational Cost: add. (14573->133), mult. (30139->180), div. (0->0), fcn. (20435->10), ass. (0->74)
t72 = sin(qJ(3));
t76 = cos(qJ(3));
t91 = qJD(2) * qJD(3);
t89 = t76 * t91;
t53 = t72 * qJDD(2) + t89;
t54 = t76 * qJDD(2) - t72 * t91;
t93 = qJD(2) * t72;
t57 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t93;
t92 = qJD(2) * t76;
t58 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t92;
t78 = qJD(2) ^ 2;
t71 = sin(qJ(4));
t75 = cos(qJ(4));
t46 = (t71 * t76 + t72 * t75) * qJD(2);
t30 = -t46 * qJD(4) - t71 * t53 + t75 * t54;
t45 = (-t71 * t72 + t75 * t76) * qJD(2);
t31 = t45 * qJD(4) + t75 * t53 + t71 * t54;
t66 = qJD(3) + qJD(4);
t41 = -t66 * mrSges(5,2) + t45 * mrSges(5,3);
t42 = t66 * mrSges(5,1) - t46 * mrSges(5,3);
t59 = qJD(3) * pkin(3) - pkin(7) * t93;
t67 = t76 ^ 2;
t69 = sin(pkin(9));
t94 = cos(pkin(9));
t56 = -t94 * g(1) - t69 * g(2);
t68 = -g(3) + qJDD(1);
t73 = sin(qJ(2));
t77 = cos(qJ(2));
t86 = -t73 * t56 + t77 * t68;
t83 = -qJDD(2) * pkin(2) - t86;
t81 = -t54 * pkin(3) + t59 * t93 + (-pkin(7) * t67 - pkin(6)) * t78 + t83;
t70 = sin(qJ(5));
t74 = cos(qJ(5));
t36 = t70 * t45 + t74 * t46;
t18 = -t36 * qJD(5) + t74 * t30 - t70 * t31;
t35 = t74 * t45 - t70 * t46;
t19 = t35 * qJD(5) + t70 * t30 + t74 * t31;
t62 = qJD(5) + t66;
t32 = -t62 * mrSges(6,2) + t35 * mrSges(6,3);
t33 = t62 * mrSges(6,1) - t36 * mrSges(6,3);
t43 = t66 * pkin(4) - t46 * pkin(8);
t44 = t45 ^ 2;
t82 = t18 * mrSges(6,1) + t35 * t32 - m(6) * (-t30 * pkin(4) - t44 * pkin(8) + t46 * t43 + t81) - t19 * mrSges(6,2) - t36 * t33;
t80 = -m(5) * t81 + t30 * mrSges(5,1) - t31 * mrSges(5,2) + t45 * t41 - t46 * t42 + t82;
t98 = (t72 * t57 - t76 * t58) * qJD(2) + m(4) * (-t78 * pkin(6) + t83) - t54 * mrSges(4,1) + t53 * mrSges(4,2) - t80;
t95 = t77 * t56 + t73 * t68;
t40 = -t78 * pkin(2) + qJDD(2) * pkin(6) + t95;
t55 = t69 * g(1) - t94 * g(2);
t87 = -t72 * t40 - t76 * t55;
t25 = (-t53 + t89) * pkin(7) + (t72 * t76 * t78 + qJDD(3)) * pkin(3) + t87;
t96 = t76 * t40 - t72 * t55;
t26 = -t67 * t78 * pkin(3) + t54 * pkin(7) - qJD(3) * t59 + t96;
t97 = t71 * t25 + t75 * t26;
t10 = m(3) * t86 + qJDD(2) * mrSges(3,1) - t78 * mrSges(3,2) - t98;
t52 = (-mrSges(4,1) * t76 + mrSges(4,2) * t72) * qJD(2);
t65 = qJDD(3) + qJDD(4);
t88 = t75 * t25 - t71 * t26;
t13 = (t45 * t66 - t31) * pkin(8) + (t45 * t46 + t65) * pkin(4) + t88;
t14 = -t44 * pkin(4) + t30 * pkin(8) - t66 * t43 + t97;
t21 = -t35 * mrSges(6,1) + t36 * mrSges(6,2);
t61 = qJDD(5) + t65;
t11 = m(6) * (t74 * t13 - t70 * t14) - t19 * mrSges(6,3) + t61 * mrSges(6,1) - t36 * t21 + t62 * t32;
t12 = m(6) * (t70 * t13 + t74 * t14) + t18 * mrSges(6,3) - t61 * mrSges(6,2) + t35 * t21 - t62 * t33;
t37 = -t45 * mrSges(5,1) + t46 * mrSges(5,2);
t7 = m(5) * t88 + t65 * mrSges(5,1) - t31 * mrSges(5,3) + t74 * t11 + t70 * t12 - t46 * t37 + t66 * t41;
t8 = m(5) * t97 - t65 * mrSges(5,2) + t30 * mrSges(5,3) - t70 * t11 + t74 * t12 + t45 * t37 - t66 * t42;
t5 = m(4) * t87 + qJDD(3) * mrSges(4,1) - t53 * mrSges(4,3) + qJD(3) * t58 - t52 * t93 + t75 * t7 + t71 * t8;
t6 = m(4) * t96 - qJDD(3) * mrSges(4,2) + t54 * mrSges(4,3) - qJD(3) * t57 + t52 * t92 - t71 * t7 + t75 * t8;
t3 = m(3) * t95 - t78 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t72 * t5 + t76 * t6;
t90 = m(2) * t68 + t77 * t10 + t73 * t3;
t85 = -t76 * t5 - t72 * t6;
t4 = (m(2) + m(3)) * t55 + t85;
t1 = m(2) * t56 - t73 * t10 + t77 * t3;
t2 = [-m(1) * g(1) + t94 * t1 - t69 * t4, t1, t3, t6, t8, t12; -m(1) * g(2) + t69 * t1 + t94 * t4, t4, t10, t5, t7, t11; -m(1) * g(3) + t90, t90, -m(3) * t55 - t85, t98, -t80, -t82;];
f_new = t2;
