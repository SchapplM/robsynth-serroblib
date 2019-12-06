% Calculate vector of cutting forces with Newton-Euler
% S5PRPRR1
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRPRR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR1_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR1_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR1_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR1_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR1_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:42:44
% EndTime: 2019-12-05 15:42:46
% DurationCPUTime: 1.09s
% Computational Cost: add. (11837->120), mult. (26948->160), div. (0->0), fcn. (19171->10), ass. (0->72)
t74 = qJD(2) ^ 2;
t66 = cos(pkin(9));
t61 = t66 ^ 2;
t64 = sin(pkin(9));
t92 = t64 ^ 2 + t61;
t97 = t92 * mrSges(4,3);
t96 = pkin(3) * t74;
t65 = sin(pkin(8));
t67 = cos(pkin(8));
t51 = t65 * g(1) - t67 * g(2);
t52 = -t67 * g(1) - t65 * g(2);
t70 = sin(qJ(2));
t73 = cos(qJ(2));
t94 = t70 * t51 + t73 * t52;
t40 = -t74 * pkin(2) + qJDD(2) * qJ(3) + t94;
t91 = pkin(6) * qJDD(2);
t63 = -g(3) + qJDD(1);
t89 = qJD(2) * qJD(3);
t93 = t66 * t63 - 0.2e1 * t64 * t89;
t25 = (t66 * t96 - t40 - t91) * t64 + t93;
t87 = t64 * t63 + (t40 + 0.2e1 * t89) * t66;
t26 = -t61 * t96 + t66 * t91 + t87;
t69 = sin(qJ(4));
t72 = cos(qJ(4));
t95 = t69 * t25 + t72 * t26;
t80 = -t64 * t69 + t66 * t72;
t45 = t80 * qJD(2);
t90 = t45 * qJD(4);
t82 = -t66 * mrSges(4,1) + t64 * mrSges(4,2);
t79 = qJDD(2) * mrSges(4,3) + t74 * t82;
t81 = t64 * t72 + t66 * t69;
t39 = t81 * qJDD(2) + t90;
t46 = t81 * qJD(2);
t85 = t72 * t25 - t69 * t26;
t13 = (-t39 + t90) * pkin(7) + (t45 * t46 + qJDD(4)) * pkin(4) + t85;
t38 = -t46 * qJD(4) + t80 * qJDD(2);
t43 = qJD(4) * pkin(4) - t46 * pkin(7);
t44 = t45 ^ 2;
t14 = -t44 * pkin(4) + t38 * pkin(7) - qJD(4) * t43 + t95;
t68 = sin(qJ(5));
t71 = cos(qJ(5));
t31 = t71 * t45 - t68 * t46;
t19 = t31 * qJD(5) + t68 * t38 + t71 * t39;
t32 = t68 * t45 + t71 * t46;
t21 = -t31 * mrSges(6,1) + t32 * mrSges(6,2);
t62 = qJD(4) + qJD(5);
t29 = -t62 * mrSges(6,2) + t31 * mrSges(6,3);
t59 = qJDD(4) + qJDD(5);
t11 = m(6) * (t71 * t13 - t68 * t14) - t19 * mrSges(6,3) + t59 * mrSges(6,1) - t32 * t21 + t62 * t29;
t18 = -t32 * qJD(5) + t71 * t38 - t68 * t39;
t30 = t62 * mrSges(6,1) - t32 * mrSges(6,3);
t12 = m(6) * (t68 * t13 + t71 * t14) + t18 * mrSges(6,3) - t59 * mrSges(6,2) + t31 * t21 - t62 * t30;
t35 = -t45 * mrSges(5,1) + t46 * mrSges(5,2);
t41 = -qJD(4) * mrSges(5,2) + t45 * mrSges(5,3);
t8 = m(5) * t85 + qJDD(4) * mrSges(5,1) - t39 * mrSges(5,3) + qJD(4) * t41 + t71 * t11 + t68 * t12 - t46 * t35;
t42 = qJD(4) * mrSges(5,1) - t46 * mrSges(5,3);
t9 = m(5) * t95 - qJDD(4) * mrSges(5,2) + t38 * mrSges(5,3) - qJD(4) * t42 - t68 * t11 + t71 * t12 + t45 * t35;
t6 = m(4) * t93 + t69 * t9 + t72 * t8 + (-m(4) * t40 - t79) * t64;
t7 = m(4) * t87 + t79 * t66 - t69 * t8 + t72 * t9;
t88 = m(3) * t63 + t66 * t6 + t64 * t7;
t86 = m(2) * t63 + t88;
t84 = t73 * t51 - t70 * t52;
t83 = qJDD(3) - t84;
t76 = (-pkin(3) * t66 - pkin(2)) * qJDD(2) + (-t92 * pkin(6) - qJ(3)) * t74 + t83;
t78 = t18 * mrSges(6,1) + t31 * t29 - m(6) * (-t38 * pkin(4) - t44 * pkin(7) + t46 * t43 + t76) - t19 * mrSges(6,2) - t32 * t30;
t77 = -m(5) * t76 + t38 * mrSges(5,1) - t39 * mrSges(5,2) + t45 * t41 - t46 * t42 + t78;
t75 = m(4) * (-qJDD(2) * pkin(2) - t74 * qJ(3) + t83) - t77;
t10 = (-mrSges(3,2) + t97) * t74 - t75 + (mrSges(3,1) - t82) * qJDD(2) + m(3) * t84;
t3 = m(3) * t94 - t74 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t64 * t6 + t66 * t7;
t2 = m(2) * t52 - t70 * t10 + t73 * t3;
t1 = m(2) * t51 + t73 * t10 + t70 * t3;
t4 = [-m(1) * g(1) - t65 * t1 + t67 * t2, t2, t3, t7, t9, t12; -m(1) * g(2) + t67 * t1 + t65 * t2, t1, t10, t6, t8, t11; -m(1) * g(3) + t86, t86, t88, t82 * qJDD(2) - t74 * t97 + t75, -t77, -t78;];
f_new = t4;
