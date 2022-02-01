% Calculate vector of cutting forces with Newton-Euler
% S5RPPRR3
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% m [6x1]
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
% Datum: 2022-01-23 09:15
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPPRR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR3_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR3_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR3_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:14:12
% EndTime: 2022-01-23 09:14:14
% DurationCPUTime: 1.18s
% Computational Cost: add. (12909->127), mult. (28559->166), div. (0->0), fcn. (19171->10), ass. (0->73)
t74 = qJD(1) ^ 2;
t66 = cos(pkin(9));
t61 = t66 ^ 2;
t64 = sin(pkin(9));
t93 = t64 ^ 2 + t61;
t98 = t93 * mrSges(4,3);
t97 = pkin(3) * t74;
t70 = sin(qJ(1));
t73 = cos(qJ(1));
t87 = t70 * g(1) - t73 * g(2);
t51 = qJDD(1) * pkin(1) + t87;
t84 = -t73 * g(1) - t70 * g(2);
t52 = -t74 * pkin(1) + t84;
t65 = sin(pkin(8));
t67 = cos(pkin(8));
t95 = t65 * t51 + t67 * t52;
t35 = -t74 * pkin(2) + qJDD(1) * qJ(3) + t95;
t92 = pkin(6) * qJDD(1);
t63 = -g(3) + qJDD(2);
t90 = qJD(1) * qJD(3);
t94 = t66 * t63 - 0.2e1 * t64 * t90;
t25 = (t66 * t97 - t35 - t92) * t64 + t94;
t88 = t64 * t63 + (t35 + 0.2e1 * t90) * t66;
t26 = -t61 * t97 + t66 * t92 + t88;
t69 = sin(qJ(4));
t72 = cos(qJ(4));
t96 = t69 * t25 + t72 * t26;
t80 = -t64 * t69 + t66 * t72;
t45 = t80 * qJD(1);
t91 = t45 * qJD(4);
t82 = -t66 * mrSges(4,1) + t64 * mrSges(4,2);
t79 = qJDD(1) * mrSges(4,3) + t74 * t82;
t81 = t64 * t72 + t66 * t69;
t40 = t81 * qJDD(1) + t91;
t46 = t81 * qJD(1);
t86 = t72 * t25 - t69 * t26;
t13 = (-t40 + t91) * pkin(7) + (t45 * t46 + qJDD(4)) * pkin(4) + t86;
t39 = -t46 * qJD(4) + t80 * qJDD(1);
t43 = qJD(4) * pkin(4) - t46 * pkin(7);
t44 = t45 ^ 2;
t14 = -t44 * pkin(4) + t39 * pkin(7) - qJD(4) * t43 + t96;
t68 = sin(qJ(5));
t71 = cos(qJ(5));
t33 = t71 * t45 - t68 * t46;
t19 = t33 * qJD(5) + t68 * t39 + t71 * t40;
t34 = t68 * t45 + t71 * t46;
t21 = -t33 * mrSges(6,1) + t34 * mrSges(6,2);
t62 = qJD(4) + qJD(5);
t29 = -t62 * mrSges(6,2) + t33 * mrSges(6,3);
t59 = qJDD(4) + qJDD(5);
t11 = m(6) * (t71 * t13 - t68 * t14) - t19 * mrSges(6,3) + t59 * mrSges(6,1) - t34 * t21 + t62 * t29;
t18 = -t34 * qJD(5) + t71 * t39 - t68 * t40;
t30 = t62 * mrSges(6,1) - t34 * mrSges(6,3);
t12 = m(6) * (t68 * t13 + t71 * t14) + t18 * mrSges(6,3) - t59 * mrSges(6,2) + t33 * t21 - t62 * t30;
t37 = -t45 * mrSges(5,1) + t46 * mrSges(5,2);
t41 = -qJD(4) * mrSges(5,2) + t45 * mrSges(5,3);
t8 = m(5) * t86 + qJDD(4) * mrSges(5,1) - t40 * mrSges(5,3) + qJD(4) * t41 + t71 * t11 + t68 * t12 - t46 * t37;
t42 = qJD(4) * mrSges(5,1) - t46 * mrSges(5,3);
t9 = m(5) * t96 - qJDD(4) * mrSges(5,2) + t39 * mrSges(5,3) - qJD(4) * t42 - t68 * t11 + t71 * t12 + t45 * t37;
t6 = m(4) * t94 + t69 * t9 + t72 * t8 + (-m(4) * t35 - t79) * t64;
t7 = m(4) * t88 + t79 * t66 - t69 * t8 + t72 * t9;
t89 = m(3) * t63 + t66 * t6 + t64 * t7;
t85 = t67 * t51 - t65 * t52;
t83 = qJDD(3) - t85;
t76 = (-pkin(3) * t66 - pkin(2)) * qJDD(1) + (-t93 * pkin(6) - qJ(3)) * t74 + t83;
t78 = t18 * mrSges(6,1) + t33 * t29 - m(6) * (-t39 * pkin(4) - t44 * pkin(7) + t46 * t43 + t76) - t19 * mrSges(6,2) - t34 * t30;
t77 = -m(5) * t76 + t39 * mrSges(5,1) - t40 * mrSges(5,2) + t45 * t41 - t46 * t42 + t78;
t75 = m(4) * (-qJDD(1) * pkin(2) - t74 * qJ(3) + t83) - t77;
t10 = (mrSges(3,1) - t82) * qJDD(1) - t75 + m(3) * t85 + (-mrSges(3,2) + t98) * t74;
t3 = m(3) * t95 - t74 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t64 * t6 + t66 * t7;
t2 = m(2) * t84 - t74 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t65 * t10 + t67 * t3;
t1 = m(2) * t87 + qJDD(1) * mrSges(2,1) - t74 * mrSges(2,2) + t67 * t10 + t65 * t3;
t4 = [-m(1) * g(1) - t70 * t1 + t73 * t2, t2, t3, t7, t9, t12; -m(1) * g(2) + t73 * t1 + t70 * t2, t1, t10, t6, t8, t11; (-m(1) - m(2)) * g(3) + t89, -m(2) * g(3) + t89, t89, t82 * qJDD(1) - t74 * t98 + t75, -t77, -t78;];
f_new = t4;
