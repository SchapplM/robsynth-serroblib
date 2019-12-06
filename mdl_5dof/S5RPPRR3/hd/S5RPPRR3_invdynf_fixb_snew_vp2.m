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
% Datum: 2019-12-05 17:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:41:43
% EndTime: 2019-12-05 17:41:45
% DurationCPUTime: 1.15s
% Computational Cost: add. (12909->127), mult. (28559->166), div. (0->0), fcn. (19171->10), ass. (0->73)
t76 = qJD(1) ^ 2;
t68 = cos(pkin(9));
t63 = t68 ^ 2;
t66 = sin(pkin(9));
t94 = t66 ^ 2 + t63;
t100 = t94 * mrSges(4,3);
t99 = pkin(3) * t76;
t72 = sin(qJ(1));
t75 = cos(qJ(1));
t95 = t75 * g(2) + t72 * g(3);
t51 = qJDD(1) * pkin(1) + t95;
t88 = t72 * g(2) - t75 * g(3);
t52 = -t76 * pkin(1) + t88;
t67 = sin(pkin(8));
t69 = cos(pkin(8));
t97 = t67 * t51 + t69 * t52;
t35 = -t76 * pkin(2) + qJDD(1) * qJ(3) + t97;
t93 = pkin(6) * qJDD(1);
t65 = -g(1) + qJDD(2);
t91 = qJD(1) * qJD(3);
t96 = t68 * t65 - 0.2e1 * t66 * t91;
t25 = (t68 * t99 - t35 - t93) * t66 + t96;
t89 = t66 * t65 + (t35 + 0.2e1 * t91) * t68;
t26 = -t63 * t99 + t68 * t93 + t89;
t71 = sin(qJ(4));
t74 = cos(qJ(4));
t98 = t71 * t25 + t74 * t26;
t82 = -t66 * t71 + t68 * t74;
t45 = t82 * qJD(1);
t92 = t45 * qJD(4);
t83 = t66 * t74 + t68 * t71;
t40 = t83 * qJDD(1) + t92;
t46 = t83 * qJD(1);
t87 = t74 * t25 - t71 * t26;
t13 = (-t40 + t92) * pkin(7) + (t45 * t46 + qJDD(4)) * pkin(4) + t87;
t39 = -t46 * qJD(4) + t82 * qJDD(1);
t43 = qJD(4) * pkin(4) - t46 * pkin(7);
t44 = t45 ^ 2;
t14 = -t44 * pkin(4) + t39 * pkin(7) - qJD(4) * t43 + t98;
t70 = sin(qJ(5));
t73 = cos(qJ(5));
t33 = t73 * t45 - t70 * t46;
t19 = t33 * qJD(5) + t70 * t39 + t73 * t40;
t34 = t70 * t45 + t73 * t46;
t21 = -t33 * mrSges(6,1) + t34 * mrSges(6,2);
t64 = qJD(4) + qJD(5);
t29 = -t64 * mrSges(6,2) + t33 * mrSges(6,3);
t61 = qJDD(4) + qJDD(5);
t11 = m(6) * (t73 * t13 - t70 * t14) - t19 * mrSges(6,3) + t61 * mrSges(6,1) - t34 * t21 + t64 * t29;
t18 = -t34 * qJD(5) + t73 * t39 - t70 * t40;
t30 = t64 * mrSges(6,1) - t34 * mrSges(6,3);
t12 = m(6) * (t70 * t13 + t73 * t14) + t18 * mrSges(6,3) - t61 * mrSges(6,2) + t33 * t21 - t64 * t30;
t37 = -t45 * mrSges(5,1) + t46 * mrSges(5,2);
t41 = -qJD(4) * mrSges(5,2) + t45 * mrSges(5,3);
t8 = m(5) * t87 + qJDD(4) * mrSges(5,1) - t40 * mrSges(5,3) + qJD(4) * t41 + t73 * t11 + t70 * t12 - t46 * t37;
t84 = -t68 * mrSges(4,1) + t66 * mrSges(4,2);
t81 = qJDD(1) * mrSges(4,3) + t76 * t84;
t42 = qJD(4) * mrSges(5,1) - t46 * mrSges(5,3);
t9 = m(5) * t98 - qJDD(4) * mrSges(5,2) + t39 * mrSges(5,3) - qJD(4) * t42 - t70 * t11 + t73 * t12 + t45 * t37;
t6 = m(4) * t96 + t71 * t9 + t74 * t8 + (-m(4) * t35 - t81) * t66;
t7 = m(4) * t89 + t81 * t68 - t71 * t8 + t74 * t9;
t90 = m(3) * t65 + t68 * t6 + t66 * t7;
t86 = t69 * t51 - t67 * t52;
t85 = qJDD(3) - t86;
t78 = (-pkin(3) * t68 - pkin(2)) * qJDD(1) + (-t94 * pkin(6) - qJ(3)) * t76 + t85;
t80 = t18 * mrSges(6,1) + t33 * t29 - m(6) * (-t39 * pkin(4) - t44 * pkin(7) + t46 * t43 + t78) - t19 * mrSges(6,2) - t34 * t30;
t79 = -m(5) * t78 + t39 * mrSges(5,1) - t40 * mrSges(5,2) + t45 * t41 - t46 * t42 + t80;
t77 = m(4) * (-qJDD(1) * pkin(2) - t76 * qJ(3) + t85) - t79;
t10 = -t77 + (-mrSges(3,2) + t100) * t76 + (mrSges(3,1) - t84) * qJDD(1) + m(3) * t86;
t3 = m(3) * t97 - t76 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t66 * t6 + t68 * t7;
t2 = m(2) * t95 + qJDD(1) * mrSges(2,1) - t76 * mrSges(2,2) + t69 * t10 + t67 * t3;
t1 = m(2) * t88 - t76 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t67 * t10 + t69 * t3;
t4 = [(-m(1) - m(2)) * g(1) + t90, t1, t3, t7, t9, t12; -m(1) * g(2) - t72 * t1 - t75 * t2, t2, t10, t6, t8, t11; -m(1) * g(3) + t75 * t1 - t72 * t2, -m(2) * g(1) + t90, t90, t84 * qJDD(1) - t76 * t100 + t77, -t79, -t80;];
f_new = t4;
