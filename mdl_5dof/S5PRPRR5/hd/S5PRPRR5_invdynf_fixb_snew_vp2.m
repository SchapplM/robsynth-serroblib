% Calculate vector of cutting forces with Newton-Euler
% S5PRPRR5
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
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRPRR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR5_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR5_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR5_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR5_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR5_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:53:32
% EndTime: 2019-12-05 15:53:36
% DurationCPUTime: 1.10s
% Computational Cost: add. (11728->121), mult. (26622->160), div. (0->0), fcn. (18917->10), ass. (0->72)
t73 = qJD(2) ^ 2;
t66 = cos(pkin(9));
t61 = t66 ^ 2;
t64 = sin(pkin(9));
t92 = t64 ^ 2 + t61;
t97 = t92 * mrSges(4,3);
t96 = pkin(3) * t73;
t65 = sin(pkin(8));
t91 = cos(pkin(8));
t53 = -t91 * g(1) - t65 * g(2);
t63 = -g(3) + qJDD(1);
t69 = sin(qJ(2));
t72 = cos(qJ(2));
t93 = t72 * t53 + t69 * t63;
t40 = -t73 * pkin(2) + qJDD(2) * qJ(3) + t93;
t90 = pkin(6) * qJDD(2);
t52 = t65 * g(1) - t91 * g(2);
t88 = qJD(2) * qJD(3);
t94 = -t66 * t52 - 0.2e1 * t64 * t88;
t25 = (t66 * t96 - t40 - t90) * t64 + t94;
t86 = -t64 * t52 + (t40 + 0.2e1 * t88) * t66;
t26 = -t61 * t96 + t66 * t90 + t86;
t68 = sin(qJ(4));
t71 = cos(qJ(4));
t95 = t68 * t25 + t71 * t26;
t79 = -t64 * t68 + t66 * t71;
t45 = t79 * qJD(2);
t89 = t45 * qJD(4);
t80 = t64 * t71 + t66 * t68;
t46 = t80 * qJD(2);
t36 = -t46 * qJD(4) + t79 * qJDD(2);
t37 = t80 * qJDD(2) + t89;
t41 = -qJD(4) * mrSges(5,2) + t45 * mrSges(5,3);
t42 = qJD(4) * mrSges(5,1) - t46 * mrSges(5,3);
t84 = -t69 * t53 + t72 * t63;
t82 = qJDD(3) - t84;
t75 = (-pkin(3) * t66 - pkin(2)) * qJDD(2) + (-t92 * pkin(6) - qJ(3)) * t73 + t82;
t67 = sin(qJ(5));
t70 = cos(qJ(5));
t32 = t67 * t45 + t70 * t46;
t18 = -t32 * qJD(5) + t70 * t36 - t67 * t37;
t31 = t70 * t45 - t67 * t46;
t19 = t31 * qJD(5) + t67 * t36 + t70 * t37;
t62 = qJD(4) + qJD(5);
t28 = -t62 * mrSges(6,2) + t31 * mrSges(6,3);
t29 = t62 * mrSges(6,1) - t32 * mrSges(6,3);
t43 = qJD(4) * pkin(4) - t46 * pkin(7);
t44 = t45 ^ 2;
t77 = t18 * mrSges(6,1) + t31 * t28 - m(6) * (-t36 * pkin(4) - t44 * pkin(7) + t46 * t43 + t75) - t19 * mrSges(6,2) - t32 * t29;
t76 = -m(5) * t75 + t36 * mrSges(5,1) - t37 * mrSges(5,2) + t45 * t41 - t46 * t42 + t77;
t74 = m(4) * (-qJDD(2) * pkin(2) - t73 * qJ(3) + t82) - t76;
t81 = -t66 * mrSges(4,1) + t64 * mrSges(4,2);
t10 = (-mrSges(3,2) + t97) * t73 + (mrSges(3,1) - t81) * qJDD(2) + m(3) * t84 - t74;
t85 = t71 * t25 - t68 * t26;
t13 = (-t37 + t89) * pkin(7) + (t45 * t46 + qJDD(4)) * pkin(4) + t85;
t14 = -t44 * pkin(4) + t36 * pkin(7) - qJD(4) * t43 + t95;
t21 = -t31 * mrSges(6,1) + t32 * mrSges(6,2);
t59 = qJDD(4) + qJDD(5);
t11 = m(6) * (t70 * t13 - t67 * t14) - t19 * mrSges(6,3) + t59 * mrSges(6,1) - t32 * t21 + t62 * t28;
t12 = m(6) * (t67 * t13 + t70 * t14) + t18 * mrSges(6,3) - t59 * mrSges(6,2) + t31 * t21 - t62 * t29;
t34 = -t45 * mrSges(5,1) + t46 * mrSges(5,2);
t7 = m(5) * t85 + qJDD(4) * mrSges(5,1) - t37 * mrSges(5,3) + qJD(4) * t41 + t70 * t11 + t67 * t12 - t46 * t34;
t78 = qJDD(2) * mrSges(4,3) + t73 * t81;
t8 = m(5) * t95 - qJDD(4) * mrSges(5,2) + t36 * mrSges(5,3) - qJD(4) * t42 - t67 * t11 + t70 * t12 + t45 * t34;
t5 = m(4) * t94 + t68 * t8 + t71 * t7 + (-m(4) * t40 - t78) * t64;
t6 = m(4) * t86 + t78 * t66 - t68 * t7 + t71 * t8;
t3 = m(3) * t93 - t73 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t64 * t5 + t66 * t6;
t87 = m(2) * t63 + t72 * t10 + t69 * t3;
t83 = -t66 * t5 - t64 * t6;
t4 = (m(2) + m(3)) * t52 + t83;
t1 = m(2) * t53 - t69 * t10 + t72 * t3;
t2 = [-m(1) * g(1) + t91 * t1 - t65 * t4, t1, t3, t6, t8, t12; -m(1) * g(2) + t65 * t1 + t91 * t4, t4, t10, t5, t7, t11; -m(1) * g(3) + t87, t87, -m(3) * t52 - t83, t81 * qJDD(2) - t73 * t97 + t74, -t76, -t77;];
f_new = t2;
