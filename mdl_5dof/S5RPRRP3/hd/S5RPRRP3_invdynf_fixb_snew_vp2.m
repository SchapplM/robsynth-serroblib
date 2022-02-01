% Calculate vector of cutting forces with Newton-Euler
% S5RPRRP3
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRRP3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP3_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP3_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP3_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:29:37
% EndTime: 2022-01-23 09:29:39
% DurationCPUTime: 0.80s
% Computational Cost: add. (6692->134), mult. (13296->172), div. (0->0), fcn. (7843->8), ass. (0->65)
t68 = sin(qJ(3));
t71 = cos(qJ(3));
t89 = qJD(1) * qJD(3);
t83 = t71 * t89;
t51 = t68 * qJDD(1) + t83;
t52 = t71 * qJDD(1) - t68 * t89;
t91 = qJD(1) * t68;
t53 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t91;
t90 = qJD(1) * t71;
t54 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t90;
t73 = qJD(1) ^ 2;
t67 = sin(qJ(4));
t70 = cos(qJ(4));
t46 = (t67 * t71 + t68 * t70) * qJD(1);
t26 = -t46 * qJD(4) - t67 * t51 + t70 * t52;
t45 = (-t67 * t68 + t70 * t71) * qJD(1);
t27 = t45 * qJD(4) + t70 * t51 + t67 * t52;
t62 = qJD(3) + qJD(4);
t36 = -t62 * mrSges(6,2) + t45 * mrSges(6,3);
t37 = -t62 * mrSges(5,2) + t45 * mrSges(5,3);
t40 = t62 * mrSges(5,1) - t46 * mrSges(5,3);
t55 = qJD(3) * pkin(3) - pkin(7) * t91;
t63 = t71 ^ 2;
t69 = sin(qJ(1));
t72 = cos(qJ(1));
t84 = t69 * g(1) - t72 * g(2);
t48 = qJDD(1) * pkin(1) + t84;
t79 = -t72 * g(1) - t69 * g(2);
t50 = -t73 * pkin(1) + t79;
t65 = sin(pkin(8));
t66 = cos(pkin(8));
t80 = t66 * t48 - t65 * t50;
t76 = -qJDD(1) * pkin(2) - t80;
t74 = -t52 * pkin(3) + t55 * t91 + (-pkin(7) * t63 - pkin(6)) * t73 + t76;
t38 = t62 * pkin(4) - t46 * qJ(5);
t39 = t62 * mrSges(6,1) - t46 * mrSges(6,3);
t41 = t45 ^ 2;
t85 = m(6) * (-t26 * pkin(4) - t41 * qJ(5) + t46 * t38 + qJDD(5) + t74) + t27 * mrSges(6,2) + t46 * t39;
t75 = m(5) * t74 + t27 * mrSges(5,2) - (mrSges(5,1) + mrSges(6,1)) * t26 + t46 * t40 - (t37 + t36) * t45 + t85;
t101 = (t68 * t53 - t71 * t54) * qJD(1) - t52 * mrSges(4,1) + t51 * mrSges(4,2) + m(4) * (-t73 * pkin(6) + t76) + t75;
t92 = t65 * t48 + t66 * t50;
t32 = -t73 * pkin(2) + qJDD(1) * pkin(6) + t92;
t64 = -g(3) + qJDD(2);
t81 = -t68 * t32 + t71 * t64;
t18 = (-t51 + t83) * pkin(7) + (t68 * t71 * t73 + qJDD(3)) * pkin(3) + t81;
t94 = t71 * t32 + t68 * t64;
t19 = -t63 * t73 * pkin(3) + t52 * pkin(7) - qJD(3) * t55 + t94;
t95 = t67 * t18 + t70 * t19;
t34 = -t45 * mrSges(5,1) + t46 * mrSges(5,2);
t61 = qJDD(3) + qJDD(4);
t33 = -t45 * mrSges(6,1) + t46 * mrSges(6,2);
t86 = m(6) * (-t41 * pkin(4) + t26 * qJ(5) + 0.2e1 * qJD(5) * t45 - t62 * t38 + t95) + t45 * t33 + t26 * mrSges(6,3);
t10 = m(5) * t95 + t26 * mrSges(5,3) + t45 * t34 + (-t40 - t39) * t62 + (-mrSges(5,2) - mrSges(6,2)) * t61 + t86;
t49 = (-mrSges(4,1) * t71 + mrSges(4,2) * t68) * qJD(1);
t82 = t70 * t18 - t67 * t19;
t87 = m(6) * (-0.2e1 * qJD(5) * t46 + (t45 * t62 - t27) * qJ(5) + (t45 * t46 + t61) * pkin(4) + t82) + t62 * t36 + t61 * mrSges(6,1);
t9 = m(5) * t82 + t61 * mrSges(5,1) + t62 * t37 + (-t34 - t33) * t46 + (-mrSges(5,3) - mrSges(6,3)) * t27 + t87;
t6 = m(4) * t81 + qJDD(3) * mrSges(4,1) - t51 * mrSges(4,3) + qJD(3) * t54 + t67 * t10 - t49 * t91 + t70 * t9;
t7 = m(4) * t94 - qJDD(3) * mrSges(4,2) + t52 * mrSges(4,3) - qJD(3) * t53 + t70 * t10 + t49 * t90 - t67 * t9;
t88 = m(3) * t64 + t71 * t6 + t68 * t7;
t8 = m(3) * t80 + qJDD(1) * mrSges(3,1) - t73 * mrSges(3,2) - t101;
t3 = m(3) * t92 - t73 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t68 * t6 + t71 * t7;
t2 = m(2) * t79 - t73 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t66 * t3 - t65 * t8;
t1 = m(2) * t84 + qJDD(1) * mrSges(2,1) - t73 * mrSges(2,2) + t65 * t3 + t66 * t8;
t4 = [-m(1) * g(1) - t69 * t1 + t72 * t2, t2, t3, t7, t10, -t61 * mrSges(6,2) - t62 * t39 + t86; -m(1) * g(2) + t72 * t1 + t69 * t2, t1, t8, t6, t9, -t27 * mrSges(6,3) - t46 * t33 + t87; (-m(1) - m(2)) * g(3) + t88, -m(2) * g(3) + t88, t88, t101, t75, -t26 * mrSges(6,1) - t45 * t36 + t85;];
f_new = t4;
