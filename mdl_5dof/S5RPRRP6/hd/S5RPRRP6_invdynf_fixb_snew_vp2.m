% Calculate vector of cutting forces with Newton-Euler
% S5RPRRP6
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
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRRP6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP6_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP6_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP6_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP6_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP6_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:42:12
% EndTime: 2019-12-31 18:42:14
% DurationCPUTime: 0.68s
% Computational Cost: add. (6052->133), mult. (11466->167), div. (0->0), fcn. (6562->8), ass. (0->66)
t70 = qJD(1) ^ 2;
t65 = sin(qJ(1));
t68 = cos(qJ(1));
t80 = t65 * g(1) - t68 * g(2);
t48 = qJDD(1) * pkin(1) + t80;
t74 = -t68 * g(1) - t65 * g(2);
t50 = -t70 * pkin(1) + t74;
t61 = sin(pkin(8));
t62 = cos(pkin(8));
t75 = t62 * t48 - t61 * t50;
t24 = -qJDD(1) * pkin(2) - t70 * pkin(6) - t75;
t64 = sin(qJ(3));
t67 = cos(qJ(3));
t85 = qJD(1) * qJD(3);
t78 = t67 * t85;
t52 = t64 * qJDD(1) + t78;
t79 = t64 * t85;
t53 = t67 * qJDD(1) - t79;
t87 = qJD(1) * t64;
t54 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t87;
t86 = t67 * qJD(1);
t55 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t86;
t63 = sin(qJ(4));
t66 = cos(qJ(4));
t46 = t66 * qJD(3) - t63 * t87;
t30 = t46 * qJD(4) + t63 * qJDD(3) + t66 * t52;
t47 = t63 * qJD(3) + t66 * t87;
t32 = -t46 * mrSges(6,1) + t47 * mrSges(6,2);
t33 = -t46 * mrSges(5,1) + t47 * mrSges(5,2);
t56 = qJD(4) - t86;
t35 = -t56 * mrSges(5,2) + t46 * mrSges(5,3);
t45 = qJDD(4) - t53;
t17 = (-t52 - t78) * pkin(7) + (-t53 + t79) * pkin(3) + t24;
t51 = (-pkin(3) * t67 - pkin(7) * t64) * qJD(1);
t69 = qJD(3) ^ 2;
t88 = t61 * t48 + t62 * t50;
t25 = -t70 * pkin(2) + qJDD(1) * pkin(6) + t88;
t60 = -g(3) + qJDD(2);
t90 = t67 * t25 + t64 * t60;
t20 = -t69 * pkin(3) + qJDD(3) * pkin(7) + t51 * t86 + t90;
t77 = t66 * t17 - t63 * t20;
t34 = -t56 * mrSges(6,2) + t46 * mrSges(6,3);
t83 = m(6) * (-0.2e1 * qJD(5) * t47 + (t46 * t56 - t30) * qJ(5) + (t46 * t47 + t45) * pkin(4) + t77) + t56 * t34 + t45 * mrSges(6,1);
t7 = m(5) * t77 + t45 * mrSges(5,1) + t56 * t35 + (-t33 - t32) * t47 + (-mrSges(5,3) - mrSges(6,3)) * t30 + t83;
t29 = -t47 * qJD(4) + t66 * qJDD(3) - t63 * t52;
t37 = t56 * mrSges(6,1) - t47 * mrSges(6,3);
t38 = t56 * mrSges(5,1) - t47 * mrSges(5,3);
t36 = t56 * pkin(4) - t47 * qJ(5);
t44 = t46 ^ 2;
t91 = t63 * t17 + t66 * t20;
t82 = m(6) * (-t44 * pkin(4) + t29 * qJ(5) + 0.2e1 * qJD(5) * t46 - t56 * t36 + t91) + t46 * t32 + t29 * mrSges(6,3);
t8 = m(5) * t91 + t29 * mrSges(5,3) + t46 * t33 + (-t38 - t37) * t56 + (-mrSges(5,2) - mrSges(6,2)) * t45 + t82;
t94 = m(4) * t24 - t53 * mrSges(4,1) + t52 * mrSges(4,2) + t63 * t8 + t66 * t7 + (t64 * t54 - t67 * t55) * qJD(1);
t76 = -t64 * t25 + t67 * t60;
t19 = -qJDD(3) * pkin(3) - t69 * pkin(7) + t51 * t87 - t76;
t81 = m(6) * (-t29 * pkin(4) - t44 * qJ(5) + t47 * t36 + qJDD(5) + t19) + t30 * mrSges(6,2) + t47 * t37;
t93 = m(5) * t19 + t30 * mrSges(5,2) - (mrSges(5,1) + mrSges(6,1)) * t29 + t47 * t38 - (t35 + t34) * t46 + t81;
t49 = (-mrSges(4,1) * t67 + mrSges(4,2) * t64) * qJD(1);
t10 = m(4) * t76 + qJDD(3) * mrSges(4,1) - t52 * mrSges(4,3) + qJD(3) * t55 - t49 * t87 - t93;
t6 = m(4) * t90 - qJDD(3) * mrSges(4,2) + t53 * mrSges(4,3) - qJD(3) * t54 + t49 * t86 - t63 * t7 + t66 * t8;
t84 = m(3) * t60 + t67 * t10 + t64 * t6;
t4 = m(3) * t75 + qJDD(1) * mrSges(3,1) - t70 * mrSges(3,2) - t94;
t3 = m(3) * t88 - t70 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t64 * t10 + t67 * t6;
t2 = m(2) * t74 - t70 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t62 * t3 - t61 * t4;
t1 = m(2) * t80 + qJDD(1) * mrSges(2,1) - t70 * mrSges(2,2) + t61 * t3 + t62 * t4;
t5 = [-m(1) * g(1) - t65 * t1 + t68 * t2, t2, t3, t6, t8, -t45 * mrSges(6,2) - t56 * t37 + t82; -m(1) * g(2) + t68 * t1 + t65 * t2, t1, t4, t10, t7, -t30 * mrSges(6,3) - t47 * t32 + t83; (-m(1) - m(2)) * g(3) + t84, -m(2) * g(3) + t84, t84, t94, t93, -t29 * mrSges(6,1) - t46 * t34 + t81;];
f_new = t5;
