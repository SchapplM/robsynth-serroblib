% Calculate vector of cutting forces with Newton-Euler
% S5RRPRP2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPRP2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP2_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP2_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP2_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP2_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP2_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:29
% EndTime: 2019-12-31 19:49:30
% DurationCPUTime: 0.51s
% Computational Cost: add. (5394->109), mult. (6924->140), div. (0->0), fcn. (3373->8), ass. (0->56)
t46 = qJD(1) + qJD(2);
t44 = t46 ^ 2;
t45 = qJDD(1) + qJDD(2);
t53 = sin(qJ(1));
t56 = cos(qJ(1));
t64 = t53 * g(1) - t56 * g(2);
t34 = qJDD(1) * pkin(1) + t64;
t58 = qJD(1) ^ 2;
t60 = -t56 * g(1) - t53 * g(2);
t35 = -t58 * pkin(1) + t60;
t52 = sin(qJ(2));
t55 = cos(qJ(2));
t62 = t55 * t34 - t52 * t35;
t21 = t45 * pkin(2) + t62;
t69 = t52 * t34 + t55 * t35;
t22 = -t44 * pkin(2) + t69;
t49 = sin(pkin(8));
t50 = cos(pkin(8));
t63 = t50 * t21 - t49 * t22;
t16 = -t45 * pkin(3) - t44 * pkin(7) - t63;
t51 = sin(qJ(4));
t54 = cos(qJ(4));
t66 = qJD(4) * t46;
t28 = t51 * t45 + t54 * t66;
t29 = t54 * t45 - t51 * t66;
t77 = t46 * t51;
t36 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t77;
t37 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t77;
t76 = t46 * t54;
t39 = mrSges(6,2) * t76 + qJD(4) * mrSges(6,3);
t67 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t76 + t39;
t72 = m(6) * (-t29 * pkin(4) - t28 * qJ(5) + (-0.2e1 * qJD(5) * t51 + (pkin(4) * t51 - qJ(5) * t54) * qJD(4)) * t46 + t16) - t29 * mrSges(6,1);
t82 = (-t67 * t54 + (t36 - t37) * t51) * t46 + m(5) * t16 - t29 * mrSges(5,1) + (mrSges(5,2) - mrSges(6,3)) * t28 + t72;
t79 = -m(2) - m(3);
t70 = t49 * t21 + t50 * t22;
t17 = -t44 * pkin(3) + t45 * pkin(7) + t70;
t25 = (-pkin(4) * t54 - qJ(5) * t51) * t46;
t57 = qJD(4) ^ 2;
t48 = -g(3) + qJDD(3);
t75 = t54 * t48;
t78 = m(6) * (-qJDD(4) * pkin(4) - t57 * qJ(5) - t75 + qJDD(5) + (t25 * t46 + t17) * t51);
t73 = mrSges(5,3) + mrSges(6,2);
t71 = t54 * t17 + t51 * t48;
t26 = (-mrSges(6,1) * t54 - mrSges(6,3) * t51) * t46;
t27 = (-mrSges(5,1) * t54 + mrSges(5,2) * t51) * t46;
t10 = m(5) * (-t51 * t17 + t75) - t78 + (-t26 - t27) * t77 - t73 * t28 + (mrSges(5,1) + mrSges(6,1)) * qJDD(4) + t67 * qJD(4);
t61 = m(6) * (-t57 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t25 * t76 + t71) + t26 * t76 + qJD(4) * t37 + qJDD(4) * mrSges(6,3);
t9 = m(5) * t71 - qJDD(4) * mrSges(5,2) - qJD(4) * t36 + t27 * t76 + t73 * t29 + t61;
t65 = m(4) * t48 + t54 * t10 + t51 * t9;
t6 = m(4) * t63 + t45 * mrSges(4,1) - t44 * mrSges(4,2) - t82;
t5 = m(4) * t70 - t44 * mrSges(4,1) - t45 * mrSges(4,2) - t51 * t10 + t54 * t9;
t4 = m(3) * t69 - t44 * mrSges(3,1) - t45 * mrSges(3,2) - t49 * t6 + t50 * t5;
t3 = m(3) * t62 + t45 * mrSges(3,1) - t44 * mrSges(3,2) + t49 * t5 + t50 * t6;
t2 = m(2) * t60 - t58 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t52 * t3 + t55 * t4;
t1 = m(2) * t64 + qJDD(1) * mrSges(2,1) - t58 * mrSges(2,2) + t55 * t3 + t52 * t4;
t7 = [-m(1) * g(1) - t53 * t1 + t56 * t2, t2, t4, t5, t9, t29 * mrSges(6,2) + t61; -m(1) * g(2) + t56 * t1 + t53 * t2, t1, t3, t6, t10, -t28 * mrSges(6,3) + (-t51 * t37 - t54 * t39) * t46 + t72; (-m(1) + t79) * g(3) + t65, t79 * g(3) + t65, -m(3) * g(3) + t65, t65, t82, -qJDD(4) * mrSges(6,1) + t28 * mrSges(6,2) - qJD(4) * t39 + t26 * t77 + t78;];
f_new = t7;
