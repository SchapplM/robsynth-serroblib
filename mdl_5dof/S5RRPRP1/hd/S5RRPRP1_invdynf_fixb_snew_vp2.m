% Calculate vector of cutting forces with Newton-Euler
% S5RRPRP1
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
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPRP1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP1_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP1_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP1_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:19:37
% EndTime: 2022-01-20 10:19:38
% DurationCPUTime: 0.54s
% Computational Cost: add. (5469->110), mult. (7066->137), div. (0->0), fcn. (3444->8), ass. (0->58)
t47 = qJDD(1) + qJDD(2);
t54 = sin(qJ(4));
t57 = cos(qJ(4));
t48 = qJD(1) + qJD(2);
t73 = qJD(4) * t48;
t67 = t57 * t73;
t29 = t54 * t47 + t67;
t30 = t57 * t47 - t54 * t73;
t81 = t48 * t54;
t39 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t81;
t80 = t48 * t57;
t40 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t80;
t41 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t80;
t46 = t48 ^ 2;
t56 = sin(qJ(1));
t59 = cos(qJ(1));
t66 = t56 * g(1) - t59 * g(2);
t35 = qJDD(1) * pkin(1) + t66;
t60 = qJD(1) ^ 2;
t63 = -t59 * g(1) - t56 * g(2);
t36 = -t60 * pkin(1) + t63;
t55 = sin(qJ(2));
t58 = cos(qJ(2));
t64 = t58 * t35 - t55 * t36;
t21 = t47 * pkin(2) + t64;
t75 = t55 * t35 + t58 * t36;
t22 = -t46 * pkin(2) + t75;
t52 = sin(pkin(8));
t53 = cos(pkin(8));
t65 = t53 * t21 - t52 * t22;
t62 = -t47 * pkin(3) - t65;
t37 = qJD(4) * pkin(4) - qJ(5) * t81;
t38 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t81;
t50 = t57 ^ 2;
t68 = m(6) * (t37 * t81 - t30 * pkin(4) + qJDD(5) + (-qJ(5) * t50 - pkin(7)) * t46 + t62) + t38 * t81 + t29 * mrSges(6,2);
t86 = -(-t54 * t39 + (t40 + t41) * t57) * t48 - (mrSges(5,1) + mrSges(6,1)) * t30 + m(5) * (-t46 * pkin(7) + t62) + t29 * mrSges(5,2) + t68;
t83 = -m(2) - m(3);
t82 = pkin(4) * t46;
t76 = t52 * t21 + t53 * t22;
t17 = -t46 * pkin(3) + t47 * pkin(7) + t76;
t51 = -g(3) + qJDD(3);
t77 = t57 * t17 + t54 * t51;
t72 = qJD(5) * t48;
t28 = (-mrSges(5,1) * t57 + mrSges(5,2) * t54) * t48;
t27 = (-mrSges(6,1) * t57 + mrSges(6,2) * t54) * t48;
t69 = m(6) * (t30 * qJ(5) - qJD(4) * t37 - t50 * t82 + 0.2e1 * t57 * t72 + t77) + t27 * t80 + t30 * mrSges(6,3);
t10 = m(5) * t77 + t30 * mrSges(5,3) + t28 * t80 + (-mrSges(5,2) - mrSges(6,2)) * qJDD(4) + (-t39 - t38) * qJD(4) + t69;
t43 = t57 * t51;
t70 = m(6) * (qJDD(4) * pkin(4) + t43 + (-t29 + t67) * qJ(5) + (t57 * t82 - t17 - 0.2e1 * t72) * t54) + qJD(4) * t40 + qJDD(4) * mrSges(6,1);
t9 = m(5) * (-t54 * t17 + t43) + qJDD(4) * mrSges(5,1) + qJD(4) * t41 + (-t27 - t28) * t81 + (-mrSges(5,3) - mrSges(6,3)) * t29 + t70;
t71 = m(4) * t51 + t54 * t10 + t57 * t9;
t6 = m(4) * t65 + t47 * mrSges(4,1) - t46 * mrSges(4,2) - t86;
t5 = m(4) * t76 - t46 * mrSges(4,1) - t47 * mrSges(4,2) + t57 * t10 - t54 * t9;
t4 = m(3) * t75 - t46 * mrSges(3,1) - t47 * mrSges(3,2) + t53 * t5 - t52 * t6;
t3 = m(3) * t64 + t47 * mrSges(3,1) - t46 * mrSges(3,2) + t52 * t5 + t53 * t6;
t2 = m(2) * t63 - t60 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t55 * t3 + t58 * t4;
t1 = m(2) * t66 + qJDD(1) * mrSges(2,1) - t60 * mrSges(2,2) + t58 * t3 + t55 * t4;
t7 = [-m(1) * g(1) - t56 * t1 + t59 * t2, t2, t4, t5, t10, -qJDD(4) * mrSges(6,2) - qJD(4) * t38 + t69; -m(1) * g(2) + t59 * t1 + t56 * t2, t1, t3, t6, t9, -t29 * mrSges(6,3) - t27 * t81 + t70; (-m(1) + t83) * g(3) + t71, t83 * g(3) + t71, -m(3) * g(3) + t71, t71, t86, -t30 * mrSges(6,1) - t40 * t80 + t68;];
f_new = t7;
