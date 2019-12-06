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
% Datum: 2019-12-05 18:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:22:03
% EndTime: 2019-12-05 18:22:04
% DurationCPUTime: 0.50s
% Computational Cost: add. (5469->110), mult. (7066->137), div. (0->0), fcn. (3444->8), ass. (0->58)
t49 = qJDD(1) + qJDD(2);
t56 = sin(qJ(4));
t59 = cos(qJ(4));
t50 = qJD(1) + qJD(2);
t74 = qJD(4) * t50;
t68 = t59 * t74;
t29 = t56 * t49 + t68;
t30 = t59 * t49 - t56 * t74;
t83 = t50 * t56;
t39 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t83;
t82 = t50 * t59;
t40 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t82;
t41 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t82;
t48 = t50 ^ 2;
t58 = sin(qJ(1));
t61 = cos(qJ(1));
t75 = t61 * g(2) + t58 * g(3);
t35 = qJDD(1) * pkin(1) + t75;
t62 = qJD(1) ^ 2;
t67 = t58 * g(2) - t61 * g(3);
t36 = -t62 * pkin(1) + t67;
t57 = sin(qJ(2));
t60 = cos(qJ(2));
t65 = t60 * t35 - t57 * t36;
t21 = t49 * pkin(2) + t65;
t77 = t57 * t35 + t60 * t36;
t22 = -t48 * pkin(2) + t77;
t54 = sin(pkin(8));
t55 = cos(pkin(8));
t66 = t55 * t21 - t54 * t22;
t64 = -t49 * pkin(3) - t66;
t37 = qJD(4) * pkin(4) - qJ(5) * t83;
t38 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t83;
t52 = t59 ^ 2;
t69 = m(6) * (t37 * t83 - t30 * pkin(4) + qJDD(5) + (-qJ(5) * t52 - pkin(7)) * t48 + t64) + t38 * t83 + t29 * mrSges(6,2);
t88 = -(-t56 * t39 + (t40 + t41) * t59) * t50 - (mrSges(5,1) + mrSges(6,1)) * t30 + m(5) * (-t48 * pkin(7) + t64) + t29 * mrSges(5,2) + t69;
t85 = -m(2) - m(3);
t84 = pkin(4) * t48;
t78 = t54 * t21 + t55 * t22;
t17 = -t48 * pkin(3) + t49 * pkin(7) + t78;
t53 = -g(1) + qJDD(3);
t79 = t59 * t17 + t56 * t53;
t73 = qJD(5) * t50;
t28 = (-mrSges(5,1) * t59 + mrSges(5,2) * t56) * t50;
t27 = (-mrSges(6,1) * t59 + mrSges(6,2) * t56) * t50;
t70 = m(6) * (t30 * qJ(5) - qJD(4) * t37 - t52 * t84 + 0.2e1 * t59 * t73 + t79) + t27 * t82 + t30 * mrSges(6,3);
t10 = m(5) * t79 + t30 * mrSges(5,3) + t28 * t82 + (-mrSges(5,2) - mrSges(6,2)) * qJDD(4) + (-t39 - t38) * qJD(4) + t70;
t43 = t59 * t53;
t71 = m(6) * (qJDD(4) * pkin(4) + t43 + (-t29 + t68) * qJ(5) + (t59 * t84 - t17 - 0.2e1 * t73) * t56) + qJD(4) * t40 + qJDD(4) * mrSges(6,1);
t9 = m(5) * (-t56 * t17 + t43) + qJDD(4) * mrSges(5,1) + qJD(4) * t41 + (-t27 - t28) * t83 + (-mrSges(5,3) - mrSges(6,3)) * t29 + t71;
t72 = m(4) * t53 + t56 * t10 + t59 * t9;
t6 = m(4) * t66 + t49 * mrSges(4,1) - t48 * mrSges(4,2) - t88;
t5 = m(4) * t78 - t48 * mrSges(4,1) - t49 * mrSges(4,2) + t59 * t10 - t56 * t9;
t4 = m(3) * t77 - t48 * mrSges(3,1) - t49 * mrSges(3,2) + t55 * t5 - t54 * t6;
t3 = m(3) * t65 + t49 * mrSges(3,1) - t48 * mrSges(3,2) + t54 * t5 + t55 * t6;
t2 = m(2) * t75 + qJDD(1) * mrSges(2,1) - t62 * mrSges(2,2) + t60 * t3 + t57 * t4;
t1 = m(2) * t67 - t62 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t57 * t3 + t60 * t4;
t7 = [(-m(1) + t85) * g(1) + t72, t1, t4, t5, t10, -qJDD(4) * mrSges(6,2) - qJD(4) * t38 + t70; -m(1) * g(2) - t58 * t1 - t61 * t2, t2, t3, t6, t9, -t29 * mrSges(6,3) - t27 * t83 + t71; -m(1) * g(3) + t61 * t1 - t58 * t2, t85 * g(1) + t72, -m(3) * g(1) + t72, t72, t88, -t30 * mrSges(6,1) - t40 * t82 + t69;];
f_new = t7;
