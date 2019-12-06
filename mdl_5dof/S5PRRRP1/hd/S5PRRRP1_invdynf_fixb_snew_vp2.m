% Calculate vector of cutting forces with Newton-Euler
% S5PRRRP1
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRRRP1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP1_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP1_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP1_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP1_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:39:52
% EndTime: 2019-12-05 16:39:53
% DurationCPUTime: 0.49s
% Computational Cost: add. (4536->102), mult. (6136->131), div. (0->0), fcn. (3444->8), ass. (0->57)
t48 = qJDD(2) + qJDD(3);
t55 = sin(qJ(4));
t58 = cos(qJ(4));
t49 = qJD(2) + qJD(3);
t74 = qJD(4) * t49;
t68 = t58 * t74;
t29 = t55 * t48 + t68;
t30 = t58 * t48 - t55 * t74;
t82 = t49 * t55;
t37 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t82;
t81 = t49 * t58;
t38 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t81;
t39 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t81;
t47 = t49 ^ 2;
t53 = sin(pkin(8));
t54 = cos(pkin(8));
t40 = t53 * g(1) - t54 * g(2);
t41 = -t54 * g(1) - t53 * g(2);
t57 = sin(qJ(2));
t60 = cos(qJ(2));
t65 = t60 * t40 - t57 * t41;
t21 = qJDD(2) * pkin(2) + t65;
t61 = qJD(2) ^ 2;
t76 = t57 * t40 + t60 * t41;
t22 = -t61 * pkin(2) + t76;
t56 = sin(qJ(3));
t59 = cos(qJ(3));
t66 = t59 * t21 - t56 * t22;
t63 = -t48 * pkin(3) - t66;
t35 = qJD(4) * pkin(4) - qJ(5) * t82;
t36 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t82;
t51 = t58 ^ 2;
t69 = m(6) * (t35 * t82 - t30 * pkin(4) + qJDD(5) + (-qJ(5) * t51 - pkin(7)) * t47 + t63) + t36 * t82 + t29 * mrSges(6,2);
t86 = -(-t55 * t37 + (t38 + t39) * t58) * t49 - (mrSges(5,1) + mrSges(6,1)) * t30 + m(5) * (-t47 * pkin(7) + t63) + t29 * mrSges(5,2) + t69;
t83 = pkin(4) * t47;
t77 = t56 * t21 + t59 * t22;
t17 = -t47 * pkin(3) + t48 * pkin(7) + t77;
t52 = -g(3) + qJDD(1);
t78 = t58 * t17 + t55 * t52;
t73 = qJD(5) * t49;
t28 = (-mrSges(5,1) * t58 + mrSges(5,2) * t55) * t49;
t27 = (-mrSges(6,1) * t58 + mrSges(6,2) * t55) * t49;
t70 = m(6) * (t30 * qJ(5) - qJD(4) * t35 - t51 * t83 + 0.2e1 * t58 * t73 + t78) + t27 * t81 + t30 * mrSges(6,3);
t10 = m(5) * t78 + t30 * mrSges(5,3) + t28 * t81 + (-mrSges(5,2) - mrSges(6,2)) * qJDD(4) + (-t37 - t36) * qJD(4) + t70;
t43 = t58 * t52;
t71 = m(6) * (qJDD(4) * pkin(4) + t43 + (-t29 + t68) * qJ(5) + (t58 * t83 - t17 - 0.2e1 * t73) * t55) + qJD(4) * t38 + qJDD(4) * mrSges(6,1);
t9 = m(5) * (-t55 * t17 + t43) + qJDD(4) * mrSges(5,1) + qJD(4) * t39 + (-t27 - t28) * t82 + (-mrSges(5,3) - mrSges(6,3)) * t29 + t71;
t72 = m(4) * t52 + t55 * t10 + t58 * t9;
t67 = m(3) * t52 + t72;
t64 = m(2) * t52 + t67;
t6 = m(4) * t66 + t48 * mrSges(4,1) - t47 * mrSges(4,2) - t86;
t5 = m(4) * t77 - t47 * mrSges(4,1) - t48 * mrSges(4,2) + t58 * t10 - t55 * t9;
t4 = m(3) * t76 - t61 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t59 * t5 - t56 * t6;
t3 = m(3) * t65 + qJDD(2) * mrSges(3,1) - t61 * mrSges(3,2) + t56 * t5 + t59 * t6;
t2 = m(2) * t41 - t57 * t3 + t60 * t4;
t1 = m(2) * t40 + t60 * t3 + t57 * t4;
t7 = [-m(1) * g(1) - t53 * t1 + t54 * t2, t2, t4, t5, t10, -qJDD(4) * mrSges(6,2) - qJD(4) * t36 + t70; -m(1) * g(2) + t54 * t1 + t53 * t2, t1, t3, t6, t9, -t29 * mrSges(6,3) - t27 * t82 + t71; -m(1) * g(3) + t64, t64, t67, t72, t86, -t30 * mrSges(6,1) - t38 * t81 + t69;];
f_new = t7;
