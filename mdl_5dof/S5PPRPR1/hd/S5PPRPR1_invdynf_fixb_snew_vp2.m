% Calculate vector of cutting forces with Newton-Euler
% S5PPRPR1
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PPRPR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR1_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR1_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR1_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR1_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR1_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:00:56
% EndTime: 2019-12-05 15:00:57
% DurationCPUTime: 0.54s
% Computational Cost: add. (5152->86), mult. (9944->120), div. (0->0), fcn. (6839->10), ass. (0->58)
t58 = qJD(3) ^ 2;
t51 = cos(pkin(9));
t46 = t51 ^ 2;
t48 = sin(pkin(9));
t74 = t48 ^ 2 + t46;
t78 = t74 * mrSges(5,3);
t77 = pkin(4) * t58;
t50 = sin(pkin(7));
t53 = cos(pkin(7));
t41 = -t53 * g(1) - t50 * g(2);
t47 = -g(3) + qJDD(1);
t49 = sin(pkin(8));
t52 = cos(pkin(8));
t32 = -t49 * t41 + t52 * t47;
t33 = t52 * t41 + t49 * t47;
t55 = sin(qJ(3));
t57 = cos(qJ(3));
t76 = t55 * t32 + t57 * t33;
t67 = t50 * g(1) - t53 * g(2);
t40 = qJDD(2) - t67;
t72 = qJD(3) * qJD(4);
t75 = t51 * t40 - 0.2e1 * t48 * t72;
t73 = pkin(6) * qJDD(3);
t21 = -t58 * pkin(3) + qJDD(3) * qJ(4) + t76;
t15 = (t51 * t77 - t21 - t73) * t48 + t75;
t69 = t48 * t40 + (t21 + 0.2e1 * t72) * t51;
t16 = -t46 * t77 + t51 * t73 + t69;
t54 = sin(qJ(5));
t56 = cos(qJ(5));
t62 = -t48 * t54 + t51 * t56;
t34 = t62 * qJD(3);
t63 = t48 * t56 + t51 * t54;
t35 = t63 * qJD(3);
t23 = -t34 * mrSges(6,1) + t35 * mrSges(6,2);
t26 = t34 * qJD(5) + t63 * qJDD(3);
t30 = -qJD(5) * mrSges(6,2) + t34 * mrSges(6,3);
t13 = m(6) * (t56 * t15 - t54 * t16) - t26 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t35 * t23 + qJD(5) * t30;
t25 = -t35 * qJD(5) + t62 * qJDD(3);
t31 = qJD(5) * mrSges(6,1) - t35 * mrSges(6,3);
t14 = m(6) * (t54 * t15 + t56 * t16) + t25 * mrSges(6,3) - qJDD(5) * mrSges(6,2) + t34 * t23 - qJD(5) * t31;
t64 = -t51 * mrSges(5,1) + t48 * mrSges(5,2);
t61 = qJDD(3) * mrSges(5,3) + t58 * t64;
t10 = m(5) * t75 + t54 * t14 + t56 * t13 + (-m(5) * t21 - t61) * t48;
t11 = m(5) * t69 - t54 * t13 + t56 * t14 + t61 * t51;
t71 = m(4) * t40 + t51 * t10 + t48 * t11;
t68 = t57 * t32 - t55 * t33;
t66 = qJDD(4) - t68;
t60 = t25 * mrSges(6,1) + t34 * t30 - m(6) * ((-pkin(4) * t51 - pkin(3)) * qJDD(3) + (-t74 * pkin(6) - qJ(4)) * t58 + t66) - t35 * t31 - t26 * mrSges(6,2);
t59 = m(5) * (-qJDD(3) * pkin(3) - t58 * qJ(4) + t66) - t60;
t12 = m(4) * t68 + (-mrSges(4,2) + t78) * t58 + (mrSges(4,1) - t64) * qJDD(3) - t59;
t6 = m(4) * t76 - t58 * mrSges(4,1) - qJDD(3) * mrSges(4,2) - t48 * t10 + t51 * t11;
t4 = m(3) * t32 + t57 * t12 + t55 * t6;
t5 = m(3) * t33 - t55 * t12 + t57 * t6;
t70 = m(2) * t47 + t52 * t4 + t49 * t5;
t65 = m(3) * t40 + t71;
t7 = m(2) * t67 - t65;
t1 = m(2) * t41 - t49 * t4 + t52 * t5;
t2 = [-m(1) * g(1) + t53 * t1 - t50 * t7, t1, t5, t6, t11, t14; -m(1) * g(2) + t50 * t1 + t53 * t7, t7, t4, t12, t10, t13; -m(1) * g(3) + t70, t70, t65, t71, qJDD(3) * t64 - t58 * t78 + t59, -t60;];
f_new = t2;
