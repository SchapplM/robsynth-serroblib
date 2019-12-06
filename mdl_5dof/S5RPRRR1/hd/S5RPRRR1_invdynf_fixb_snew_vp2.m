% Calculate vector of cutting forces with Newton-Euler
% S5RPRRR1
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% Datum: 2019-12-05 18:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRRR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(1,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR1_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_invdynf_fixb_snew_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR1_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR1_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR1_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:08:39
% EndTime: 2019-12-05 18:08:40
% DurationCPUTime: 0.49s
% Computational Cost: add. (4586->113), mult. (8239->145), div. (0->0), fcn. (5828->8), ass. (0->60)
t69 = -m(2) - m(3);
t51 = sin(qJ(1));
t55 = cos(qJ(1));
t63 = -t55 * g(1) - t51 * g(2);
t32 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t63;
t50 = sin(qJ(3));
t54 = cos(qJ(3));
t28 = -t50 * g(3) + t54 * t32;
t56 = qJD(1) ^ 2;
t62 = t51 * g(1) - t55 * g(2);
t36 = -t56 * qJ(2) + qJDD(2) - t62;
t49 = sin(qJ(4));
t53 = cos(qJ(4));
t17 = t49 * t28 - t53 * t36;
t66 = qJD(1) * t50;
t39 = t53 * qJD(3) - t49 * t66;
t64 = qJD(1) * qJD(3);
t42 = t50 * qJDD(1) + t54 * t64;
t23 = t39 * qJD(4) + t49 * qJDD(3) + t53 * t42;
t40 = t49 * qJD(3) + t53 * t66;
t24 = -t39 * mrSges(5,1) + t40 * mrSges(5,2);
t65 = t54 * qJD(1);
t46 = qJD(4) - t65;
t29 = -t46 * mrSges(5,2) + t39 * mrSges(5,3);
t43 = t54 * qJDD(1) - t50 * t64;
t38 = qJDD(4) - t43;
t48 = sin(qJ(5));
t52 = cos(qJ(5));
t26 = t52 * t40 + t48 * t46;
t14 = -t26 * qJD(5) - t48 * t23 + t52 * t38;
t25 = -t48 * t40 + t52 * t46;
t15 = t25 * qJD(5) + t52 * t23 + t48 * t38;
t35 = qJD(5) - t39;
t19 = -t35 * mrSges(6,2) + t25 * mrSges(6,3);
t20 = t35 * mrSges(6,1) - t26 * mrSges(6,3);
t60 = t14 * mrSges(6,1) - t15 * mrSges(6,2) + t25 * t19 - t26 * t20;
t11 = t38 * mrSges(5,1) - t23 * mrSges(5,3) - t40 * t24 + t46 * t29 + (-m(5) - m(6)) * t17 + t60;
t41 = (-mrSges(4,1) * t54 + mrSges(4,2) * t50) * qJD(1);
t44 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t66;
t16 = -t25 * mrSges(6,1) + t26 * mrSges(6,2);
t18 = t53 * t28 + t49 * t36;
t22 = -t40 * qJD(4) + t53 * qJDD(3) - t49 * t42;
t21 = qJDD(5) - t22;
t27 = t54 * g(3) + t50 * t32;
t12 = m(6) * (-t48 * t18 + t52 * t27) - t15 * mrSges(6,3) + t21 * mrSges(6,1) - t26 * t16 + t35 * t19;
t13 = m(6) * (t52 * t18 + t48 * t27) + t14 * mrSges(6,3) - t21 * mrSges(6,2) + t25 * t16 - t35 * t20;
t30 = t46 * mrSges(5,1) - t40 * mrSges(5,3);
t9 = m(5) * t18 - t38 * mrSges(5,2) + t22 * mrSges(5,3) - t48 * t12 + t52 * t13 + t39 * t24 - t46 * t30;
t5 = m(4) * t28 - qJDD(3) * mrSges(4,2) + t43 * mrSges(4,3) - qJD(3) * t44 - t49 * t11 + t41 * t65 + t53 * t9;
t45 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t65;
t57 = t22 * mrSges(5,1) - t23 * mrSges(5,2) - t52 * t12 - t48 * t13 + t39 * t29 - t40 * t30;
t7 = -t41 * t66 + qJDD(3) * mrSges(4,1) - t42 * mrSges(4,3) + qJD(3) * t45 + (-m(4) - m(5)) * t27 + t57;
t68 = t50 * t5 + t54 * t7;
t67 = mrSges(2,1) + mrSges(3,1);
t61 = m(3) * t32 + qJDD(1) * mrSges(3,3) + t54 * t5 - t50 * t7;
t59 = -m(4) * t36 + t43 * mrSges(4,1) - t42 * mrSges(4,2) - t53 * t11 - t44 * t66 + t45 * t65 - t49 * t9;
t58 = -m(3) * t36 + t59;
t2 = m(2) * t62 + (-mrSges(2,2) + mrSges(3,3)) * t56 + t67 * qJDD(1) + t58;
t1 = m(2) * t63 - qJDD(1) * mrSges(2,2) - t67 * t56 + t61;
t3 = [-m(1) * g(1) + t55 * t1 - t51 * t2, t1, -(t56 * mrSges(3,1)) + t61, t5, t9, t13; -m(1) * g(2) + t51 * t1 + t55 * t2, t2, -m(3) * g(3) + t68, t7, t11, t12; (-m(1) + t69) * g(3) + t68, t69 * g(3) + t68, -qJDD(1) * mrSges(3,1) - t56 * mrSges(3,3) - t58, -t59, m(5) * t27 - t57, m(6) * t17 - t60;];
f_new = t3;
