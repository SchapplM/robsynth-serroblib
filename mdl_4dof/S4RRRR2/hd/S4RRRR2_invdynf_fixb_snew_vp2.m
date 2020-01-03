% Calculate vector of cutting forces with Newton-Euler
% S4RRRR2
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x5]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RRRR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR2_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR2_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR2_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR2_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR2_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:09
% EndTime: 2019-12-31 17:23:10
% DurationCPUTime: 0.44s
% Computational Cost: add. (4763->100), mult. (6130->136), div. (0->0), fcn. (3457->8), ass. (0->57)
t41 = qJDD(1) + qJDD(2);
t46 = sin(qJ(3));
t50 = cos(qJ(3));
t43 = qJD(1) + qJD(2);
t62 = qJD(3) * t43;
t28 = t46 * t41 + t50 * t62;
t29 = t50 * t41 - t46 * t62;
t66 = t43 * t46;
t35 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t66;
t65 = t43 * t50;
t36 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t65;
t39 = t43 ^ 2;
t45 = sin(qJ(4));
t49 = cos(qJ(4));
t26 = (t45 * t50 + t46 * t49) * t43;
t16 = -t26 * qJD(4) - t45 * t28 + t49 * t29;
t25 = (-t45 * t46 + t49 * t50) * t43;
t17 = t25 * qJD(4) + t49 * t28 + t45 * t29;
t42 = qJD(3) + qJD(4);
t23 = -t42 * mrSges(5,2) + t25 * mrSges(5,3);
t24 = t42 * mrSges(5,1) - t26 * mrSges(5,3);
t37 = qJD(3) * pkin(3) - pkin(7) * t66;
t44 = t50 ^ 2;
t48 = sin(qJ(1));
t52 = cos(qJ(1));
t61 = t48 * g(1) - t52 * g(2);
t33 = qJDD(1) * pkin(1) + t61;
t53 = qJD(1) ^ 2;
t58 = -t52 * g(1) - t48 * g(2);
t34 = -t53 * pkin(1) + t58;
t47 = sin(qJ(2));
t51 = cos(qJ(2));
t59 = t51 * t33 - t47 * t34;
t56 = -t41 * pkin(2) - t59;
t55 = t16 * mrSges(5,1) + t25 * t23 - m(5) * (t37 * t66 - t29 * pkin(3) + (-pkin(7) * t44 - pkin(6)) * t39 + t56) - t17 * mrSges(5,2) - t26 * t24;
t70 = (t46 * t35 - t50 * t36) * t43 + m(4) * (-t39 * pkin(6) + t56) - t29 * mrSges(4,1) + t28 * mrSges(4,2) - t55;
t69 = -m(2) - m(3);
t63 = t47 * t33 + t51 * t34;
t22 = -t39 * pkin(2) + t41 * pkin(6) + t63;
t64 = t46 * t22;
t67 = pkin(3) * t39;
t11 = qJDD(3) * pkin(3) - t28 * pkin(7) - t64 + (pkin(7) * t62 + t46 * t67 - g(3)) * t50;
t60 = -t46 * g(3) + t50 * t22;
t12 = t29 * pkin(7) - qJD(3) * t37 - t44 * t67 + t60;
t20 = -t25 * mrSges(5,1) + t26 * mrSges(5,2);
t40 = qJDD(3) + qJDD(4);
t10 = m(5) * (t45 * t11 + t49 * t12) + t16 * mrSges(5,3) - t40 * mrSges(5,2) + t25 * t20 - t42 * t24;
t27 = (-mrSges(4,1) * t50 + mrSges(4,2) * t46) * t43;
t9 = m(5) * (t49 * t11 - t45 * t12) - t17 * mrSges(5,3) + t40 * mrSges(5,1) - t26 * t20 + t42 * t23;
t6 = m(4) * (-t50 * g(3) - t64) - t28 * mrSges(4,3) + qJDD(3) * mrSges(4,1) - t27 * t66 + qJD(3) * t36 + t45 * t10 + t49 * t9;
t7 = m(4) * t60 - qJDD(3) * mrSges(4,2) + t29 * mrSges(4,3) - qJD(3) * t35 + t49 * t10 + t27 * t65 - t45 * t9;
t68 = t46 * t7 + t50 * t6;
t8 = m(3) * t59 + t41 * mrSges(3,1) - t39 * mrSges(3,2) - t70;
t3 = m(3) * t63 - t39 * mrSges(3,1) - t41 * mrSges(3,2) - t46 * t6 + t50 * t7;
t2 = m(2) * t58 - t53 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t51 * t3 - t47 * t8;
t1 = m(2) * t61 + qJDD(1) * mrSges(2,1) - t53 * mrSges(2,2) + t47 * t3 + t51 * t8;
t4 = [-m(1) * g(1) - t48 * t1 + t52 * t2, t2, t3, t7, t10; -m(1) * g(2) + t52 * t1 + t48 * t2, t1, t8, t6, t9; (-m(1) + t69) * g(3) + t68, t69 * g(3) + t68, -m(3) * g(3) + t68, t70, -t55;];
f_new = t4;
