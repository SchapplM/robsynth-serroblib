% Calculate vector of inverse dynamics joint torques for
% S4PRRR2
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(2,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR2_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR2_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_invdynJ_fixb_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR2_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR2_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR2_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:22
% EndTime: 2019-07-18 13:27:22
% DurationCPUTime: 0.30s
% Computational Cost: add. (318->93), mult. (594->123), div. (0->0), fcn. (272->10), ass. (0->52)
t37 = cos(qJ(3));
t58 = t37 * pkin(1);
t30 = qJDD(2) + qJDD(3);
t25 = qJDD(4) + t30;
t34 = sin(qJ(3));
t52 = qJD(3) * t37;
t15 = (qJD(2) * t52 + qJDD(2) * t34) * pkin(1);
t33 = sin(qJ(4));
t36 = cos(qJ(4));
t31 = qJD(2) + qJD(3);
t53 = qJD(2) * pkin(1);
t48 = t37 * t53;
t16 = t31 * pkin(2) + t48;
t49 = t34 * t53;
t7 = t33 * t16 + t36 * t49;
t14 = -qJD(3) * t49 + qJDD(2) * t58;
t8 = t30 * pkin(2) + t14;
t3 = -qJD(4) * t7 - t33 * t15 + t36 * t8;
t57 = t3 * mrSges(5,1) + Ifges(5,3) * t25;
t56 = t33 * t34;
t55 = t34 * t36;
t32 = qJ(2) + qJ(3);
t29 = qJ(4) + t32;
t20 = sin(t29);
t17 = t20 * mrSges(5,2);
t27 = sin(t32);
t54 = -t27 * mrSges(4,2) - t17;
t51 = qJD(4) * t33;
t50 = qJD(4) * t36;
t46 = m(5) * pkin(2) + mrSges(4,1);
t6 = t36 * t16 - t33 * t49;
t2 = qJD(4) * t6 + t36 * t15 + t33 * t8;
t21 = cos(t29);
t45 = g(3) * t21 - t2;
t44 = t14 * mrSges(4,1) + Ifges(4,3) * t30 + t57;
t43 = g(1) * t21 + g(3) * t20;
t28 = cos(t32);
t42 = g(1) * t28 + g(3) * t27;
t35 = sin(qJ(2));
t38 = cos(qJ(2));
t41 = g(1) * t38 + g(3) * t35;
t40 = -t33 * t37 - t55;
t39 = t36 * t37 - t56;
t26 = qJD(4) + t31;
t22 = pkin(2) + t58;
t12 = pkin(1) * t55 + t33 * t22;
t11 = -pkin(1) * t56 + t36 * t22;
t10 = t39 * t53;
t9 = t40 * t53;
t5 = -t22 * t51 + (t40 * qJD(3) - t34 * t50) * pkin(1);
t4 = t22 * t50 + (t39 * qJD(3) - t34 * t51) * pkin(1);
t1 = [(g(2) + qJDD(1)) * (m(2) + m(3) + m(4) + m(5)); Ifges(3,3) * qJDD(2) - t15 * mrSges(4,2) + m(5) * (t3 * t11 + t2 * t12 + t7 * t4 + t6 * t5) + (t11 * t25 + t5 * t26) * mrSges(5,1) + (t35 * mrSges(3,1) + t20 * mrSges(5,1) + t38 * mrSges(3,2) + t28 * mrSges(4,2) + t46 * t27) * g(3) + (t38 * mrSges(3,1) + t21 * mrSges(5,1) - t35 * mrSges(3,2) + t46 * t28 + t54) * g(1) + (-t12 * t25 - t4 * t26 + t45) * mrSges(5,2) + ((-t34 * t30 - t31 * t52) * mrSges(4,2) + (-qJD(3) * t34 * t31 + t37 * t30) * mrSges(4,1) + t41 * m(5) + (t14 * t37 + t15 * t34 + t41) * m(4)) * pkin(1) + t44; -m(5) * (t7 * t10 + t6 * t9) + t54 * g(1) + (-t9 * t26 + t43) * mrSges(5,1) + (t10 * t26 + t45) * mrSges(5,2) + (t31 * t49 + t42) * mrSges(4,1) + (g(3) * t28 + t31 * t48 - t15) * mrSges(4,2) + ((-t33 * t25 - t26 * t50) * mrSges(5,2) + (t36 * t25 - t26 * t51) * mrSges(5,1) + (t2 * t33 + t3 * t36 + t7 * t50 - t6 * t51 + t42) * m(5)) * pkin(2) + t44; -g(1) * t17 + (t7 * t26 + t43) * mrSges(5,1) + (t6 * t26 + t45) * mrSges(5,2) + t57;];
tau  = t1;
