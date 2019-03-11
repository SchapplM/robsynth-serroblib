% Calculate vector of inverse dynamics joint torques for
% S4PRRR1
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-03-08 18:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR1_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR1_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR1_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR1_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR1_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:25:14
% EndTime: 2019-03-08 18:25:15
% DurationCPUTime: 0.28s
% Computational Cost: add. (350->95), mult. (594->124), div. (0->0), fcn. (272->10), ass. (0->53)
t59 = m(5) * pkin(3);
t34 = pkin(7) + qJ(2);
t32 = qJ(3) + t34;
t24 = qJ(4) + t32;
t19 = sin(t24);
t58 = g(1) * t19;
t22 = sin(t32);
t57 = g(1) * t22;
t39 = cos(qJ(3));
t56 = t39 * pkin(2);
t33 = qJDD(2) + qJDD(3);
t28 = qJDD(4) + t33;
t37 = sin(qJ(3));
t50 = qJD(3) * t39;
t15 = (qJD(2) * t50 + qJDD(2) * t37) * pkin(2);
t36 = sin(qJ(4));
t38 = cos(qJ(4));
t35 = qJD(2) + qJD(3);
t51 = pkin(2) * qJD(2);
t46 = t39 * t51;
t16 = t35 * pkin(3) + t46;
t47 = t37 * t51;
t7 = t36 * t16 + t38 * t47;
t14 = -qJD(3) * t47 + qJDD(2) * t56;
t8 = t33 * pkin(3) + t14;
t3 = -qJD(4) * t7 - t36 * t15 + t38 * t8;
t55 = t3 * mrSges(5,1) + Ifges(5,3) * t28;
t54 = t36 * t37;
t53 = t37 * t38;
t20 = cos(t24);
t17 = t20 * mrSges(5,1);
t23 = cos(t32);
t52 = -t23 * mrSges(4,1) - t17;
t49 = qJD(4) * t36;
t48 = qJD(4) * t38;
t44 = t14 * mrSges(4,1) + Ifges(4,3) * t33 + t55;
t29 = sin(t34);
t30 = cos(t34);
t43 = g(1) * t29 - g(2) * t30;
t42 = -t36 * t39 - t53;
t41 = t38 * t39 - t54;
t6 = t38 * t16 - t36 * t47;
t2 = qJD(4) * t6 + t38 * t15 + t36 * t8;
t40 = g(1) * t20 + g(2) * t19 - t2;
t31 = qJD(4) + t35;
t25 = pkin(3) + t56;
t12 = pkin(2) * t53 + t36 * t25;
t11 = -pkin(2) * t54 + t38 * t25;
t10 = t41 * t51;
t9 = t42 * t51;
t5 = -t25 * t49 + (t42 * qJD(3) - t37 * t48) * pkin(2);
t4 = t25 * t48 + (t41 * qJD(3) - t37 * t49) * pkin(2);
t1 = [(-g(3) + qJDD(1)) * (m(2) + m(3) + m(4) + m(5)); Ifges(3,3) * qJDD(2) - t15 * mrSges(4,2) + m(5) * (t3 * t11 + t2 * t12 + t7 * t4 + t6 * t5) + (t11 * t28 + t5 * t31) * mrSges(5,1) + (-t30 * mrSges(3,1) + t29 * mrSges(3,2) + t22 * mrSges(4,2) - t23 * t59 + t52) * g(2) + (t29 * mrSges(3,1) + t19 * mrSges(5,1) + t30 * mrSges(3,2) + t23 * mrSges(4,2) + (mrSges(4,1) + t59) * t22) * g(1) + (-t12 * t28 - t4 * t31 + t40) * mrSges(5,2) + ((-t37 * t33 - t35 * t50) * mrSges(4,2) + (-qJD(3) * t37 * t35 + t39 * t33) * mrSges(4,1) + t43 * m(5) + (t14 * t39 + t15 * t37 + t43) * m(4)) * pkin(2) + t44; -m(5) * (t7 * t10 + t6 * t9) + t52 * g(2) + (-t9 * t31 + t58) * mrSges(5,1) + (t10 * t31 + t40) * mrSges(5,2) + (t35 * t47 + t57) * mrSges(4,1) + (g(1) * t23 + g(2) * t22 + t35 * t46 - t15) * mrSges(4,2) + ((-t36 * t28 - t31 * t48) * mrSges(5,2) + (t38 * t28 - t31 * t49) * mrSges(5,1) + (-g(2) * t23 + t2 * t36 + t3 * t38 + t7 * t48 - t6 * t49 + t57) * m(5)) * pkin(3) + t44; -g(2) * t17 + (t7 * t31 + t58) * mrSges(5,1) + (t6 * t31 + t40) * mrSges(5,2) + t55;];
tau  = t1;
