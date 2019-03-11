% Calculate vector of inverse dynamics joint torques for
% S3RRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% qJDD [3x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% m_mdh [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [3x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S3RRR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_invdynJ_fixb_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRR1_invdynJ_fixb_slag_vp2: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RRR1_invdynJ_fixb_slag_vp2: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RRR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_invdynJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RRR1_invdynJ_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RRR1_invdynJ_fixb_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3RRR1_invdynJ_fixb_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:07:57
% EndTime: 2019-03-08 18:07:58
% DurationCPUTime: 0.28s
% Computational Cost: add. (311->90), mult. (586->123), div. (0->0), fcn. (272->10), ass. (0->52)
t57 = m(4) * pkin(2);
t32 = qJ(1) + qJ(2);
t29 = qJ(3) + t32;
t20 = sin(t29);
t56 = g(1) * t20;
t27 = sin(t32);
t55 = g(1) * t27;
t37 = cos(qJ(2));
t54 = t37 * pkin(1);
t30 = qJDD(1) + qJDD(2);
t25 = qJDD(3) + t30;
t34 = sin(qJ(2));
t48 = qJD(2) * t37;
t15 = (qJD(1) * t48 + qJDD(1) * t34) * pkin(1);
t33 = sin(qJ(3));
t36 = cos(qJ(3));
t31 = qJD(1) + qJD(2);
t49 = pkin(1) * qJD(1);
t44 = t37 * t49;
t16 = t31 * pkin(2) + t44;
t45 = t34 * t49;
t7 = t33 * t16 + t36 * t45;
t14 = -qJD(2) * t45 + qJDD(1) * t54;
t8 = t30 * pkin(2) + t14;
t3 = -qJD(3) * t7 - t33 * t15 + t36 * t8;
t53 = t3 * mrSges(4,1) + Ifges(4,3) * t25;
t52 = t33 * t34;
t51 = t34 * t36;
t21 = cos(t29);
t17 = t21 * mrSges(4,1);
t28 = cos(t32);
t50 = -t28 * mrSges(3,1) - t17;
t47 = qJD(3) * t33;
t46 = qJD(3) * t36;
t43 = t14 * mrSges(3,1) + Ifges(3,3) * t30 + t53;
t35 = sin(qJ(1));
t38 = cos(qJ(1));
t42 = g(1) * t35 - g(2) * t38;
t41 = -t33 * t37 - t51;
t40 = t36 * t37 - t52;
t6 = t36 * t16 - t33 * t45;
t2 = qJD(3) * t6 + t36 * t15 + t33 * t8;
t39 = g(1) * t21 + g(2) * t20 - t2;
t26 = qJD(3) + t31;
t22 = pkin(2) + t54;
t12 = pkin(1) * t51 + t33 * t22;
t11 = -pkin(1) * t52 + t36 * t22;
t10 = t40 * t49;
t9 = t41 * t49;
t5 = -t22 * t47 + (t41 * qJD(2) - t34 * t46) * pkin(1);
t4 = t22 * t46 + (t40 * qJD(2) - t34 * t47) * pkin(1);
t1 = [Ifges(2,3) * qJDD(1) - t15 * mrSges(3,2) + m(4) * (t3 * t11 + t2 * t12 + t7 * t4 + t6 * t5) + (t11 * t25 + t5 * t26) * mrSges(4,1) + (-t38 * mrSges(2,1) + t35 * mrSges(2,2) + t27 * mrSges(3,2) - t28 * t57 + t50) * g(2) + (t35 * mrSges(2,1) + t20 * mrSges(4,1) + t38 * mrSges(2,2) + t28 * mrSges(3,2) + (mrSges(3,1) + t57) * t27) * g(1) + (-t12 * t25 - t4 * t26 + t39) * mrSges(4,2) + ((-t34 * t30 - t31 * t48) * mrSges(3,2) + (-qJD(2) * t34 * t31 + t37 * t30) * mrSges(3,1) + t42 * m(4) + (t14 * t37 + t15 * t34 + t42) * m(3)) * pkin(1) + t43; -m(4) * (t7 * t10 + t6 * t9) + t50 * g(2) + (-t9 * t26 + t56) * mrSges(4,1) + (t10 * t26 + t39) * mrSges(4,2) + (t31 * t45 + t55) * mrSges(3,1) + (g(1) * t28 + g(2) * t27 + t31 * t44 - t15) * mrSges(3,2) + ((-t33 * t25 - t26 * t46) * mrSges(4,2) + (t36 * t25 - t26 * t47) * mrSges(4,1) + (-g(2) * t28 + t2 * t33 + t3 * t36 + t7 * t46 - t6 * t47 + t55) * m(4)) * pkin(2) + t43; -g(2) * t17 + (t7 * t26 + t56) * mrSges(4,1) + (t6 * t26 + t39) * mrSges(4,2) + t53;];
tau  = t1;
