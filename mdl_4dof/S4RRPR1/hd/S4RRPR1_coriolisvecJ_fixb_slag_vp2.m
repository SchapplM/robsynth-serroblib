% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR1_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:34:54
% EndTime: 2019-03-08 18:34:55
% DurationCPUTime: 0.26s
% Computational Cost: add. (417->68), mult. (1056->106), div. (0->0), fcn. (582->6), ass. (0->48)
t31 = sin(pkin(7));
t36 = cos(qJ(2));
t32 = cos(pkin(7));
t34 = sin(qJ(2));
t48 = t32 * t34;
t38 = pkin(1) * (-t31 * t36 - t48);
t19 = qJD(1) * t38;
t49 = t31 * t34;
t37 = pkin(1) * (t32 * t36 - t49);
t21 = qJD(1) * t37;
t33 = sin(qJ(4));
t35 = cos(qJ(4));
t27 = pkin(2) * t32 + pkin(3);
t55 = pkin(2) * t31;
t39 = t27 * t35 - t33 * t55;
t58 = t39 * qJD(4) - t19 * t33 - t21 * t35;
t40 = t27 * t33 + t35 * t55;
t57 = -t40 * qJD(4) - t19 * t35 + t21 * t33;
t56 = mrSges(3,1) * t34 + mrSges(3,2) * t36;
t20 = qJD(2) * t38;
t14 = qJD(1) * t20;
t22 = qJD(2) * t37;
t15 = qJD(1) * t22;
t30 = qJD(1) + qJD(2);
t47 = pkin(1) * qJD(1);
t25 = t30 * pkin(2) + t36 * t47;
t44 = t34 * t47;
t11 = t32 * t25 - t31 * t44;
t10 = pkin(3) * t30 + t11;
t12 = t25 * t31 + t32 * t44;
t7 = t10 * t33 + t12 * t35;
t3 = -t7 * qJD(4) + t14 * t35 - t15 * t33;
t1 = t3 * mrSges(5,1);
t54 = t14 * mrSges(4,1) + t1;
t29 = qJD(4) + t30;
t51 = t29 * mrSges(5,1);
t50 = t30 * mrSges(4,1);
t28 = pkin(1) * t36 + pkin(2);
t43 = -pkin(1) * t49 + t32 * t28;
t6 = t10 * t35 - t12 * t33;
t18 = pkin(3) + t43;
t23 = pkin(1) * t48 + t28 * t31;
t42 = t18 * t35 - t23 * t33;
t41 = t18 * t33 + t23 * t35;
t5 = -t41 * qJD(4) + t20 * t35 - t22 * t33;
t4 = t42 * qJD(4) + t20 * t33 + t22 * t35;
t2 = t6 * qJD(4) + t14 * t33 + t15 * t35;
t8 = [t20 * t50 + t5 * t51 + (-t29 * t4 - t2) * mrSges(5,2) + (-t22 * t30 - t15) * mrSges(4,2) + m(5) * (t2 * t41 + t3 * t42 + t7 * t4 + t6 * t5) + m(4) * (t11 * t20 + t12 * t22 + t14 * t43 + t15 * t23) + t54 + t56 * pkin(1) * qJD(2) * (-qJD(1) - t30); -t19 * t50 - t2 * mrSges(5,2) + (t21 * t30 - t15) * mrSges(4,2) + (t57 * mrSges(5,1) - t58 * mrSges(5,2)) * t29 + t54 + t56 * t47 * (-qJD(2) + t30) + (t2 * t40 + t3 * t39 + t57 * t6 + t58 * t7) * m(5) + ((t14 * t32 + t15 * t31) * pkin(2) - t11 * t19 - t12 * t21) * m(4); 0; t7 * t51 + t1 + (t29 * t6 - t2) * mrSges(5,2);];
tauc  = t8(:);
