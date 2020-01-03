% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRPR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR5_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:22:58
% EndTime: 2019-12-31 16:22:59
% DurationCPUTime: 0.38s
% Computational Cost: add. (323->95), mult. (879->151), div. (0->0), fcn. (532->6), ass. (0->48)
t52 = -Ifges(5,1) / 0.2e1;
t32 = cos(qJ(4));
t41 = qJD(2) * t32;
t27 = Ifges(5,4) * t41;
t51 = -t27 / 0.2e1;
t30 = sin(qJ(4));
t33 = cos(qJ(2));
t24 = qJD(2) * pkin(2) + t33 * qJD(1);
t28 = sin(pkin(7));
t29 = cos(pkin(7));
t31 = sin(qJ(2));
t43 = qJD(1) * t31;
t8 = t28 * t24 + t29 * t43;
t6 = qJD(2) * pkin(5) + t8;
t3 = t32 * qJD(3) - t30 * t6;
t50 = qJD(4) * t3;
t4 = t30 * qJD(3) + t32 * t6;
t49 = t32 * t4;
t18 = t28 * t31 - t29 * t33;
t19 = t28 * t33 + t29 * t31;
t12 = t19 * qJD(2);
t9 = qJD(1) * t12;
t48 = t9 * t18;
t47 = Ifges(5,4) * t30;
t45 = Ifges(5,5) * qJD(4);
t44 = Ifges(5,6) * qJD(4);
t42 = qJD(2) * t30;
t40 = t45 / 0.2e1;
t39 = -t44 / 0.2e1;
t38 = -t3 * t30 + t49;
t22 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t42;
t23 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t41;
t37 = -t30 * t22 + t32 * t23;
t7 = t29 * t24 - t28 * t43;
t5 = -qJD(2) * pkin(3) - t7;
t36 = t3 * mrSges(5,3) + t42 * t52 + t51 - t45 / 0.2e1 - t5 * mrSges(5,2);
t35 = t4 * mrSges(5,3) + t44 / 0.2e1 + (Ifges(5,2) * t32 + t47) * qJD(2) / 0.2e1 - t5 * mrSges(5,1);
t14 = t18 * qJD(2);
t26 = -t29 * pkin(2) - pkin(3);
t25 = t28 * pkin(2) + pkin(5);
t20 = (-t32 * mrSges(5,1) + t30 * mrSges(5,2)) * qJD(2);
t17 = (mrSges(5,1) * t30 + mrSges(5,2) * t32) * qJD(4) * qJD(2);
t13 = t18 * qJD(1);
t11 = t19 * qJD(1);
t10 = qJD(1) * t14;
t2 = -qJD(4) * t4 + t30 * t10;
t1 = -t32 * t10 + t50;
t15 = [t18 * t17 + (-t31 * mrSges(3,1) - t33 * mrSges(3,2)) * qJD(2) ^ 2 + (-qJD(2) * mrSges(4,1) + t20) * t12 + (-t22 * t32 - t23 * t30) * t19 * qJD(4) - (-qJD(2) * mrSges(4,2) + t37) * t14 + m(4) * (-t10 * t19 - t7 * t12 - t8 * t14 + t48) + m(5) * (t5 * t12 + t48 - t38 * t14 + (t1 * t32 - t2 * t30 + (-t3 * t32 - t4 * t30) * qJD(4)) * t19); -t11 * t20 + t26 * t17 + (-t13 * qJD(2) + t10) * mrSges(4,2) + (t11 * qJD(2) - t9) * mrSges(4,1) + (-t9 * mrSges(5,1) + t1 * mrSges(5,3) + t13 * t23 + (-t25 * t22 + 0.3e1 / 0.2e1 * t27 + t40 - t36) * qJD(4)) * t32 - m(5) * (t5 * t11 - t13 * t49) + m(5) * (t9 * t26 + (t1 - t50) * t25 * t32) + (t9 * mrSges(5,2) + (t39 + (-m(5) * t4 - t23) * t25 + (-0.3e1 / 0.2e1 * t47 + (0.3e1 / 0.2e1 * Ifges(5,1) - 0.3e1 / 0.2e1 * Ifges(5,2)) * t32) * qJD(2) - t35) * qJD(4) + (-t25 * m(5) - mrSges(5,3)) * t2 - (t3 * m(5) + t22) * t13) * t30 + (t7 * t11 + t8 * t13 + (-t10 * t28 - t29 * t9) * pkin(2)) * m(4); m(5) * (t1 * t30 + t2 * t32) + (m(5) * t38 + (-t30 ^ 2 - t32 ^ 2) * qJD(2) * mrSges(5,3) + t37) * qJD(4); t2 * mrSges(5,1) - t1 * mrSges(5,2) + t4 * t22 - t3 * t23 + ((t40 + t51 + t36) * t32 + (t39 + (t47 / 0.2e1 + (t52 + Ifges(5,2) / 0.2e1) * t32) * qJD(2) + t35) * t30) * qJD(2);];
tauc = t15(:);
