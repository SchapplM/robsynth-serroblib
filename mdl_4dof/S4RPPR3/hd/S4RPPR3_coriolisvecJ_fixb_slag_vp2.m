% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPPR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR3_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:50
% EndTime: 2019-12-31 16:37:51
% DurationCPUTime: 0.36s
% Computational Cost: add. (418->88), mult. (1093->138), div. (0->0), fcn. (702->6), ass. (0->49)
t37 = sin(pkin(7));
t39 = cos(pkin(7));
t48 = qJD(1) * (t37 ^ 2 + t39 ^ 2);
t59 = mrSges(4,3) * t48;
t41 = sin(qJ(4));
t42 = cos(qJ(4));
t29 = -t41 * t37 + t42 * t39;
t24 = t29 * qJD(4);
t57 = t24 / 0.2e1;
t30 = t42 * t37 + t41 * t39;
t25 = t30 * qJD(4);
t56 = -t25 / 0.2e1;
t32 = sin(pkin(6)) * pkin(1) + qJ(3);
t55 = pkin(5) + t32;
t23 = t30 * qJD(1);
t54 = Ifges(5,4) * t23;
t18 = qJD(1) * t24;
t53 = t29 * t18;
t19 = qJD(1) * t25;
t52 = t30 * t19;
t31 = qJD(1) * t32;
t17 = t37 * qJD(2) + t39 * t31;
t50 = pkin(5) * qJD(1);
t49 = t19 * mrSges(5,1) + t18 * mrSges(5,2);
t34 = t39 * qJD(2);
t11 = t34 + (-t31 - t50) * t37;
t12 = t39 * t50 + t17;
t3 = t42 * t11 - t41 * t12;
t4 = t41 * t11 + t42 * t12;
t47 = -(-t37 * t31 + t34) * t37 + t17 * t39;
t26 = t55 * t37;
t27 = t55 * t39;
t9 = -t42 * t26 - t41 * t27;
t10 = -t41 * t26 + t42 * t27;
t46 = -cos(pkin(6)) * pkin(1) - pkin(2) - t39 * pkin(3);
t45 = t30 * qJD(3);
t44 = t29 * qJD(3);
t22 = t29 * qJD(1);
t21 = t46 * qJD(1) + qJD(3);
t20 = Ifges(5,4) * t22;
t14 = qJD(4) * mrSges(5,1) - t23 * mrSges(5,3);
t13 = -qJD(4) * mrSges(5,2) + t22 * mrSges(5,3);
t8 = Ifges(5,1) * t23 + Ifges(5,5) * qJD(4) + t20;
t7 = Ifges(5,2) * t22 + Ifges(5,6) * qJD(4) + t54;
t6 = -t10 * qJD(4) - t45;
t5 = t9 * qJD(4) + t44;
t2 = -qJD(1) * t45 - t4 * qJD(4);
t1 = qJD(1) * t44 + t3 * qJD(4);
t15 = [t21 * (t25 * mrSges(5,1) + t24 * mrSges(5,2)) + m(5) * (t1 * t10 + t2 * t9 + t3 * t6 + t4 * t5) + t5 * t13 + t6 * t14 + qJD(4) * (Ifges(5,5) * t24 - Ifges(5,6) * t25) / 0.2e1 + t8 * t57 + t7 * t56 + t46 * t49 + (-t19 * t29 + t22 * t56) * Ifges(5,2) + (t18 * t30 + t23 * t57) * Ifges(5,1) + (m(4) * (t32 * t48 + t47) + 0.2e1 * t59) * qJD(3) + (t1 * t29 - t10 * t19 - t9 * t18 - t2 * t30 - t3 * t24 - t4 * t25) * mrSges(5,3) + (t22 * t57 + t23 * t56 - t52 + t53) * Ifges(5,4); m(5) * (t1 * t30 + t2 * t29 + t4 * t24 - t3 * t25) + t24 * t13 - t25 * t14 + (-t52 - t53) * mrSges(5,3); -t22 * t13 + t23 * t14 - m(5) * (t4 * t22 - t3 * t23) + t49 + (-m(4) * t47 - t59) * qJD(1); Ifges(5,5) * t18 - Ifges(5,6) * t19 - t1 * mrSges(5,2) + t2 * mrSges(5,1) - t21 * (t23 * mrSges(5,1) + t22 * mrSges(5,2)) - t23 * (Ifges(5,1) * t22 - t54) / 0.2e1 + t23 * t7 / 0.2e1 - qJD(4) * (Ifges(5,5) * t22 - Ifges(5,6) * t23) / 0.2e1 - t3 * t13 + t4 * t14 + (t3 * t22 + t4 * t23) * mrSges(5,3) - (-Ifges(5,2) * t23 + t20 + t8) * t22 / 0.2e1;];
tauc = t15(:);
