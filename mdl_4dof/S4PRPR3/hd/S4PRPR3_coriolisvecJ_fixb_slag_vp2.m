% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRPR3
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
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRPR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR3_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:48
% EndTime: 2019-12-31 16:20:49
% DurationCPUTime: 0.33s
% Computational Cost: add. (345->85), mult. (966->136), div. (0->0), fcn. (629->4), ass. (0->47)
t35 = sin(pkin(7));
t36 = cos(pkin(7));
t43 = (t35 ^ 2 + t36 ^ 2) * qJD(2);
t55 = mrSges(4,3) * t43;
t37 = sin(qJ(4));
t38 = cos(qJ(4));
t23 = -t37 * t35 + t38 * t36;
t21 = t23 * qJD(4);
t53 = t21 / 0.2e1;
t24 = t38 * t35 + t37 * t36;
t22 = t24 * qJD(4);
t52 = -t22 / 0.2e1;
t20 = t24 * qJD(2);
t51 = Ifges(5,4) * t20;
t14 = qJD(2) * t21;
t50 = t23 * t14;
t15 = qJD(2) * t22;
t49 = t24 * t15;
t48 = pkin(5) + qJ(3);
t46 = qJD(2) * qJ(3);
t26 = t35 * qJD(1) + t36 * t46;
t45 = -t36 * pkin(3) - pkin(2);
t28 = t48 * t35;
t44 = t15 * mrSges(5,1) + t14 * mrSges(5,2);
t32 = t36 * qJD(1);
t16 = -qJD(2) * t28 + t32;
t17 = t36 * qJD(2) * pkin(5) + t26;
t5 = t38 * t16 - t37 * t17;
t6 = t37 * t16 + t38 * t17;
t42 = -(-t35 * t46 + t32) * t35 + t26 * t36;
t29 = t48 * t36;
t9 = -t38 * t28 - t37 * t29;
t10 = -t37 * t28 + t38 * t29;
t41 = t24 * qJD(3);
t40 = t23 * qJD(3);
t27 = t45 * qJD(2) + qJD(3);
t19 = t23 * qJD(2);
t18 = Ifges(5,4) * t19;
t12 = qJD(4) * mrSges(5,1) - t20 * mrSges(5,3);
t11 = -qJD(4) * mrSges(5,2) + t19 * mrSges(5,3);
t8 = Ifges(5,1) * t20 + Ifges(5,5) * qJD(4) + t18;
t7 = Ifges(5,2) * t19 + Ifges(5,6) * qJD(4) + t51;
t4 = -t10 * qJD(4) - t41;
t3 = t9 * qJD(4) + t40;
t2 = -qJD(2) * t41 - t6 * qJD(4);
t1 = qJD(2) * t40 + t5 * qJD(4);
t13 = [t21 * t11 - t22 * t12 + m(5) * (t1 * t24 + t2 * t23 + t6 * t21 - t5 * t22) + (-t49 - t50) * mrSges(5,3); qJD(4) * (Ifges(5,5) * t21 - Ifges(5,6) * t22) / 0.2e1 + t8 * t53 + t7 * t52 + t4 * t12 + t3 * t11 + t45 * t44 + t27 * (t22 * mrSges(5,1) + t21 * mrSges(5,2)) + m(5) * (t1 * t10 + t2 * t9 + t6 * t3 + t5 * t4) + (-t15 * t23 + t19 * t52) * Ifges(5,2) + (t14 * t24 + t20 * t53) * Ifges(5,1) + (m(4) * (qJ(3) * t43 + t42) + 0.2e1 * t55) * qJD(3) + (t1 * t23 - t10 * t15 - t9 * t14 - t2 * t24 - t5 * t21 - t6 * t22) * mrSges(5,3) + (t19 * t53 + t20 * t52 - t49 + t50) * Ifges(5,4); -t19 * t11 + t20 * t12 - m(5) * (t6 * t19 - t5 * t20) + t44 + (-m(4) * t42 - t55) * qJD(2); Ifges(5,5) * t14 - Ifges(5,6) * t15 - t1 * mrSges(5,2) + t2 * mrSges(5,1) - t27 * (t20 * mrSges(5,1) + t19 * mrSges(5,2)) - t20 * (Ifges(5,1) * t19 - t51) / 0.2e1 + t20 * t7 / 0.2e1 - qJD(4) * (Ifges(5,5) * t19 - Ifges(5,6) * t20) / 0.2e1 - t5 * t11 + t6 * t12 + (t5 * t19 + t6 * t20) * mrSges(5,3) - (-Ifges(5,2) * t20 + t18 + t8) * t19 / 0.2e1;];
tauc = t13(:);
