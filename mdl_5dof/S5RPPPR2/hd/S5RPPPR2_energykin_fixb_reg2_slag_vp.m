% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPPR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 08:59:31
% EndTime: 2022-01-23 08:59:31
% DurationCPUTime: 0.21s
% Computational Cost: add. (319->50), mult. (897->124), div. (0->0), fcn. (581->8), ass. (0->44)
t35 = sin(pkin(9));
t38 = cos(pkin(9));
t40 = cos(pkin(7));
t37 = sin(pkin(7));
t39 = cos(pkin(8));
t52 = t37 * t39;
t20 = (t35 * t52 + t38 * t40) * qJD(1);
t42 = qJD(1) ^ 2;
t55 = t42 / 0.2e1;
t23 = qJD(2) + (-pkin(2) * t40 - qJ(3) * t37 - pkin(1)) * qJD(1);
t36 = sin(pkin(8));
t49 = qJ(2) * qJD(1);
t45 = t40 * t49;
t15 = t36 * t23 + t39 * t45;
t50 = qJD(1) * t40;
t10 = -qJ(4) * t50 + t15;
t28 = t37 * t49 + qJD(3);
t51 = qJD(1) * t37;
t17 = (pkin(3) * t36 - qJ(4) * t39) * t51 + t28;
t7 = t38 * t10 + t35 * t17;
t54 = cos(qJ(5));
t33 = t37 ^ 2;
t53 = t33 * t42;
t48 = t37 * t42 * t40;
t47 = t36 * t51;
t46 = t53 / 0.2e1;
t14 = t39 * t23 - t36 * t45;
t6 = -t35 * t10 + t38 * t17;
t9 = pkin(3) * t50 + qJD(4) - t14;
t41 = sin(qJ(5));
t34 = t40 ^ 2;
t32 = -qJD(1) * pkin(1) + qJD(2);
t29 = t34 * t55;
t25 = t36 ^ 2 * t46;
t22 = (-t35 * t40 + t38 * t52) * qJD(1);
t18 = qJD(5) + t20;
t13 = t54 * t22 + t41 * t47;
t11 = t41 * t22 - t54 * t47;
t5 = pkin(6) * t47 + t7;
t4 = -pkin(4) * t47 - t6;
t3 = t20 * pkin(4) - t22 * pkin(6) + t9;
t2 = t41 * t3 + t54 * t5;
t1 = t54 * t3 - t41 * t5;
t8 = [0, 0, 0, 0, 0, t55, 0, 0, 0, 0, t46, t48, 0, t29, 0, 0, -t32 * t50, t32 * t51, (t33 + t34) * t42 * qJ(2), t32 ^ 2 / 0.2e1 + (t34 / 0.2e1 + t33 / 0.2e1) * qJ(2) ^ 2 * t42, t39 ^ 2 * t46, -t39 * t36 * t53, -t39 * t48, t25, t36 * t48, t29, (t28 * t36 * t37 - t14 * t40) * qJD(1), (t15 * t40 + t28 * t52) * qJD(1), (-t14 * t39 - t15 * t36) * t51, t15 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * t47, t20 ^ 2 / 0.2e1, -t20 * t47, t25, t9 * t20 + t6 * t47, t9 * t22 - t7 * t47, -t7 * t20 - t6 * t22, t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t13 * t11, t13 * t18, t11 ^ 2 / 0.2e1, -t11 * t18, t18 ^ 2 / 0.2e1, t1 * t18 + t4 * t11, t4 * t13 - t2 * t18, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1;];
T_reg = t8;
