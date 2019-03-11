% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRP2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP2_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:01:15
% EndTime: 2019-03-09 06:01:15
% DurationCPUTime: 0.16s
% Computational Cost: add. (514->56), mult. (1123->128), div. (0->0), fcn. (708->8), ass. (0->45)
t53 = qJD(1) ^ 2;
t45 = t53 / 0.2e1;
t46 = sin(pkin(10));
t37 = (pkin(1) * t46 + pkin(7)) * qJD(1);
t50 = sin(qJ(3));
t52 = cos(qJ(3));
t29 = t50 * qJD(2) + t52 * t37;
t26 = qJD(3) * pkin(8) + t29;
t47 = cos(pkin(10));
t55 = -pkin(1) * t47 - pkin(2);
t27 = (-pkin(3) * t52 - pkin(8) * t50 + t55) * qJD(1);
t49 = sin(qJ(4));
t51 = cos(qJ(4));
t13 = t51 * t26 + t49 * t27;
t58 = qJD(1) * t50;
t33 = -t51 * qJD(3) + t49 * t58;
t11 = -t33 * pkin(9) + t13;
t48 = sin(qJ(5));
t59 = cos(qJ(5));
t12 = -t49 * t26 + t51 * t27;
t35 = t49 * qJD(3) + t51 * t58;
t57 = t52 * qJD(1);
t41 = -qJD(4) + t57;
t7 = -t41 * pkin(4) - t35 * pkin(9) + t12;
t4 = t59 * t11 + t48 * t7;
t60 = pkin(1) * t53;
t56 = qJD(1) * qJD(3);
t3 = -t48 * t11 + t59 * t7;
t28 = t52 * qJD(2) - t50 * t37;
t25 = -qJD(3) * pkin(3) - t28;
t16 = t33 * pkin(4) + t25;
t39 = -qJD(5) + t41;
t38 = t55 * qJD(1);
t36 = t39 ^ 2 / 0.2e1;
t21 = -t48 * t33 + t59 * t35;
t19 = t59 * t33 + t48 * t35;
t18 = t21 ^ 2 / 0.2e1;
t17 = t19 ^ 2 / 0.2e1;
t15 = t21 * t39;
t14 = t19 * t39;
t10 = t21 * t19;
t8 = t19 * pkin(5) + qJD(6) + t16;
t2 = -t19 * qJ(6) + t4;
t1 = -t39 * pkin(5) - t21 * qJ(6) + t3;
t5 = [0, 0, 0, 0, 0, t45, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t47 * t60, -t46 * t60, 0, qJD(2) ^ 2 / 0.2e1 + (t46 ^ 2 / 0.2e1 + t47 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t53, t50 ^ 2 * t45, t52 * t53 * t50, t50 * t56, t52 ^ 2 * t45, t52 * t56, qJD(3) ^ 2 / 0.2e1, t28 * qJD(3) - t38 * t57, -t29 * qJD(3) + t38 * t58 (-t28 * t50 + t29 * t52) * qJD(1), t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1, t35 ^ 2 / 0.2e1, -t35 * t33, -t35 * t41, t33 ^ 2 / 0.2e1, t33 * t41, t41 ^ 2 / 0.2e1, -t12 * t41 + t25 * t33, t13 * t41 + t25 * t35, -t12 * t35 - t13 * t33, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1, t18, -t10, -t15, t17, t14, t36, t16 * t19 - t3 * t39, t16 * t21 + t4 * t39, -t4 * t19 - t3 * t21, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1, t18, -t10, -t15, t17, t14, t36, -t1 * t39 + t8 * t19, t2 * t39 + t8 * t21, -t1 * t21 - t2 * t19, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1;];
T_reg  = t5;
