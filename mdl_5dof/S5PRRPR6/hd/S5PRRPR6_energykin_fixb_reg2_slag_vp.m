% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRPR6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR6_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:32:45
% EndTime: 2019-12-05 16:32:45
% DurationCPUTime: 0.19s
% Computational Cost: add. (272->45), mult. (681->116), div. (0->0), fcn. (466->10), ass. (0->41)
t41 = qJD(2) ^ 2;
t52 = t41 / 0.2e1;
t51 = cos(qJ(5));
t38 = sin(qJ(2));
t34 = sin(pkin(5));
t49 = qJD(1) * t34;
t25 = qJD(2) * pkin(7) + t38 * t49;
t37 = sin(qJ(3));
t39 = cos(qJ(3));
t35 = cos(pkin(5));
t48 = qJD(1) * t35;
t17 = t39 * t25 + t37 * t48;
t15 = qJD(3) * qJ(4) + t17;
t40 = cos(qJ(2));
t44 = t40 * t49;
t18 = -t44 + (-pkin(3) * t39 - qJ(4) * t37 - pkin(2)) * qJD(2);
t33 = sin(pkin(10));
t50 = cos(pkin(10));
t6 = t50 * t15 + t33 * t18;
t47 = qJD(2) * t37;
t46 = t39 * qJD(2);
t45 = qJD(2) * qJD(3);
t43 = qJD(2) * t49;
t5 = -t33 * t15 + t50 * t18;
t16 = -t37 * t25 + t39 * t48;
t12 = -qJD(3) * pkin(3) + qJD(4) - t16;
t42 = qJD(1) ^ 2;
t36 = sin(qJ(5));
t30 = t39 ^ 2 * t52;
t28 = -qJD(5) + t46;
t26 = -qJD(2) * pkin(2) - t44;
t24 = t33 * qJD(3) + t50 * t47;
t22 = -t50 * qJD(3) + t33 * t47;
t10 = -t36 * t22 + t51 * t24;
t8 = t51 * t22 + t36 * t24;
t7 = t22 * pkin(4) + t12;
t4 = -t22 * pkin(8) + t6;
t3 = -pkin(4) * t46 - t24 * pkin(8) + t5;
t2 = t36 * t3 + t51 * t4;
t1 = t51 * t3 - t36 * t4;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t42 / 0.2e1, 0, 0, 0, 0, 0, t52, t40 * t43, -t38 * t43, 0, (t35 ^ 2 / 0.2e1 + (t38 ^ 2 / 0.2e1 + t40 ^ 2 / 0.2e1) * t34 ^ 2) * t42, t37 ^ 2 * t52, t37 * t41 * t39, t37 * t45, t30, t39 * t45, qJD(3) ^ 2 / 0.2e1, t16 * qJD(3) - t26 * t46, -t17 * qJD(3) + t26 * t47, (-t16 * t37 + t17 * t39) * qJD(2), t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1, t24 ^ 2 / 0.2e1, -t24 * t22, -t24 * t46, t22 ^ 2 / 0.2e1, t22 * t46, t30, t12 * t22 - t5 * t46, t12 * t24 + t6 * t46, -t6 * t22 - t5 * t24, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, -t10 * t28, t8 ^ 2 / 0.2e1, t8 * t28, t28 ^ 2 / 0.2e1, -t1 * t28 + t7 * t8, t7 * t10 + t2 * t28, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg = t9;
