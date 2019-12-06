% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPRR6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR6_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:57:52
% EndTime: 2019-12-05 15:57:52
% DurationCPUTime: 0.18s
% Computational Cost: add. (242->44), mult. (642->113), div. (0->0), fcn. (460->10), ass. (0->40)
t38 = qJD(2) ^ 2;
t48 = t38 / 0.2e1;
t36 = sin(qJ(2));
t31 = sin(pkin(5));
t45 = qJD(1) * t31;
t24 = qJD(2) * qJ(3) + t36 * t45;
t32 = cos(pkin(10));
t33 = cos(pkin(5));
t44 = qJD(1) * t33;
t26 = t32 * t44;
t30 = sin(pkin(10));
t10 = t26 + (-pkin(7) * qJD(2) - t24) * t30;
t16 = t32 * t24 + t30 * t44;
t42 = qJD(2) * t32;
t11 = pkin(7) * t42 + t16;
t35 = sin(qJ(4));
t47 = cos(qJ(4));
t6 = t35 * t10 + t47 * t11;
t46 = cos(qJ(5));
t43 = qJD(2) * t30;
t41 = qJD(2) * t45;
t19 = t35 * t43 - t47 * t42;
t37 = cos(qJ(2));
t40 = -t37 * t45 + qJD(3);
t5 = t47 * t10 - t35 * t11;
t17 = (-pkin(3) * t32 - pkin(2)) * qJD(2) + t40;
t39 = qJD(1) ^ 2;
t34 = sin(qJ(5));
t23 = -qJD(2) * pkin(2) + t40;
t21 = (t47 * t30 + t32 * t35) * qJD(2);
t18 = qJD(5) + t19;
t15 = -t30 * t24 + t26;
t14 = t34 * qJD(4) + t46 * t21;
t12 = -t46 * qJD(4) + t34 * t21;
t7 = t19 * pkin(4) - t21 * pkin(8) + t17;
t4 = qJD(4) * pkin(8) + t6;
t3 = -qJD(4) * pkin(4) - t5;
t2 = t34 * t7 + t46 * t4;
t1 = -t34 * t4 + t46 * t7;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t39 / 0.2e1, 0, 0, 0, 0, 0, t48, t37 * t41, -t36 * t41, 0, (t33 ^ 2 / 0.2e1 + (t36 ^ 2 / 0.2e1 + t37 ^ 2 / 0.2e1) * t31 ^ 2) * t39, t30 ^ 2 * t48, t30 * t38 * t32, 0, t32 ^ 2 * t48, 0, 0, -t23 * t42, t23 * t43, (-t15 * t30 + t16 * t32) * qJD(2), t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1, t21 ^ 2 / 0.2e1, -t21 * t19, t21 * qJD(4), t19 ^ 2 / 0.2e1, -t19 * qJD(4), qJD(4) ^ 2 / 0.2e1, t5 * qJD(4) + t17 * t19, -t6 * qJD(4) + t17 * t21, -t6 * t19 - t5 * t21, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t14 ^ 2 / 0.2e1, -t14 * t12, t14 * t18, t12 ^ 2 / 0.2e1, -t12 * t18, t18 ^ 2 / 0.2e1, t1 * t18 + t3 * t12, t3 * t14 - t2 * t18, -t1 * t14 - t2 * t12, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t8;
