% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:41
% EndTime: 2019-12-05 17:47:42
% DurationCPUTime: 0.18s
% Computational Cost: add. (275->40), mult. (599->98), div. (0->0), fcn. (340->6), ass. (0->34)
t37 = qJD(1) ^ 2;
t28 = t37 / 0.2e1;
t42 = cos(qJ(5));
t36 = cos(qJ(3));
t22 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t38 = -qJ(4) * qJD(1) + t22;
t14 = qJD(3) * pkin(3) + t36 * t38;
t35 = sin(qJ(3));
t16 = t38 * t35;
t32 = sin(pkin(8));
t33 = cos(pkin(8));
t6 = t32 * t14 + t33 * t16;
t41 = t37 * qJ(2);
t40 = qJD(3) * t22;
t39 = qJD(1) * qJD(3);
t20 = qJD(4) + (pkin(3) * t35 + qJ(2)) * qJD(1);
t5 = t33 * t14 - t32 * t16;
t34 = sin(qJ(5));
t31 = t36 ^ 2;
t30 = t35 ^ 2;
t27 = qJD(3) ^ 2 / 0.2e1;
t26 = qJD(3) + qJD(5);
t25 = qJ(2) ^ 2 * t28;
t24 = -pkin(1) * qJD(1) + qJD(2);
t19 = (-t32 * t35 + t33 * t36) * qJD(1);
t17 = (-t32 * t36 - t33 * t35) * qJD(1);
t10 = -t17 * pkin(4) + t20;
t9 = t34 * t17 + t42 * t19;
t7 = -t42 * t17 + t34 * t19;
t4 = t17 * pkin(7) + t6;
t3 = qJD(3) * pkin(4) - t19 * pkin(7) + t5;
t2 = t34 * t3 + t42 * t4;
t1 = t42 * t3 - t34 * t4;
t8 = [0, 0, 0, 0, 0, t28, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, 0, t24 * qJD(1), t41, t25 + t24 ^ 2 / 0.2e1, t31 * t28, -t36 * t37 * t35, t36 * t39, t30 * t28, -t35 * t39, t27, t35 * t41 + t36 * t40, -t35 * t40 + t36 * t41, (-t30 - t31) * t22 * qJD(1), t25 + (t30 / 0.2e1 + t31 / 0.2e1) * t22 ^ 2, t19 ^ 2 / 0.2e1, t19 * t17, t19 * qJD(3), t17 ^ 2 / 0.2e1, t17 * qJD(3), t27, t5 * qJD(3) - t20 * t17, -t6 * qJD(3) + t20 * t19, t6 * t17 - t5 * t19, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1, t9 ^ 2 / 0.2e1, -t9 * t7, t9 * t26, t7 ^ 2 / 0.2e1, -t7 * t26, t26 ^ 2 / 0.2e1, t1 * t26 + t10 * t7, t10 * t9 - t2 * t26, -t1 * t9 - t2 * t7, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg = t8;
