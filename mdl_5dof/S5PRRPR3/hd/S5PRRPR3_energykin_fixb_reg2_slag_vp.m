% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRPR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:39
% EndTime: 2019-12-05 16:19:39
% DurationCPUTime: 0.16s
% Computational Cost: add. (228->40), mult. (588->99), div. (0->0), fcn. (384->6), ass. (0->33)
t34 = qJD(2) ^ 2;
t41 = t34 / 0.2e1;
t40 = cos(qJ(5));
t33 = cos(qJ(3));
t39 = t33 * t34;
t26 = t33 * qJD(1);
t32 = sin(qJ(3));
t37 = qJD(2) * t32;
t14 = qJD(3) * pkin(3) + t26 + (-pkin(6) - qJ(4)) * t37;
t36 = qJD(2) * t33;
t21 = pkin(6) * t36 + t32 * qJD(1);
t16 = qJ(4) * t36 + t21;
t30 = sin(pkin(9));
t38 = cos(pkin(9));
t6 = t30 * t14 + t38 * t16;
t35 = qJD(2) * qJD(3);
t5 = t38 * t14 - t30 * t16;
t22 = qJD(4) + (-pkin(3) * t33 - pkin(2)) * qJD(2);
t31 = sin(qJ(5));
t29 = qJD(1) ^ 2 / 0.2e1;
t28 = qJD(3) ^ 2 / 0.2e1;
t27 = qJD(3) + qJD(5);
t20 = -pkin(6) * t37 + t26;
t19 = (t30 * t33 + t32 * t38) * qJD(2);
t17 = t30 * t37 - t36 * t38;
t10 = t17 * pkin(4) + t22;
t9 = -t31 * t17 + t19 * t40;
t7 = t17 * t40 + t31 * t19;
t4 = -t17 * pkin(7) + t6;
t3 = qJD(3) * pkin(4) - t19 * pkin(7) + t5;
t2 = t31 * t3 + t4 * t40;
t1 = t3 * t40 - t31 * t4;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, 0, 0, 0, 0, t41, 0, 0, 0, t29, t32 ^ 2 * t41, t32 * t39, t32 * t35, t33 ^ 2 * t41, t33 * t35, t28, pkin(2) * t39 + t20 * qJD(3), -t34 * pkin(2) * t32 - t21 * qJD(3), (-t20 * t32 + t21 * t33) * qJD(2), t21 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1 + pkin(2) ^ 2 * t41, t19 ^ 2 / 0.2e1, -t19 * t17, t19 * qJD(3), t17 ^ 2 / 0.2e1, -t17 * qJD(3), t28, t5 * qJD(3) + t22 * t17, -t6 * qJD(3) + t22 * t19, -t6 * t17 - t5 * t19, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t9 ^ 2 / 0.2e1, -t9 * t7, t9 * t27, t7 ^ 2 / 0.2e1, -t7 * t27, t27 ^ 2 / 0.2e1, t1 * t27 + t10 * t7, t10 * t9 - t2 * t27, -t1 * t9 - t2 * t7, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg = t8;
