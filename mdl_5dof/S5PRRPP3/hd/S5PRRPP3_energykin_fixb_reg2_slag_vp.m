% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRPP3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP3_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:13:31
% EndTime: 2019-12-05 16:13:31
% DurationCPUTime: 0.15s
% Computational Cost: add. (161->38), mult. (393->90), div. (0->0), fcn. (217->6), ass. (0->36)
t28 = qJD(2) ^ 2;
t41 = t28 / 0.2e1;
t25 = sin(qJ(2));
t17 = qJD(2) * pkin(6) + t25 * qJD(1);
t26 = cos(qJ(3));
t10 = qJD(3) * qJ(4) + t26 * t17;
t23 = sin(pkin(8));
t39 = cos(pkin(8));
t24 = sin(qJ(3));
t27 = cos(qJ(2));
t35 = t27 * qJD(1);
t7 = -t35 + (-pkin(3) * t26 - qJ(4) * t24 - pkin(2)) * qJD(2);
t5 = t39 * t10 + t23 * t7;
t38 = qJD(2) * t24;
t12 = -t39 * qJD(3) + t23 * t38;
t14 = t23 * qJD(3) + t39 * t38;
t40 = t14 * t12;
t37 = qJD(2) * t26;
t36 = qJD(3) * t17;
t34 = t12 ^ 2 / 0.2e1;
t33 = qJD(1) * qJD(2);
t32 = qJD(2) * qJD(3);
t31 = t12 * t37;
t30 = t14 * t37;
t9 = -qJD(3) * pkin(3) + t24 * t17 + qJD(4);
t4 = -t23 * t10 + t39 * t7;
t29 = qJD(1) ^ 2;
t22 = t26 ^ 2;
t21 = t24 ^ 2;
t19 = t22 * t41;
t18 = -qJD(2) * pkin(2) - t35;
t11 = t14 ^ 2 / 0.2e1;
t3 = t12 * pkin(4) - t14 * qJ(5) + t9;
t2 = -qJ(5) * t37 + t5;
t1 = pkin(4) * t37 + qJD(5) - t4;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t29 / 0.2e1, 0, 0, 0, 0, 0, t41, t27 * t33, -t25 * t33, 0, (t25 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1) * t29, t21 * t41, t24 * t28 * t26, t24 * t32, t19, t26 * t32, qJD(3) ^ 2 / 0.2e1, -t18 * t37 - t24 * t36, t18 * t38 - t26 * t36, (t21 + t22) * t17 * qJD(2), t18 ^ 2 / 0.2e1 + (t22 / 0.2e1 + t21 / 0.2e1) * t17 ^ 2, t11, -t40, -t30, t34, t31, t19, t9 * t12 - t4 * t37, t9 * t14 + t5 * t37, -t5 * t12 - t4 * t14, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t11, -t30, t40, t19, -t31, t34, t1 * t37 + t3 * t12, t1 * t14 - t2 * t12, -t3 * t14 - t2 * t37, t2 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg = t6;
