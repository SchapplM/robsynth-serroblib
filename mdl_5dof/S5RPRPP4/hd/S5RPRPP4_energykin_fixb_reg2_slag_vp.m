% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPP4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:49
% EndTime: 2019-12-31 18:14:49
% DurationCPUTime: 0.14s
% Computational Cost: add. (181->36), mult. (404->80), div. (0->0), fcn. (196->4), ass. (0->33)
t31 = qJD(1) ^ 2;
t23 = t31 / 0.2e1;
t29 = sin(qJ(3));
t17 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t32 = -qJ(4) * qJD(1) + t17;
t11 = t32 * t29;
t27 = sin(pkin(7));
t28 = cos(pkin(7));
t30 = cos(qJ(3));
t8 = qJD(3) * pkin(3) + t32 * t30;
t4 = t28 * t11 + t27 * t8;
t12 = (-t27 * t30 - t28 * t29) * qJD(1);
t37 = qJD(1) * t29;
t14 = t28 * t30 * qJD(1) - t27 * t37;
t39 = t14 * t12;
t38 = t31 * qJ(2);
t36 = qJD(3) * t12;
t35 = qJD(3) * t17;
t34 = t12 ^ 2 / 0.2e1;
t33 = qJD(1) * qJD(3);
t15 = pkin(3) * t37 + qJD(1) * qJ(2) + qJD(4);
t3 = -t27 * t11 + t28 * t8;
t26 = t30 ^ 2;
t25 = t29 ^ 2;
t22 = qJD(3) ^ 2 / 0.2e1;
t21 = qJ(2) ^ 2 * t23;
t20 = -pkin(1) * qJD(1) + qJD(2);
t10 = t14 * qJD(3);
t9 = t14 ^ 2 / 0.2e1;
t5 = -t12 * pkin(4) - t14 * qJ(5) + t15;
t2 = qJD(3) * qJ(5) + t4;
t1 = -qJD(3) * pkin(4) + qJD(5) - t3;
t6 = [0, 0, 0, 0, 0, t23, 0, 0, 0, 0, t23, 0, 0, 0, 0, 0, 0, t20 * qJD(1), t38, t21 + t20 ^ 2 / 0.2e1, t26 * t23, -t30 * t31 * t29, t30 * t33, t25 * t23, -t29 * t33, t22, t29 * t38 + t30 * t35, -t29 * t35 + t30 * t38, (-t25 - t26) * t17 * qJD(1), t21 + (t25 / 0.2e1 + t26 / 0.2e1) * t17 ^ 2, t9, t39, t10, t34, t36, t22, t3 * qJD(3) - t15 * t12, -t4 * qJD(3) + t15 * t14, t4 * t12 - t3 * t14, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t9, t10, -t39, t22, -t36, t34, -t1 * qJD(3) - t5 * t12, t1 * t14 + t2 * t12, t2 * qJD(3) - t5 * t14, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg = t6;
