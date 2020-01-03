% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP12_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP12_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:57:19
% EndTime: 2019-12-31 18:57:19
% DurationCPUTime: 0.12s
% Computational Cost: add. (179->39), mult. (382->82), div. (0->0), fcn. (180->4), ass. (0->33)
t33 = qJD(1) ^ 2;
t27 = t33 / 0.2e1;
t31 = sin(qJ(3));
t32 = cos(qJ(3));
t12 = (pkin(3) * t31 - pkin(7) * t32 + qJ(2)) * qJD(1);
t22 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t13 = qJD(3) * pkin(7) + t31 * t22;
t30 = sin(qJ(4));
t38 = cos(qJ(4));
t4 = t30 * t12 + t38 * t13;
t37 = t33 * qJ(2);
t36 = qJD(1) * t32;
t35 = qJD(3) * t22;
t34 = qJD(1) * qJD(3);
t3 = t38 * t12 - t30 * t13;
t14 = -qJD(3) * pkin(3) - t32 * t22;
t29 = t32 ^ 2;
t28 = t31 ^ 2;
t25 = qJ(2) ^ 2 * t27;
t24 = -pkin(1) * qJD(1) + qJD(2);
t23 = t31 * qJD(1) + qJD(4);
t20 = t23 ^ 2 / 0.2e1;
t19 = t30 * qJD(3) + t38 * t36;
t17 = -t38 * qJD(3) + t30 * t36;
t16 = t19 ^ 2 / 0.2e1;
t15 = t17 ^ 2 / 0.2e1;
t8 = t19 * t23;
t7 = t17 * t23;
t6 = t19 * t17;
t5 = t17 * pkin(4) + qJD(5) + t14;
t2 = -t17 * qJ(5) + t4;
t1 = t23 * pkin(4) - t19 * qJ(5) + t3;
t9 = [0, 0, 0, 0, 0, t27, 0, 0, 0, 0, t27, 0, 0, 0, 0, 0, 0, t24 * qJD(1), t37, t25 + t24 ^ 2 / 0.2e1, t29 * t27, -t32 * t33 * t31, t32 * t34, t28 * t27, -t31 * t34, qJD(3) ^ 2 / 0.2e1, t31 * t37 + t32 * t35, -t31 * t35 + t32 * t37, (-t28 - t29) * t22 * qJD(1), t25 + (t28 / 0.2e1 + t29 / 0.2e1) * t22 ^ 2, t16, -t6, t8, t15, -t7, t20, t14 * t17 + t3 * t23, t14 * t19 - t4 * t23, -t4 * t17 - t3 * t19, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t16, -t6, t8, t15, -t7, t20, t1 * t23 + t5 * t17, t5 * t19 - t2 * t23, -t1 * t19 - t2 * t17, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t9;
