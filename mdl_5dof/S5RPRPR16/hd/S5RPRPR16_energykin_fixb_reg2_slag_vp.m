% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRPR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR16_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR16_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:39:31
% EndTime: 2019-12-31 18:39:31
% DurationCPUTime: 0.14s
% Computational Cost: add. (148->42), mult. (314->87), div. (0->0), fcn. (119->4), ass. (0->37)
t29 = qJD(1) ^ 2;
t22 = t29 / 0.2e1;
t41 = cos(qJ(5));
t27 = sin(qJ(3));
t38 = qJD(1) * t27;
t40 = pkin(3) * t38 + qJD(1) * qJ(2);
t39 = t29 * qJ(2);
t13 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t37 = qJD(3) * t13;
t28 = cos(qJ(3));
t36 = t28 * qJD(1);
t35 = qJD(3) * qJ(4);
t34 = qJD(1) * qJD(3);
t33 = t28 * t29 * t27;
t32 = t27 * t34;
t31 = t28 * t34;
t30 = pkin(4) * qJD(1) - t13;
t26 = sin(qJ(5));
t25 = t28 ^ 2;
t24 = t27 ^ 2;
t21 = qJD(3) ^ 2 / 0.2e1;
t20 = qJ(2) ^ 2 * t22;
t19 = -qJD(1) * pkin(1) + qJD(2);
t17 = t25 * t22;
t16 = t24 * t22;
t15 = qJD(5) + t36;
t11 = qJD(3) * t41 + t26 * t38;
t9 = t26 * qJD(3) - t38 * t41;
t8 = -t27 * t13 - t35;
t7 = -qJ(4) * t36 + t40;
t6 = -qJD(3) * pkin(3) - t28 * t13 + qJD(4);
t5 = -t27 * t30 + t35;
t4 = (pkin(7) * t27 - qJ(4) * t28) * qJD(1) + t40;
t3 = qJD(4) + t30 * t28 + (-pkin(3) - pkin(7)) * qJD(3);
t2 = t26 * t3 + t4 * t41;
t1 = -t26 * t4 + t3 * t41;
t10 = [0, 0, 0, 0, 0, t22, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, t19 * qJD(1), t39, t20 + t19 ^ 2 / 0.2e1, t17, -t33, t31, t16, -t32, t21, t27 * t39 + t28 * t37, -t27 * t37 + t28 * t39, (-t24 - t25) * t13 * qJD(1), t20 + (t24 / 0.2e1 + t25 / 0.2e1) * t13 ^ 2, t21, -t31, t32, t17, -t33, t16, (t27 * t8 + t28 * t6) * qJD(1), t6 * qJD(3) - t38 * t7, -t8 * qJD(3) - t36 * t7, t7 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t11 ^ 2 / 0.2e1, -t11 * t9, t11 * t15, t9 ^ 2 / 0.2e1, -t9 * t15, t15 ^ 2 / 0.2e1, t1 * t15 + t5 * t9, t5 * t11 - t2 * t15, -t1 * t11 - t2 * t9, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t10;
