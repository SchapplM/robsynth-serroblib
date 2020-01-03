% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRP6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP6_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:55:15
% EndTime: 2019-12-31 17:55:16
% DurationCPUTime: 0.11s
% Computational Cost: add. (165->32), mult. (377->76), div. (0->0), fcn. (192->4), ass. (0->31)
t31 = qJD(1) ^ 2;
t25 = t31 / 0.2e1;
t29 = sin(qJ(4));
t30 = cos(qJ(4));
t27 = sin(pkin(7));
t17 = qJD(2) + (-pkin(1) - qJ(3)) * qJD(1);
t32 = -pkin(6) * qJD(1) + t17;
t8 = t32 * t27;
t28 = cos(pkin(7));
t9 = t32 * t28;
t5 = t29 * t9 + t30 * t8;
t12 = (-t27 * t30 - t28 * t29) * qJD(1);
t35 = qJD(1) * t28;
t36 = qJD(1) * t27;
t14 = -t29 * t36 + t30 * t35;
t37 = t14 * t12;
t34 = qJD(4) * t12;
t33 = t12 ^ 2 / 0.2e1;
t20 = qJD(1) * qJ(2) + qJD(3);
t15 = pkin(3) * t36 + t20;
t4 = -t29 * t8 + t30 * t9;
t24 = qJD(4) ^ 2 / 0.2e1;
t23 = t28 ^ 2;
t22 = t27 ^ 2;
t21 = -qJD(1) * pkin(1) + qJD(2);
t11 = t14 * qJD(4);
t10 = t14 ^ 2 / 0.2e1;
t3 = -t12 * pkin(4) - t14 * qJ(5) + t15;
t2 = qJD(4) * qJ(5) + t5;
t1 = -qJD(4) * pkin(4) + qJD(5) - t4;
t6 = [0, 0, 0, 0, 0, t25, 0, 0, 0, 0, t25, 0, 0, 0, 0, 0, 0, t21 * qJD(1), t31 * qJ(2), qJ(2) ^ 2 * t25 + t21 ^ 2 / 0.2e1, t23 * t25, -t28 * t31 * t27, 0, t22 * t25, 0, 0, t20 * t36, t20 * t35, (-t22 - t23) * t17 * qJD(1), t20 ^ 2 / 0.2e1 + (t22 / 0.2e1 + t23 / 0.2e1) * t17 ^ 2, t10, t37, t11, t33, t34, t24, t4 * qJD(4) - t15 * t12, -t5 * qJD(4) + t15 * t14, t5 * t12 - t4 * t14, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t10, t11, -t37, t24, -t34, t33, -t1 * qJD(4) - t3 * t12, t1 * t14 + t2 * t12, t2 * qJD(4) - t3 * t14, t2 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg = t6;
