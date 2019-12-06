% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR2_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:55
% EndTime: 2019-12-05 17:39:55
% DurationCPUTime: 0.14s
% Computational Cost: add. (253->36), mult. (566->97), div. (0->0), fcn. (336->6), ass. (0->31)
t36 = qJD(1) ^ 2;
t29 = t36 / 0.2e1;
t39 = cos(qJ(5));
t31 = sin(pkin(8));
t22 = qJD(2) + (-pkin(1) - qJ(3)) * qJD(1);
t37 = -pkin(6) * qJD(1) + t22;
t15 = t37 * t31;
t32 = cos(pkin(8));
t16 = t37 * t32;
t34 = sin(qJ(4));
t35 = cos(qJ(4));
t6 = t35 * t15 + t34 * t16;
t38 = qJD(1) * t31;
t24 = qJD(1) * qJ(2) + qJD(3);
t20 = pkin(3) * t38 + t24;
t5 = -t34 * t15 + t35 * t16;
t33 = sin(qJ(5));
t28 = qJD(4) + qJD(5);
t27 = t32 ^ 2;
t26 = t31 ^ 2;
t25 = -pkin(1) * qJD(1) + qJD(2);
t19 = (-t31 * t34 + t32 * t35) * qJD(1);
t17 = (-t31 * t35 - t32 * t34) * qJD(1);
t10 = -t17 * pkin(4) + t20;
t9 = t33 * t17 + t39 * t19;
t7 = -t39 * t17 + t33 * t19;
t4 = t17 * pkin(7) + t6;
t3 = qJD(4) * pkin(4) - t19 * pkin(7) + t5;
t2 = t33 * t3 + t39 * t4;
t1 = t39 * t3 - t33 * t4;
t8 = [0, 0, 0, 0, 0, t29, 0, 0, 0, 0, t29, 0, 0, 0, 0, 0, 0, t25 * qJD(1), t36 * qJ(2), qJ(2) ^ 2 * t29 + t25 ^ 2 / 0.2e1, t27 * t29, -t32 * t36 * t31, 0, t26 * t29, 0, 0, t24 * t38, t24 * t32 * qJD(1), (-t26 - t27) * t22 * qJD(1), t24 ^ 2 / 0.2e1 + (t26 / 0.2e1 + t27 / 0.2e1) * t22 ^ 2, t19 ^ 2 / 0.2e1, t19 * t17, t19 * qJD(4), t17 ^ 2 / 0.2e1, t17 * qJD(4), qJD(4) ^ 2 / 0.2e1, t5 * qJD(4) - t20 * t17, -t6 * qJD(4) + t20 * t19, t6 * t17 - t5 * t19, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1, t9 ^ 2 / 0.2e1, -t9 * t7, t9 * t28, t7 ^ 2 / 0.2e1, -t7 * t28, t28 ^ 2 / 0.2e1, t1 * t28 + t10 * t7, t10 * t9 - t2 * t28, -t1 * t9 - t2 * t7, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg = t8;
