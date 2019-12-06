% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPRP4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP4_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:36:05
% EndTime: 2019-12-05 15:36:05
% DurationCPUTime: 0.11s
% Computational Cost: add. (106->30), mult. (262->73), div. (0->0), fcn. (142->6), ass. (0->31)
t27 = qJD(2) ^ 2;
t20 = t27 / 0.2e1;
t23 = sin(qJ(4));
t25 = cos(qJ(4));
t26 = cos(qJ(2));
t12 = qJD(2) * pkin(2) + t26 * qJD(1);
t21 = sin(pkin(8));
t22 = cos(pkin(8));
t24 = sin(qJ(2));
t35 = qJD(1) * t24;
t10 = t21 * t12 + t22 * t35;
t8 = qJD(2) * pkin(6) + t10;
t5 = t23 * qJD(3) + t25 * t8;
t34 = qJD(2) * t23;
t33 = qJD(2) * t25;
t32 = qJD(1) * qJD(2);
t31 = qJD(2) * qJD(4);
t30 = t23 * t27 * t25;
t29 = t25 * t31;
t9 = t22 * t12 - t21 * t35;
t4 = t25 * qJD(3) - t23 * t8;
t28 = qJD(1) ^ 2;
t19 = qJD(4) ^ 2 / 0.2e1;
t17 = t23 * t31;
t16 = t25 ^ 2 * t20;
t15 = t23 ^ 2 * t20;
t7 = -qJD(2) * pkin(3) - t9;
t3 = (-pkin(4) * t25 - qJ(5) * t23 - pkin(3)) * qJD(2) - t9;
t2 = qJD(4) * qJ(5) + t5;
t1 = -qJD(4) * pkin(4) + qJD(5) - t4;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t28 / 0.2e1, 0, 0, 0, 0, 0, t20, t26 * t32, -t24 * t32, 0, (t24 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1) * t28, 0, 0, 0, 0, 0, t20, t9 * qJD(2), -t10 * qJD(2), 0, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t15, t30, t17, t16, t29, t19, t4 * qJD(4) - t7 * t33, -t5 * qJD(4) + t7 * t34, (-t23 * t4 + t25 * t5) * qJD(2), t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t15, t17, -t30, t19, -t29, t16, -t1 * qJD(4) - t3 * t33, (t1 * t23 + t2 * t25) * qJD(2), t2 * qJD(4) - t3 * t34, t2 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg = t6;
