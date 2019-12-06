% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRPRP1
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
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPRP1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP1_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:38
% EndTime: 2019-12-05 15:28:38
% DurationCPUTime: 0.12s
% Computational Cost: add. (134->33), mult. (366->75), div. (0->0), fcn. (222->4), ass. (0->31)
t28 = qJD(2) ^ 2;
t36 = t28 / 0.2e1;
t26 = cos(pkin(8));
t22 = t26 * qJD(1);
t25 = sin(pkin(8));
t33 = qJD(2) * t25;
t10 = t22 + (-pkin(6) - qJ(3)) * t33;
t29 = qJ(3) * qJD(2);
t16 = t25 * qJD(1) + t26 * t29;
t32 = qJD(2) * t26;
t11 = pkin(6) * t32 + t16;
t27 = sin(qJ(4));
t35 = cos(qJ(4));
t5 = t27 * t10 + t35 * t11;
t12 = t27 * t33 - t35 * t32;
t14 = (t35 * t25 + t26 * t27) * qJD(2);
t34 = t14 * t12;
t31 = qJD(4) * t12;
t30 = t12 ^ 2 / 0.2e1;
t4 = t35 * t10 - t27 * t11;
t17 = qJD(3) + (-pkin(3) * t26 - pkin(2)) * qJD(2);
t24 = qJD(1) ^ 2 / 0.2e1;
t23 = qJD(4) ^ 2 / 0.2e1;
t20 = -qJD(2) * pkin(2) + qJD(3);
t15 = -t25 * t29 + t22;
t9 = t14 * qJD(4);
t8 = t14 ^ 2 / 0.2e1;
t3 = qJD(4) * qJ(5) + t5;
t2 = t12 * pkin(4) - t14 * qJ(5) + t17;
t1 = -qJD(4) * pkin(4) + qJD(5) - t4;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, 0, 0, 0, 0, t36, 0, 0, 0, t24, t25 ^ 2 * t36, t25 * t28 * t26, 0, t26 ^ 2 * t36, 0, 0, -t20 * t32, t20 * t33, (-t15 * t25 + t16 * t26) * qJD(2), t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1, t8, -t34, t9, t30, -t31, t23, t4 * qJD(4) + t17 * t12, -t5 * qJD(4) + t17 * t14, -t5 * t12 - t4 * t14, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t8, t9, t34, t23, t31, t30, -t1 * qJD(4) + t2 * t12, t1 * t14 - t3 * t12, t3 * qJD(4) - t2 * t14, t3 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg = t6;
