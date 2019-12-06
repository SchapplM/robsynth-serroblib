% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRPRP5
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
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPRP5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP5_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:38:40
% EndTime: 2019-12-05 15:38:40
% DurationCPUTime: 0.12s
% Computational Cost: add. (145->33), mult. (387->82), div. (0->0), fcn. (233->6), ass. (0->35)
t28 = qJD(2) ^ 2;
t39 = t28 / 0.2e1;
t25 = sin(qJ(4));
t38 = cos(qJ(4));
t23 = sin(pkin(8));
t26 = sin(qJ(2));
t18 = qJD(2) * qJ(3) + t26 * qJD(1);
t31 = pkin(6) * qJD(2) + t18;
t8 = t31 * t23;
t24 = cos(pkin(8));
t9 = t31 * t24;
t5 = -t25 * t8 + t38 * t9;
t35 = qJD(2) * t24;
t36 = qJD(2) * t23;
t12 = t25 * t36 - t38 * t35;
t14 = (t38 * t23 + t24 * t25) * qJD(2);
t37 = t14 * t12;
t34 = qJD(4) * t12;
t33 = t12 ^ 2 / 0.2e1;
t32 = qJD(1) * qJD(2);
t27 = cos(qJ(2));
t30 = -t27 * qJD(1) + qJD(3);
t4 = -t25 * t9 - t38 * t8;
t15 = (-pkin(3) * t24 - pkin(2)) * qJD(2) + t30;
t29 = qJD(1) ^ 2;
t22 = qJD(4) ^ 2 / 0.2e1;
t21 = t24 ^ 2;
t20 = t23 ^ 2;
t16 = -qJD(2) * pkin(2) + t30;
t11 = t14 * qJD(4);
t10 = t14 ^ 2 / 0.2e1;
t3 = qJD(4) * qJ(5) + t5;
t2 = -qJD(4) * pkin(4) + qJD(5) - t4;
t1 = t12 * pkin(4) - t14 * qJ(5) + t15;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t29 / 0.2e1, 0, 0, 0, 0, 0, t39, t27 * t32, -t26 * t32, 0, (t26 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1) * t29, t20 * t39, t23 * t28 * t24, 0, t21 * t39, 0, 0, -t16 * t35, t16 * t36, (t20 + t21) * t18 * qJD(2), t16 ^ 2 / 0.2e1 + (t21 / 0.2e1 + t20 / 0.2e1) * t18 ^ 2, t10, -t37, t11, t33, -t34, t22, t4 * qJD(4) + t15 * t12, -t5 * qJD(4) + t15 * t14, -t5 * t12 - t4 * t14, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t10, t11, t37, t22, t34, t33, -t2 * qJD(4) + t1 * t12, -t3 * t12 + t2 * t14, t3 * qJD(4) - t1 * t14, t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1;];
T_reg = t6;
