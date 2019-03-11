% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRP8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP8_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:24:32
% EndTime: 2019-03-09 06:24:32
% DurationCPUTime: 0.16s
% Computational Cost: add. (480->52), mult. (926->114), div. (0->0), fcn. (553->6), ass. (0->43)
t41 = sin(qJ(4));
t42 = sin(qJ(3));
t43 = cos(qJ(4));
t44 = cos(qJ(3));
t25 = (t41 * t44 + t42 * t43) * qJD(1);
t45 = qJD(1) ^ 2;
t36 = t45 / 0.2e1;
t27 = (-t41 * t42 + t43 * t44) * qJD(1);
t28 = (pkin(3) * t42 + qJ(2)) * qJD(1);
t13 = t25 * pkin(4) - t27 * pkin(9) + t28;
t40 = sin(qJ(5));
t54 = cos(qJ(5));
t30 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t47 = -pkin(8) * qJD(1) + t30;
t21 = qJD(3) * pkin(3) + t47 * t44;
t23 = t47 * t42;
t11 = t41 * t21 + t43 * t23;
t35 = qJD(3) + qJD(4);
t8 = t35 * pkin(9) + t11;
t5 = t40 * t13 + t54 * t8;
t15 = t40 * t27 - t54 * t35;
t17 = t54 * t27 + t40 * t35;
t53 = t17 * t15;
t24 = qJD(5) + t25;
t52 = t24 * t15;
t51 = t45 * qJ(2);
t50 = qJD(3) * t30;
t49 = t15 ^ 2 / 0.2e1;
t48 = qJD(1) * qJD(3);
t10 = t43 * t21 - t41 * t23;
t7 = -t35 * pkin(4) - t10;
t4 = t54 * t13 - t40 * t8;
t39 = t44 ^ 2;
t38 = t42 ^ 2;
t34 = qJ(2) ^ 2 * t36;
t33 = -pkin(1) * qJD(1) + qJD(2);
t22 = t24 ^ 2 / 0.2e1;
t14 = t17 ^ 2 / 0.2e1;
t12 = t17 * t24;
t3 = t15 * pkin(5) - t17 * qJ(6) + t7;
t2 = t24 * qJ(6) + t5;
t1 = -t24 * pkin(5) + qJD(6) - t4;
t6 = [0, 0, 0, 0, 0, t36, 0, 0, 0, 0, t36, 0, 0, 0, 0, 0, 0, t33 * qJD(1), t51, t34 + t33 ^ 2 / 0.2e1, t39 * t36, -t44 * t45 * t42, t44 * t48, t38 * t36, -t42 * t48, qJD(3) ^ 2 / 0.2e1, t42 * t51 + t44 * t50, -t42 * t50 + t44 * t51 (-t38 - t39) * t30 * qJD(1), t34 + (t38 / 0.2e1 + t39 / 0.2e1) * t30 ^ 2, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t35, t25 ^ 2 / 0.2e1, -t25 * t35, t35 ^ 2 / 0.2e1, t10 * t35 + t28 * t25, -t11 * t35 + t28 * t27, -t10 * t27 - t11 * t25, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1, t14, -t53, t12, t49, -t52, t22, t7 * t15 + t4 * t24, t7 * t17 - t5 * t24, -t5 * t15 - t4 * t17, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t14, t12, t53, t22, t52, t49, -t1 * t24 + t3 * t15, t1 * t17 - t2 * t15, -t3 * t17 + t2 * t24, t2 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
