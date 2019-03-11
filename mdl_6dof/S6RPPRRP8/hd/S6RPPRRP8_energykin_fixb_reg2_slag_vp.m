% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRRP8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP8_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:16:23
% EndTime: 2019-03-09 02:16:23
% DurationCPUTime: 0.17s
% Computational Cost: add. (418->48), mult. (887->110), div. (0->0), fcn. (549->6), ass. (0->40)
t39 = sin(pkin(9));
t40 = cos(pkin(9));
t42 = sin(qJ(4));
t43 = cos(qJ(4));
t25 = (t39 * t43 + t40 * t42) * qJD(1);
t44 = qJD(1) ^ 2;
t37 = t44 / 0.2e1;
t27 = (-t39 * t42 + t40 * t43) * qJD(1);
t32 = qJD(1) * qJ(2) + qJD(3);
t48 = qJD(1) * t39;
t28 = pkin(3) * t48 + t32;
t10 = t25 * pkin(4) - t27 * pkin(8) + t28;
t41 = sin(qJ(5));
t51 = cos(qJ(5));
t30 = qJD(2) + (-pkin(1) - qJ(3)) * qJD(1);
t46 = -pkin(7) * qJD(1) + t30;
t22 = t46 * t39;
t23 = t46 * t40;
t12 = t43 * t22 + t42 * t23;
t9 = qJD(4) * pkin(8) + t12;
t4 = t41 * t10 + t51 * t9;
t15 = -t51 * qJD(4) + t41 * t27;
t17 = t41 * qJD(4) + t51 * t27;
t50 = t17 * t15;
t24 = qJD(5) + t25;
t49 = t24 * t15;
t47 = t15 ^ 2 / 0.2e1;
t11 = -t42 * t22 + t43 * t23;
t8 = -qJD(4) * pkin(4) - t11;
t3 = t51 * t10 - t41 * t9;
t36 = t40 ^ 2;
t35 = t39 ^ 2;
t33 = -pkin(1) * qJD(1) + qJD(2);
t21 = t24 ^ 2 / 0.2e1;
t14 = t17 ^ 2 / 0.2e1;
t13 = t17 * t24;
t5 = t15 * pkin(5) - t17 * qJ(6) + t8;
t2 = t24 * qJ(6) + t4;
t1 = -t24 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t37, 0, 0, 0, 0, t37, 0, 0, 0, 0, 0, 0, t33 * qJD(1), t44 * qJ(2), qJ(2) ^ 2 * t37 + t33 ^ 2 / 0.2e1, t36 * t37, -t40 * t44 * t39, 0, t35 * t37, 0, 0, t32 * t48, t32 * t40 * qJD(1) (-t35 - t36) * t30 * qJD(1), t32 ^ 2 / 0.2e1 + (t35 / 0.2e1 + t36 / 0.2e1) * t30 ^ 2, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * qJD(4), t25 ^ 2 / 0.2e1, -t25 * qJD(4), qJD(4) ^ 2 / 0.2e1, t11 * qJD(4) + t28 * t25, -t12 * qJD(4) + t28 * t27, -t11 * t27 - t12 * t25, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1, t14, -t50, t13, t47, -t49, t21, t8 * t15 + t3 * t24, t8 * t17 - t4 * t24, -t4 * t15 - t3 * t17, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1, t14, t13, t50, t21, t49, t47, -t1 * t24 + t5 * t15, t1 * t17 - t2 * t15, -t5 * t17 + t2 * t24, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
