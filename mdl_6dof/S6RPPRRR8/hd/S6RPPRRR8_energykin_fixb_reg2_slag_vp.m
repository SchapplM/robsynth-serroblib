% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRRR8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR8_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:36:19
% EndTime: 2019-03-09 02:36:19
% DurationCPUTime: 0.17s
% Computational Cost: add. (573->54), mult. (1210->130), div. (0->0), fcn. (797->8), ass. (0->41)
t44 = sin(pkin(10));
t45 = cos(pkin(10));
t48 = sin(qJ(4));
t49 = cos(qJ(4));
t30 = (t44 * t49 + t45 * t48) * qJD(1);
t50 = qJD(1) ^ 2;
t42 = t50 / 0.2e1;
t55 = cos(qJ(5));
t54 = cos(qJ(6));
t35 = qJD(2) + (-pkin(1) - qJ(3)) * qJD(1);
t52 = -pkin(7) * qJD(1) + t35;
t27 = t52 * t44;
t28 = t52 * t45;
t18 = t49 * t27 + t48 * t28;
t15 = qJD(4) * pkin(8) + t18;
t32 = (-t44 * t48 + t45 * t49) * qJD(1);
t37 = qJD(1) * qJ(2) + qJD(3);
t53 = qJD(1) * t44;
t33 = pkin(3) * t53 + t37;
t16 = t30 * pkin(4) - t32 * pkin(8) + t33;
t47 = sin(qJ(5));
t6 = t55 * t15 + t47 * t16;
t5 = -t47 * t15 + t55 * t16;
t17 = -t48 * t27 + t49 * t28;
t14 = -qJD(4) * pkin(4) - t17;
t29 = qJD(5) + t30;
t46 = sin(qJ(6));
t41 = t45 ^ 2;
t40 = t44 ^ 2;
t38 = -pkin(1) * qJD(1) + qJD(2);
t26 = qJD(6) + t29;
t22 = t47 * qJD(4) + t55 * t32;
t20 = -t55 * qJD(4) + t47 * t32;
t10 = -t46 * t20 + t54 * t22;
t8 = t54 * t20 + t46 * t22;
t7 = t20 * pkin(5) + t14;
t4 = -t20 * pkin(9) + t6;
t3 = t29 * pkin(5) - t22 * pkin(9) + t5;
t2 = t46 * t3 + t54 * t4;
t1 = t54 * t3 - t46 * t4;
t9 = [0, 0, 0, 0, 0, t42, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, t38 * qJD(1), t50 * qJ(2), qJ(2) ^ 2 * t42 + t38 ^ 2 / 0.2e1, t41 * t42, -t45 * t50 * t44, 0, t40 * t42, 0, 0, t37 * t53, t37 * t45 * qJD(1) (-t40 - t41) * t35 * qJD(1), t37 ^ 2 / 0.2e1 + (t40 / 0.2e1 + t41 / 0.2e1) * t35 ^ 2, t32 ^ 2 / 0.2e1, -t32 * t30, t32 * qJD(4), t30 ^ 2 / 0.2e1, -t30 * qJD(4), qJD(4) ^ 2 / 0.2e1, t17 * qJD(4) + t33 * t30, -t18 * qJD(4) + t33 * t32, -t17 * t32 - t18 * t30, t18 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * t29, t20 ^ 2 / 0.2e1, -t20 * t29, t29 ^ 2 / 0.2e1, t14 * t20 + t5 * t29, t14 * t22 - t6 * t29, -t6 * t20 - t5 * t22, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, t10 * t26, t8 ^ 2 / 0.2e1, -t8 * t26, t26 ^ 2 / 0.2e1, t1 * t26 + t7 * t8, t7 * t10 - t2 * t26, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t9;
