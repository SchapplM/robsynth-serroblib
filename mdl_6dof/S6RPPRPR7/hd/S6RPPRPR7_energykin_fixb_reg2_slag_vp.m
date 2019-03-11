% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRPR7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR7_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:53:43
% EndTime: 2019-03-09 01:53:43
% DurationCPUTime: 0.17s
% Computational Cost: add. (563->54), mult. (1210->127), div. (0->0), fcn. (797->8), ass. (0->41)
t44 = sin(pkin(9));
t45 = cos(pkin(9));
t47 = sin(qJ(4));
t48 = cos(qJ(4));
t29 = (t44 * t48 + t45 * t47) * qJD(1);
t49 = qJD(1) ^ 2;
t41 = t49 / 0.2e1;
t55 = cos(qJ(6));
t34 = qJD(2) + (-pkin(1) - qJ(3)) * qJD(1);
t51 = -pkin(7) * qJD(1) + t34;
t26 = t51 * t44;
t27 = t51 * t45;
t18 = t48 * t26 + t47 * t27;
t15 = qJD(4) * qJ(5) + t18;
t31 = (-t44 * t47 + t45 * t48) * qJD(1);
t36 = qJD(1) * qJ(2) + qJD(3);
t53 = qJD(1) * t44;
t32 = pkin(3) * t53 + t36;
t16 = t29 * pkin(4) - t31 * qJ(5) + t32;
t43 = sin(pkin(10));
t54 = cos(pkin(10));
t6 = t54 * t15 + t43 * t16;
t52 = t29 ^ 2 / 0.2e1;
t5 = -t43 * t15 + t54 * t16;
t17 = -t47 * t26 + t48 * t27;
t14 = -qJD(4) * pkin(4) + qJD(5) - t17;
t46 = sin(qJ(6));
t40 = t45 ^ 2;
t39 = t44 ^ 2;
t37 = -pkin(1) * qJD(1) + qJD(2);
t28 = qJD(6) + t29;
t22 = t43 * qJD(4) + t54 * t31;
t20 = -t54 * qJD(4) + t43 * t31;
t10 = -t46 * t20 + t55 * t22;
t8 = t55 * t20 + t46 * t22;
t7 = t20 * pkin(5) + t14;
t4 = -t20 * pkin(8) + t6;
t3 = t29 * pkin(5) - t22 * pkin(8) + t5;
t2 = t46 * t3 + t55 * t4;
t1 = t55 * t3 - t46 * t4;
t9 = [0, 0, 0, 0, 0, t41, 0, 0, 0, 0, t41, 0, 0, 0, 0, 0, 0, t37 * qJD(1), t49 * qJ(2), qJ(2) ^ 2 * t41 + t37 ^ 2 / 0.2e1, t40 * t41, -t45 * t49 * t44, 0, t39 * t41, 0, 0, t36 * t53, t36 * t45 * qJD(1) (-t39 - t40) * t34 * qJD(1), t36 ^ 2 / 0.2e1 + (t39 / 0.2e1 + t40 / 0.2e1) * t34 ^ 2, t31 ^ 2 / 0.2e1, -t31 * t29, t31 * qJD(4), t52, -t29 * qJD(4), qJD(4) ^ 2 / 0.2e1, t17 * qJD(4) + t32 * t29, -t18 * qJD(4) + t32 * t31, -t17 * t31 - t18 * t29, t18 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * t29, t20 ^ 2 / 0.2e1, -t20 * t29, t52, t14 * t20 + t5 * t29, t14 * t22 - t6 * t29, -t6 * t20 - t5 * t22, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, t10 * t28, t8 ^ 2 / 0.2e1, -t8 * t28, t28 ^ 2 / 0.2e1, t1 * t28 + t7 * t8, t7 * t10 - t2 * t28, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t9;
