% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRP7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP7_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:22:45
% EndTime: 2019-03-09 03:22:45
% DurationCPUTime: 0.16s
% Computational Cost: add. (446->54), mult. (932->111), div. (0->0), fcn. (559->6), ass. (0->43)
t44 = sin(pkin(9));
t45 = cos(pkin(9));
t47 = sin(qJ(3));
t48 = cos(qJ(3));
t29 = (t44 * t48 + t45 * t47) * qJD(1);
t49 = qJD(1) ^ 2;
t40 = t49 / 0.2e1;
t31 = (-t44 * t47 + t45 * t48) * qJD(1);
t32 = qJD(4) + (pkin(3) * t47 + qJ(2)) * qJD(1);
t14 = t29 * pkin(4) - t31 * pkin(8) + t32;
t46 = sin(qJ(5));
t55 = cos(qJ(5));
t34 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t51 = -qJ(4) * qJD(1) + t34;
t25 = qJD(3) * pkin(3) + t51 * t48;
t27 = t51 * t47;
t13 = t44 * t25 + t45 * t27;
t9 = qJD(3) * pkin(8) + t13;
t4 = t46 * t14 + t55 * t9;
t54 = t49 * qJ(2);
t53 = qJD(3) * t34;
t52 = qJD(1) * qJD(3);
t3 = t55 * t14 - t46 * t9;
t12 = t45 * t25 - t44 * t27;
t8 = -qJD(3) * pkin(4) - t12;
t43 = t48 ^ 2;
t42 = t47 ^ 2;
t39 = qJD(3) ^ 2 / 0.2e1;
t37 = qJ(2) ^ 2 * t40;
t36 = -pkin(1) * qJD(1) + qJD(2);
t28 = qJD(5) + t29;
t26 = t28 ^ 2 / 0.2e1;
t21 = t46 * qJD(3) + t55 * t31;
t19 = -t55 * qJD(3) + t46 * t31;
t18 = t21 ^ 2 / 0.2e1;
t17 = t19 ^ 2 / 0.2e1;
t16 = t21 * t28;
t15 = t19 * t28;
t7 = t21 * t19;
t5 = t19 * pkin(5) + qJD(6) + t8;
t2 = -t19 * qJ(6) + t4;
t1 = t28 * pkin(5) - t21 * qJ(6) + t3;
t6 = [0, 0, 0, 0, 0, t40, 0, 0, 0, 0, t40, 0, 0, 0, 0, 0, 0, t36 * qJD(1), t54, t37 + t36 ^ 2 / 0.2e1, t43 * t40, -t48 * t49 * t47, t48 * t52, t42 * t40, -t47 * t52, t39, t47 * t54 + t48 * t53, -t47 * t53 + t48 * t54 (-t42 - t43) * t34 * qJD(1), t37 + (t42 / 0.2e1 + t43 / 0.2e1) * t34 ^ 2, t31 ^ 2 / 0.2e1, -t31 * t29, t31 * qJD(3), t29 ^ 2 / 0.2e1, -t29 * qJD(3), t39, t12 * qJD(3) + t32 * t29, -t13 * qJD(3) + t32 * t31, -t12 * t31 - t13 * t29, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1, t18, -t7, t16, t17, -t15, t26, t8 * t19 + t3 * t28, t8 * t21 - t4 * t28, -t4 * t19 - t3 * t21, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1, t18, -t7, t16, t17, -t15, t26, t1 * t28 + t5 * t19, -t2 * t28 + t5 * t21, -t1 * t21 - t2 * t19, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t6;
