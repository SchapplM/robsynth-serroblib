% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPPR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:15:37
% EndTime: 2019-03-09 08:15:37
% DurationCPUTime: 0.16s
% Computational Cost: add. (399->61), mult. (851->131), div. (0->0), fcn. (424->6), ass. (0->46)
t49 = qJD(1) ^ 2;
t60 = t49 / 0.2e1;
t59 = -pkin(2) - pkin(3);
t47 = sin(qJ(2));
t55 = t47 * qJD(1);
t48 = cos(qJ(2));
t56 = qJD(1) * t48;
t21 = -qJD(1) * pkin(1) - pkin(2) * t56 - qJ(3) * t55;
t16 = pkin(3) * t56 + qJD(4) - t21;
t13 = (pkin(4) * t47 + qJ(5) * t48) * qJD(1) + t16;
t53 = qJ(4) * qJD(1);
t54 = pkin(7) * t55 + qJD(3);
t51 = -t47 * t53 + t54;
t15 = (-qJ(5) + t59) * qJD(2) + t51;
t42 = sin(pkin(9));
t43 = cos(pkin(9));
t6 = t42 * t13 + t43 * t15;
t58 = cos(qJ(6));
t57 = t48 * t49;
t26 = pkin(7) * t56 + qJD(2) * qJ(3);
t52 = qJD(1) * qJD(2);
t5 = t43 * t13 - t42 * t15;
t31 = t47 * t52;
t32 = t48 * t52;
t20 = t48 * t53 - t26;
t17 = qJD(2) * pkin(4) + qJD(5) - t20;
t46 = sin(qJ(6));
t41 = t48 ^ 2;
t40 = t47 ^ 2;
t38 = qJD(2) ^ 2 / 0.2e1;
t30 = t41 * t60;
t29 = t40 * t60;
t28 = qJD(6) + t55;
t27 = t47 * t57;
t25 = -qJD(2) * pkin(2) + t54;
t23 = t42 * qJD(2) + t43 * t56;
t22 = -t43 * qJD(2) + t42 * t56;
t18 = t59 * qJD(2) + t51;
t12 = t46 * t22 - t58 * t23;
t10 = -t58 * t22 - t46 * t23;
t9 = -t22 * pkin(5) + t17;
t4 = t22 * pkin(8) + t6;
t3 = pkin(5) * t55 + t23 * pkin(8) + t5;
t2 = t46 * t3 + t58 * t4;
t1 = t58 * t3 - t46 * t4;
t7 = [0, 0, 0, 0, 0, t60, 0, 0, 0, 0, t29, t27, t31, t30, t32, t38, pkin(1) * t57 - pkin(7) * t31, -t49 * pkin(1) * t47 - pkin(7) * t32 (t40 + t41) * t49 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t41 / 0.2e1 + t40 / 0.2e1) * pkin(7) ^ 2) * t49, t29, t31, -t27, t38, -t32, t30, -t25 * qJD(2) - t21 * t56 (t25 * t47 + t26 * t48) * qJD(1), t26 * qJD(2) - t21 * t55, t26 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1, t30, t27, t32, t29, t31, t38, -t20 * qJD(2) + t16 * t55, t18 * qJD(2) - t16 * t56 (-t18 * t47 + t20 * t48) * qJD(1), t18 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1, t23 ^ 2 / 0.2e1, -t23 * t22, -t23 * t55, t22 ^ 2 / 0.2e1, t22 * t55, t29, -t17 * t22 + t5 * t55, -t17 * t23 - t6 * t55, t6 * t22 + t5 * t23, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t12 ^ 2 / 0.2e1, -t12 * t10, t12 * t28, t10 ^ 2 / 0.2e1, -t10 * t28, t28 ^ 2 / 0.2e1, t1 * t28 + t9 * t10, t9 * t12 - t2 * t28, -t1 * t12 - t2 * t10, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1;];
T_reg  = t7;
