% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPRP3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:35:32
% EndTime: 2019-03-09 08:35:32
% DurationCPUTime: 0.14s
% Computational Cost: add. (291->57), mult. (638->113), div. (0->0), fcn. (287->4), ass. (0->46)
t49 = qJD(1) ^ 2;
t59 = t49 / 0.2e1;
t58 = -pkin(2) - pkin(3);
t46 = sin(qJ(2));
t53 = qJ(4) * qJD(1);
t55 = t46 * qJD(1);
t54 = pkin(7) * t55 + qJD(3);
t51 = -t46 * t53 + t54;
t11 = (-pkin(8) + t58) * qJD(2) + t51;
t45 = sin(qJ(5));
t47 = cos(qJ(5));
t48 = cos(qJ(2));
t56 = qJD(1) * t48;
t18 = -qJD(1) * pkin(1) - pkin(2) * t56 - qJ(3) * t55;
t12 = pkin(3) * t56 + qJD(4) - t18;
t8 = (pkin(4) * t46 + pkin(8) * t48) * qJD(1) + t12;
t4 = t47 * t11 + t45 * t8;
t57 = t48 * t49;
t25 = pkin(7) * t56 + qJD(2) * qJ(3);
t52 = qJD(1) * qJD(2);
t3 = -t45 * t11 + t47 * t8;
t31 = t46 * t52;
t32 = t48 * t52;
t17 = t48 * t53 - t25;
t16 = qJD(2) * pkin(4) - t17;
t42 = t48 ^ 2;
t41 = t46 ^ 2;
t39 = qJD(2) ^ 2 / 0.2e1;
t30 = t42 * t59;
t29 = t41 * t59;
t28 = qJD(5) + t55;
t27 = t46 * t57;
t26 = t28 ^ 2 / 0.2e1;
t24 = -qJD(2) * pkin(2) + t54;
t22 = t45 * qJD(2) + t47 * t56;
t21 = -t47 * qJD(2) + t45 * t56;
t20 = t22 ^ 2 / 0.2e1;
t19 = t21 ^ 2 / 0.2e1;
t15 = t58 * qJD(2) + t51;
t14 = t22 * t28;
t13 = t21 * t28;
t9 = t22 * t21;
t7 = -t21 * pkin(5) + qJD(6) + t16;
t2 = t21 * qJ(6) + t4;
t1 = t28 * pkin(5) + t22 * qJ(6) + t3;
t5 = [0, 0, 0, 0, 0, t59, 0, 0, 0, 0, t29, t27, t31, t30, t32, t39, pkin(1) * t57 - pkin(7) * t31, -t49 * pkin(1) * t46 - pkin(7) * t32 (t41 + t42) * t49 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t42 / 0.2e1 + t41 / 0.2e1) * pkin(7) ^ 2) * t49, t29, t31, -t27, t39, -t32, t30, -t24 * qJD(2) - t18 * t56 (t24 * t46 + t25 * t48) * qJD(1), t25 * qJD(2) - t18 * t55, t25 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1, t30, t27, t32, t29, t31, t39, -t17 * qJD(2) + t12 * t55, t15 * qJD(2) - t12 * t56 (-t15 * t46 + t17 * t48) * qJD(1), t15 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1, t20, -t9, -t14, t19, t13, t26, -t16 * t21 + t3 * t28, -t16 * t22 - t4 * t28, t4 * t21 + t3 * t22, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1, t20, -t9, -t14, t19, t13, t26, t1 * t28 - t7 * t21, -t2 * t28 - t7 * t22, t1 * t22 + t2 * t21, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t5;
