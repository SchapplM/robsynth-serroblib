% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRRP2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP2_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:04:18
% EndTime: 2019-03-09 00:04:18
% DurationCPUTime: 0.15s
% Computational Cost: add. (496->56), mult. (1132->134), div. (0->0), fcn. (816->10), ass. (0->51)
t48 = qJD(2) ^ 2;
t62 = t48 / 0.2e1;
t46 = cos(qJ(3));
t47 = cos(qJ(2));
t40 = sin(pkin(6));
t57 = qJD(1) * t40;
t51 = t47 * t57;
t26 = -t51 + (-pkin(3) * t46 - pkin(2)) * qJD(2);
t43 = sin(qJ(4));
t54 = qJD(2) * t46;
t44 = sin(qJ(3));
t55 = qJD(2) * t44;
t61 = cos(qJ(4));
t28 = t43 * t55 - t61 * t54;
t30 = (t43 * t46 + t61 * t44) * qJD(2);
t12 = t28 * pkin(4) - t30 * pkin(10) + t26;
t42 = sin(qJ(5));
t60 = cos(qJ(5));
t45 = sin(qJ(2));
t32 = qJD(2) * pkin(8) + t45 * t57;
t41 = cos(pkin(6));
t56 = qJD(1) * t41;
t35 = t46 * t56;
t17 = qJD(3) * pkin(3) + t35 + (-pkin(9) * qJD(2) - t32) * t44;
t24 = t46 * t32 + t44 * t56;
t18 = pkin(9) * t54 + t24;
t10 = t43 * t17 + t61 * t18;
t39 = qJD(3) + qJD(4);
t8 = t39 * pkin(10) + t10;
t4 = t42 * t12 + t60 * t8;
t20 = t42 * t30 - t60 * t39;
t22 = t60 * t30 + t42 * t39;
t59 = t22 * t20;
t27 = qJD(5) + t28;
t58 = t27 * t20;
t53 = t20 ^ 2 / 0.2e1;
t52 = qJD(2) * qJD(3);
t50 = qJD(2) * t57;
t9 = t61 * t17 - t43 * t18;
t3 = t60 * t12 - t42 * t8;
t7 = -t39 * pkin(4) - t9;
t49 = qJD(1) ^ 2;
t33 = -qJD(2) * pkin(2) - t51;
t25 = t27 ^ 2 / 0.2e1;
t23 = -t44 * t32 + t35;
t19 = t22 ^ 2 / 0.2e1;
t13 = t22 * t27;
t5 = t20 * pkin(5) - t22 * qJ(6) + t7;
t2 = t27 * qJ(6) + t4;
t1 = -t27 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t49 / 0.2e1, 0, 0, 0, 0, 0, t62, t47 * t50, -t45 * t50, 0 (t41 ^ 2 / 0.2e1 + (t45 ^ 2 / 0.2e1 + t47 ^ 2 / 0.2e1) * t40 ^ 2) * t49, t44 ^ 2 * t62, t44 * t48 * t46, t44 * t52, t46 ^ 2 * t62, t46 * t52, qJD(3) ^ 2 / 0.2e1, t23 * qJD(3) - t33 * t54, -t24 * qJD(3) + t33 * t55 (-t23 * t44 + t24 * t46) * qJD(2), t24 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1, t30 ^ 2 / 0.2e1, -t30 * t28, t30 * t39, t28 ^ 2 / 0.2e1, -t28 * t39, t39 ^ 2 / 0.2e1, t26 * t28 + t9 * t39, -t10 * t39 + t26 * t30, -t10 * t28 - t9 * t30, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1, t19, -t59, t13, t53, -t58, t25, t7 * t20 + t3 * t27, t7 * t22 - t4 * t27, -t4 * t20 - t3 * t22, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t19, t13, t59, t25, t58, t53, -t1 * t27 + t5 * t20, t1 * t22 - t2 * t20, t2 * t27 - t5 * t22, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
