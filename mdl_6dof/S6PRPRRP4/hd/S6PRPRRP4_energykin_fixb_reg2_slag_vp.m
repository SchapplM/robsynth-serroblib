% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRRP4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP4_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:12:17
% EndTime: 2019-03-08 20:12:17
% DurationCPUTime: 0.14s
% Computational Cost: add. (436->54), mult. (1088->126), div. (0->0), fcn. (804->10), ass. (0->49)
t47 = qJD(2) ^ 2;
t60 = t47 / 0.2e1;
t41 = cos(pkin(11));
t46 = cos(qJ(2));
t40 = sin(pkin(6));
t55 = qJD(1) * t40;
t49 = -t46 * t55 + qJD(3);
t26 = (-pkin(3) * t41 - pkin(2)) * qJD(2) + t49;
t44 = sin(qJ(4));
t52 = qJD(2) * t41;
t39 = sin(pkin(11));
t53 = qJD(2) * t39;
t59 = cos(qJ(4));
t28 = t44 * t53 - t59 * t52;
t30 = (t59 * t39 + t41 * t44) * qJD(2);
t12 = t28 * pkin(4) - t30 * pkin(9) + t26;
t43 = sin(qJ(5));
t58 = cos(qJ(5));
t45 = sin(qJ(2));
t33 = qJD(2) * qJ(3) + t45 * t55;
t42 = cos(pkin(6));
t54 = qJD(1) * t42;
t35 = t41 * t54;
t17 = t35 + (-pkin(8) * qJD(2) - t33) * t39;
t24 = t41 * t33 + t39 * t54;
t18 = pkin(8) * t52 + t24;
t10 = t44 * t17 + t59 * t18;
t8 = qJD(4) * pkin(9) + t10;
t4 = t43 * t12 + t58 * t8;
t20 = -t58 * qJD(4) + t43 * t30;
t22 = t43 * qJD(4) + t58 * t30;
t57 = t22 * t20;
t27 = qJD(5) + t28;
t56 = t27 * t20;
t51 = t20 ^ 2 / 0.2e1;
t50 = qJD(2) * t55;
t9 = t59 * t17 - t44 * t18;
t3 = t58 * t12 - t43 * t8;
t7 = -qJD(4) * pkin(4) - t9;
t48 = qJD(1) ^ 2;
t32 = -qJD(2) * pkin(2) + t49;
t25 = t27 ^ 2 / 0.2e1;
t23 = -t39 * t33 + t35;
t19 = t22 ^ 2 / 0.2e1;
t13 = t22 * t27;
t5 = t20 * pkin(5) - t22 * qJ(6) + t7;
t2 = t27 * qJ(6) + t4;
t1 = -t27 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t48 / 0.2e1, 0, 0, 0, 0, 0, t60, t46 * t50, -t45 * t50, 0 (t42 ^ 2 / 0.2e1 + (t45 ^ 2 / 0.2e1 + t46 ^ 2 / 0.2e1) * t40 ^ 2) * t48, t39 ^ 2 * t60, t39 * t47 * t41, 0, t41 ^ 2 * t60, 0, 0, -t32 * t52, t32 * t53 (-t23 * t39 + t24 * t41) * qJD(2), t24 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1, t30 ^ 2 / 0.2e1, -t30 * t28, t30 * qJD(4), t28 ^ 2 / 0.2e1, -t28 * qJD(4), qJD(4) ^ 2 / 0.2e1, t9 * qJD(4) + t26 * t28, -t10 * qJD(4) + t26 * t30, -t10 * t28 - t9 * t30, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1, t19, -t57, t13, t51, -t56, t25, t7 * t20 + t3 * t27, t7 * t22 - t4 * t27, -t4 * t20 - t3 * t22, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t19, t13, t57, t25, t56, t51, -t1 * t27 + t5 * t20, t1 * t22 - t2 * t20, t2 * t27 - t5 * t22, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
