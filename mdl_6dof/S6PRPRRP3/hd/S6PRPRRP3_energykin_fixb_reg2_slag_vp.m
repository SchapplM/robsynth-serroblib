% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRPRRP3
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
% Datum: 2019-03-08 20:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRRP3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:07:44
% EndTime: 2019-03-08 20:07:44
% DurationCPUTime: 0.14s
% Computational Cost: add. (436->56), mult. (1094->126), div. (0->0), fcn. (810->10), ass. (0->49)
t51 = qJD(2) ^ 2;
t61 = t51 / 0.2e1;
t45 = cos(pkin(11));
t50 = cos(qJ(2));
t44 = sin(pkin(6));
t58 = qJD(1) * t44;
t53 = -t50 * t58 + qJD(3);
t30 = (-pkin(3) * t45 - pkin(2)) * qJD(2) + t53;
t48 = sin(qJ(4));
t55 = qJD(2) * t45;
t43 = sin(pkin(11));
t56 = qJD(2) * t43;
t60 = cos(qJ(4));
t32 = t48 * t56 - t60 * t55;
t34 = (t60 * t43 + t45 * t48) * qJD(2);
t13 = t32 * pkin(4) - t34 * pkin(9) + t30;
t47 = sin(qJ(5));
t59 = cos(qJ(5));
t49 = sin(qJ(2));
t37 = qJD(2) * qJ(3) + t49 * t58;
t46 = cos(pkin(6));
t57 = qJD(1) * t46;
t39 = t45 * t57;
t20 = t39 + (-pkin(8) * qJD(2) - t37) * t43;
t28 = t45 * t37 + t43 * t57;
t21 = pkin(8) * t55 + t28;
t10 = t48 * t20 + t60 * t21;
t8 = qJD(4) * pkin(9) + t10;
t4 = t47 * t13 + t59 * t8;
t3 = t59 * t13 - t47 * t8;
t54 = qJD(2) * t58;
t9 = t60 * t20 - t48 * t21;
t7 = -qJD(4) * pkin(4) - t9;
t52 = qJD(1) ^ 2;
t36 = -qJD(2) * pkin(2) + t53;
t31 = qJD(5) + t32;
t29 = t31 ^ 2 / 0.2e1;
t27 = -t43 * t37 + t39;
t26 = t47 * qJD(4) + t59 * t34;
t24 = -t59 * qJD(4) + t47 * t34;
t23 = t26 ^ 2 / 0.2e1;
t22 = t24 ^ 2 / 0.2e1;
t16 = t26 * t31;
t15 = t24 * t31;
t14 = t26 * t24;
t5 = t24 * pkin(5) + qJD(6) + t7;
t2 = -t24 * qJ(6) + t4;
t1 = t31 * pkin(5) - t26 * qJ(6) + t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t52 / 0.2e1, 0, 0, 0, 0, 0, t61, t50 * t54, -t49 * t54, 0 (t46 ^ 2 / 0.2e1 + (t49 ^ 2 / 0.2e1 + t50 ^ 2 / 0.2e1) * t44 ^ 2) * t52, t43 ^ 2 * t61, t43 * t51 * t45, 0, t45 ^ 2 * t61, 0, 0, -t36 * t55, t36 * t56 (-t27 * t43 + t28 * t45) * qJD(2), t28 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1 + t36 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, t34 * qJD(4), t32 ^ 2 / 0.2e1, -t32 * qJD(4), qJD(4) ^ 2 / 0.2e1, t9 * qJD(4) + t30 * t32, -t10 * qJD(4) + t30 * t34, -t10 * t32 - t9 * t34, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1, t23, -t14, t16, t22, -t15, t29, t7 * t24 + t3 * t31, t7 * t26 - t4 * t31, -t4 * t24 - t3 * t26, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t23, -t14, t16, t22, -t15, t29, t1 * t31 + t5 * t24, -t2 * t31 + t5 * t26, -t1 * t26 - t2 * t24, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t6;
