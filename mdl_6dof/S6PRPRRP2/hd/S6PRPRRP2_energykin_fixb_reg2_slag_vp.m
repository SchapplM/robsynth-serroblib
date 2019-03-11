% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRPRRP2
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
% Datum: 2019-03-08 20:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRRP2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP2_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:03:04
% EndTime: 2019-03-08 20:03:05
% DurationCPUTime: 0.12s
% Computational Cost: add. (306->48), mult. (720->112), div. (0->0), fcn. (492->10), ass. (0->45)
t45 = qJD(2) ^ 2;
t35 = t45 / 0.2e1;
t44 = cos(qJ(2));
t37 = sin(pkin(6));
t53 = qJD(1) * t37;
t25 = qJD(2) * pkin(2) + t44 * t53;
t36 = sin(pkin(11));
t38 = cos(pkin(11));
t42 = sin(qJ(2));
t48 = t42 * t53;
t17 = t38 * t25 - t36 * t48;
t41 = sin(qJ(4));
t43 = cos(qJ(4));
t12 = (-pkin(4) * t43 - pkin(9) * t41 - pkin(3)) * qJD(2) - t17;
t40 = sin(qJ(5));
t56 = cos(qJ(5));
t18 = t36 * t25 + t38 * t48;
t16 = qJD(2) * pkin(8) + t18;
t39 = cos(pkin(6));
t30 = t39 * qJD(1) + qJD(3);
t10 = t43 * t16 + t41 * t30;
t8 = qJD(4) * pkin(9) + t10;
t4 = t40 * t12 + t56 * t8;
t52 = qJD(2) * t41;
t22 = -t56 * qJD(4) + t40 * t52;
t24 = t40 * qJD(4) + t56 * t52;
t55 = t24 * t22;
t51 = t43 * qJD(2);
t31 = -qJD(5) + t51;
t54 = t31 * t22;
t50 = t22 ^ 2 / 0.2e1;
t49 = qJD(2) * qJD(4);
t47 = qJD(2) * t53;
t9 = -t41 * t16 + t43 * t30;
t7 = -qJD(4) * pkin(4) - t9;
t3 = t56 * t12 - t40 * t8;
t46 = qJD(1) ^ 2;
t29 = t31 ^ 2 / 0.2e1;
t21 = t24 ^ 2 / 0.2e1;
t19 = t24 * t31;
t15 = -qJD(2) * pkin(3) - t17;
t5 = t22 * pkin(5) - t24 * qJ(6) + t7;
t2 = -t31 * qJ(6) + t4;
t1 = t31 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t46 / 0.2e1, 0, 0, 0, 0, 0, t35, t44 * t47, -t42 * t47, 0 (t39 ^ 2 / 0.2e1 + (t42 ^ 2 / 0.2e1 + t44 ^ 2 / 0.2e1) * t37 ^ 2) * t46, 0, 0, 0, 0, 0, t35, t17 * qJD(2), -t18 * qJD(2), 0, t18 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1, t41 ^ 2 * t35, t41 * t45 * t43, t41 * t49, t43 ^ 2 * t35, t43 * t49, qJD(4) ^ 2 / 0.2e1, t9 * qJD(4) - t15 * t51, -t10 * qJD(4) + t15 * t52 (t10 * t43 - t41 * t9) * qJD(2), t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t21, -t55, -t19, t50, t54, t29, t7 * t22 - t3 * t31, t7 * t24 + t4 * t31, -t4 * t22 - t3 * t24, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t21, -t19, t55, t29, -t54, t50, t1 * t31 + t5 * t22, t1 * t24 - t2 * t22, -t2 * t31 - t5 * t24, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
