% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRPRRP1
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
% Datum: 2019-03-08 19:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRRP1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:58:28
% EndTime: 2019-03-08 19:58:28
% DurationCPUTime: 0.12s
% Computational Cost: add. (306->48), mult. (726->112), div. (0->0), fcn. (498->10), ass. (0->45)
t49 = qJD(2) ^ 2;
t39 = t49 / 0.2e1;
t48 = cos(qJ(2));
t41 = sin(pkin(6));
t56 = qJD(1) * t41;
t29 = qJD(2) * pkin(2) + t48 * t56;
t40 = sin(pkin(11));
t42 = cos(pkin(11));
t46 = sin(qJ(2));
t52 = t46 * t56;
t18 = t42 * t29 - t40 * t52;
t45 = sin(qJ(4));
t47 = cos(qJ(4));
t13 = (-pkin(4) * t47 - pkin(9) * t45 - pkin(3)) * qJD(2) - t18;
t44 = sin(qJ(5));
t57 = cos(qJ(5));
t19 = t40 * t29 + t42 * t52;
t17 = qJD(2) * pkin(8) + t19;
t43 = cos(pkin(6));
t34 = t43 * qJD(1) + qJD(3);
t10 = t47 * t17 + t45 * t34;
t8 = qJD(4) * pkin(9) + t10;
t4 = t44 * t13 + t57 * t8;
t55 = qJD(2) * t45;
t54 = t47 * qJD(2);
t53 = qJD(2) * qJD(4);
t3 = t57 * t13 - t44 * t8;
t51 = qJD(2) * t56;
t9 = -t45 * t17 + t47 * t34;
t7 = -qJD(4) * pkin(4) - t9;
t50 = qJD(1) ^ 2;
t35 = -qJD(5) + t54;
t33 = t35 ^ 2 / 0.2e1;
t28 = t44 * qJD(4) + t57 * t55;
t26 = -t57 * qJD(4) + t44 * t55;
t25 = t28 ^ 2 / 0.2e1;
t24 = t26 ^ 2 / 0.2e1;
t22 = t28 * t35;
t21 = t26 * t35;
t20 = t28 * t26;
t16 = -qJD(2) * pkin(3) - t18;
t5 = t26 * pkin(5) + qJD(6) + t7;
t2 = -t26 * qJ(6) + t4;
t1 = -t35 * pkin(5) - t28 * qJ(6) + t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t50 / 0.2e1, 0, 0, 0, 0, 0, t39, t48 * t51, -t46 * t51, 0 (t43 ^ 2 / 0.2e1 + (t46 ^ 2 / 0.2e1 + t48 ^ 2 / 0.2e1) * t41 ^ 2) * t50, 0, 0, 0, 0, 0, t39, t18 * qJD(2), -t19 * qJD(2), 0, t19 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1, t45 ^ 2 * t39, t45 * t49 * t47, t45 * t53, t47 ^ 2 * t39, t47 * t53, qJD(4) ^ 2 / 0.2e1, t9 * qJD(4) - t16 * t54, -t10 * qJD(4) + t16 * t55 (t10 * t47 - t45 * t9) * qJD(2), t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1, t25, -t20, -t22, t24, t21, t33, t7 * t26 - t3 * t35, t7 * t28 + t4 * t35, -t4 * t26 - t3 * t28, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t25, -t20, -t22, t24, t21, t33, -t1 * t35 + t5 * t26, t2 * t35 + t5 * t28, -t1 * t28 - t2 * t26, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t6;
