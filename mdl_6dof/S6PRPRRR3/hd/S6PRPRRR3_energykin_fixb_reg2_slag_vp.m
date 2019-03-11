% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRRR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:34:10
% EndTime: 2019-03-08 20:34:11
% DurationCPUTime: 0.17s
% Computational Cost: add. (619->60), mult. (1558->146), div. (0->0), fcn. (1211->12), ass. (0->50)
t52 = qJD(2) ^ 2;
t63 = t52 / 0.2e1;
t50 = sin(qJ(2));
t44 = sin(pkin(6));
t59 = qJD(1) * t44;
t36 = qJD(2) * qJ(3) + t50 * t59;
t45 = cos(pkin(12));
t46 = cos(pkin(6));
t58 = qJD(1) * t46;
t38 = t45 * t58;
t43 = sin(pkin(12));
t25 = t38 + (-pkin(8) * qJD(2) - t36) * t43;
t28 = t45 * t36 + t43 * t58;
t56 = qJD(2) * t45;
t26 = pkin(8) * t56 + t28;
t49 = sin(qJ(4));
t62 = cos(qJ(4));
t12 = t62 * t25 - t49 * t26;
t33 = (t62 * t43 + t45 * t49) * qJD(2);
t10 = qJD(4) * pkin(4) - t33 * pkin(9) + t12;
t13 = t49 * t25 + t62 * t26;
t57 = qJD(2) * t43;
t31 = t49 * t57 - t62 * t56;
t11 = -t31 * pkin(9) + t13;
t48 = sin(qJ(5));
t61 = cos(qJ(5));
t6 = t48 * t10 + t61 * t11;
t60 = cos(qJ(6));
t55 = qJD(2) * t59;
t18 = t61 * t31 + t48 * t33;
t51 = cos(qJ(2));
t54 = -t51 * t59 + qJD(3);
t5 = t61 * t10 - t48 * t11;
t30 = (-pkin(3) * t45 - pkin(2)) * qJD(2) + t54;
t21 = t31 * pkin(4) + t30;
t53 = qJD(1) ^ 2;
t47 = sin(qJ(6));
t42 = qJD(4) + qJD(5);
t35 = -qJD(2) * pkin(2) + t54;
t27 = -t43 * t36 + t38;
t20 = -t48 * t31 + t61 * t33;
t17 = qJD(6) + t18;
t16 = t60 * t20 + t47 * t42;
t14 = t47 * t20 - t60 * t42;
t7 = t18 * pkin(5) - t20 * pkin(10) + t21;
t4 = t42 * pkin(10) + t6;
t3 = -t42 * pkin(5) - t5;
t2 = t60 * t4 + t47 * t7;
t1 = -t47 * t4 + t60 * t7;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t53 / 0.2e1, 0, 0, 0, 0, 0, t63, t51 * t55, -t50 * t55, 0 (t46 ^ 2 / 0.2e1 + (t50 ^ 2 / 0.2e1 + t51 ^ 2 / 0.2e1) * t44 ^ 2) * t53, t43 ^ 2 * t63, t43 * t52 * t45, 0, t45 ^ 2 * t63, 0, 0, -t35 * t56, t35 * t57 (-t27 * t43 + t28 * t45) * qJD(2), t28 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1 + t35 ^ 2 / 0.2e1, t33 ^ 2 / 0.2e1, -t33 * t31, t33 * qJD(4), t31 ^ 2 / 0.2e1, -t31 * qJD(4), qJD(4) ^ 2 / 0.2e1, t12 * qJD(4) + t30 * t31, -t13 * qJD(4) + t30 * t33, -t12 * t33 - t13 * t31, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1, t20 ^ 2 / 0.2e1, -t20 * t18, t20 * t42, t18 ^ 2 / 0.2e1, -t18 * t42, t42 ^ 2 / 0.2e1, t21 * t18 + t5 * t42, t21 * t20 - t6 * t42, -t6 * t18 - t5 * t20, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t16 ^ 2 / 0.2e1, -t16 * t14, t16 * t17, t14 ^ 2 / 0.2e1, -t14 * t17, t17 ^ 2 / 0.2e1, t1 * t17 + t3 * t14, t3 * t16 - t2 * t17, -t1 * t16 - t2 * t14, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
