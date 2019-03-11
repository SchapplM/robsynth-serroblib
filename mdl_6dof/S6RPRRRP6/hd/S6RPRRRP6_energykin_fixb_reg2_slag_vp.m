% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRP6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP6_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:16:43
% EndTime: 2019-03-09 06:16:43
% DurationCPUTime: 0.19s
% Computational Cost: add. (794->63), mult. (1958->135), div. (0->0), fcn. (1443->8), ass. (0->49)
t58 = qJD(1) ^ 2;
t64 = t58 / 0.2e1;
t55 = sin(qJ(3));
t57 = cos(qJ(3));
t52 = cos(pkin(10));
t60 = qJD(1) * t52;
t51 = sin(pkin(10));
t61 = qJD(1) * t51;
t40 = t55 * t61 - t57 * t60;
t42 = (t51 * t57 + t52 * t55) * qJD(1);
t45 = qJD(2) + (-pkin(2) * t52 - pkin(1)) * qJD(1);
t24 = t40 * pkin(3) - t42 * pkin(8) + t45;
t62 = pkin(7) + qJ(2);
t43 = t62 * t61;
t44 = t62 * t60;
t29 = -t55 * t43 + t57 * t44;
t27 = qJD(3) * pkin(8) + t29;
t54 = sin(qJ(4));
t56 = cos(qJ(4));
t13 = t54 * t24 + t56 * t27;
t31 = -t56 * qJD(3) + t54 * t42;
t10 = -t31 * pkin(9) + t13;
t53 = sin(qJ(5));
t63 = cos(qJ(5));
t12 = t56 * t24 - t54 * t27;
t33 = t54 * qJD(3) + t56 * t42;
t36 = qJD(4) + t40;
t7 = t36 * pkin(4) - t33 * pkin(9) + t12;
t4 = t63 * t10 + t53 * t7;
t3 = -t53 * t10 + t63 * t7;
t28 = -t57 * t43 - t55 * t44;
t26 = -qJD(3) * pkin(3) - t28;
t16 = t31 * pkin(4) + t26;
t50 = t52 ^ 2;
t49 = t51 ^ 2;
t47 = -qJD(1) * pkin(1) + qJD(2);
t35 = qJD(5) + t36;
t34 = t35 ^ 2 / 0.2e1;
t21 = -t53 * t31 + t63 * t33;
t19 = t63 * t31 + t53 * t33;
t18 = t21 ^ 2 / 0.2e1;
t17 = t19 ^ 2 / 0.2e1;
t15 = t21 * t35;
t14 = t19 * t35;
t11 = t21 * t19;
t8 = t19 * pkin(5) + qJD(6) + t16;
t2 = -t19 * qJ(6) + t4;
t1 = t35 * pkin(5) - t21 * qJ(6) + t3;
t5 = [0, 0, 0, 0, 0, t64, 0, 0, 0, 0, t49 * t64, t51 * t58 * t52, 0, t50 * t64, 0, 0, -t47 * t60, t47 * t61 (t49 + t50) * t58 * qJ(2), t47 ^ 2 / 0.2e1 + (t50 / 0.2e1 + t49 / 0.2e1) * qJ(2) ^ 2 * t58, t42 ^ 2 / 0.2e1, -t42 * t40, t42 * qJD(3), t40 ^ 2 / 0.2e1, -t40 * qJD(3), qJD(3) ^ 2 / 0.2e1, t28 * qJD(3) + t45 * t40, -t29 * qJD(3) + t45 * t42, -t28 * t42 - t29 * t40, t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t45 ^ 2 / 0.2e1, t33 ^ 2 / 0.2e1, -t33 * t31, t33 * t36, t31 ^ 2 / 0.2e1, -t31 * t36, t36 ^ 2 / 0.2e1, t12 * t36 + t26 * t31, -t13 * t36 + t26 * t33, -t12 * t33 - t13 * t31, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1, t18, -t11, t15, t17, -t14, t34, t16 * t19 + t3 * t35, t16 * t21 - t4 * t35, -t4 * t19 - t3 * t21, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1, t18, -t11, t15, t17, -t14, t34, t1 * t35 + t8 * t19, -t2 * t35 + t8 * t21, -t1 * t21 - t2 * t19, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1;];
T_reg  = t5;
