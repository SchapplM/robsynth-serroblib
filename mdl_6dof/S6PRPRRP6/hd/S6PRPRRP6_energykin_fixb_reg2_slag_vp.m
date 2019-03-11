% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRRP6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP6_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:20:54
% EndTime: 2019-03-08 20:20:55
% DurationCPUTime: 0.15s
% Computational Cost: add. (251->46), mult. (534->102), div. (0->0), fcn. (320->8), ass. (0->46)
t38 = qJD(2) ^ 2;
t30 = t38 / 0.2e1;
t39 = qJD(1) ^ 2;
t54 = t39 / 0.2e1;
t34 = sin(qJ(4));
t36 = cos(qJ(4));
t35 = sin(qJ(2));
t31 = sin(pkin(6));
t50 = qJD(1) * t31;
t42 = t35 * t50;
t12 = t42 + (pkin(4) * t34 - pkin(9) * t36 + qJ(3)) * qJD(2);
t33 = sin(qJ(5));
t53 = cos(qJ(5));
t37 = cos(qJ(2));
t40 = -t37 * t50 + qJD(3);
t15 = (-pkin(2) - pkin(8)) * qJD(2) + t40;
t32 = cos(pkin(6));
t49 = qJD(1) * t32;
t10 = t34 * t15 + t36 * t49;
t8 = qJD(4) * pkin(9) + t10;
t4 = t33 * t12 + t53 * t8;
t48 = qJD(2) * t36;
t18 = -t53 * qJD(4) + t33 * t48;
t20 = t33 * qJD(4) + t53 * t48;
t52 = t20 * t18;
t46 = t34 * qJD(2);
t26 = qJD(5) + t46;
t51 = t26 * t18;
t21 = qJD(2) * qJ(3) + t42;
t47 = t21 * qJD(2);
t45 = t18 ^ 2 / 0.2e1;
t44 = t21 ^ 2 / 0.2e1;
t43 = qJD(2) * qJD(4);
t41 = qJD(2) * t50;
t9 = t36 * t15 - t34 * t49;
t7 = -qJD(4) * pkin(4) - t9;
t3 = t53 * t12 - t33 * t8;
t27 = t32 ^ 2 * t54;
t23 = t26 ^ 2 / 0.2e1;
t17 = -qJD(2) * pkin(2) + t40;
t16 = t20 ^ 2 / 0.2e1;
t13 = t20 * t26;
t5 = t18 * pkin(5) - t20 * qJ(6) + t7;
t2 = t26 * qJ(6) + t4;
t1 = -t26 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, 0, 0, 0, 0, t30, t37 * t41, -t35 * t41, 0, t27 + (t35 ^ 2 / 0.2e1 + t37 ^ 2 / 0.2e1) * t39 * t31 ^ 2, t30, 0, 0, 0, 0, 0, 0, t17 * qJD(2), t47, t27 + t44 + t17 ^ 2 / 0.2e1, t36 ^ 2 * t30, -t36 * t38 * t34, t36 * t43, t34 ^ 2 * t30, -t34 * t43, qJD(4) ^ 2 / 0.2e1, t9 * qJD(4) + t21 * t46, -t10 * qJD(4) + t36 * t47 (-t10 * t34 - t36 * t9) * qJD(2), t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t44, t16, -t52, t13, t45, -t51, t23, t7 * t18 + t3 * t26, t7 * t20 - t4 * t26, -t4 * t18 - t3 * t20, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t16, t13, t52, t23, t51, t45, -t1 * t26 + t5 * t18, t1 * t20 - t2 * t18, t2 * t26 - t5 * t20, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
