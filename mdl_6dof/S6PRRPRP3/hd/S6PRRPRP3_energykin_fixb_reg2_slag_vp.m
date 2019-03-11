% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPRP3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:38:19
% EndTime: 2019-03-08 21:38:19
% DurationCPUTime: 0.15s
% Computational Cost: add. (489->57), mult. (1150->129), div. (0->0), fcn. (815->10), ass. (0->50)
t49 = qJD(2) ^ 2;
t63 = t49 / 0.2e1;
t44 = sin(qJ(5));
t62 = cos(qJ(5));
t46 = sin(qJ(2));
t42 = sin(pkin(6));
t58 = qJD(1) * t42;
t32 = qJD(2) * pkin(8) + t46 * t58;
t45 = sin(qJ(3));
t47 = cos(qJ(3));
t43 = cos(pkin(6));
t57 = qJD(1) * t43;
t24 = t47 * t32 + t45 * t57;
t22 = qJD(3) * qJ(4) + t24;
t48 = cos(qJ(2));
t52 = t48 * t58;
t25 = -t52 + (-pkin(3) * t47 - qJ(4) * t45 - pkin(2)) * qJD(2);
t41 = sin(pkin(11));
t59 = cos(pkin(11));
t10 = -t41 * t22 + t59 * t25;
t56 = qJD(2) * t45;
t31 = t41 * qJD(3) + t59 * t56;
t55 = t47 * qJD(2);
t7 = -pkin(4) * t55 - t31 * pkin(9) + t10;
t11 = t59 * t22 + t41 * t25;
t29 = -t59 * qJD(3) + t41 * t56;
t9 = -t29 * pkin(9) + t11;
t4 = t44 * t7 + t62 * t9;
t15 = t62 * t29 + t44 * t31;
t17 = -t44 * t29 + t62 * t31;
t61 = t17 * t15;
t36 = -qJD(5) + t55;
t60 = t36 * t15;
t54 = t15 ^ 2 / 0.2e1;
t53 = qJD(2) * qJD(3);
t51 = qJD(2) * t58;
t23 = -t45 * t32 + t47 * t57;
t3 = -t44 * t9 + t62 * t7;
t19 = -qJD(3) * pkin(3) + qJD(4) - t23;
t12 = t29 * pkin(4) + t19;
t50 = qJD(1) ^ 2;
t38 = t47 ^ 2 * t63;
t34 = t36 ^ 2 / 0.2e1;
t33 = -qJD(2) * pkin(2) - t52;
t14 = t17 ^ 2 / 0.2e1;
t13 = t17 * t36;
t5 = t15 * pkin(5) - t17 * qJ(6) + t12;
t2 = -t36 * qJ(6) + t4;
t1 = t36 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t50 / 0.2e1, 0, 0, 0, 0, 0, t63, t48 * t51, -t46 * t51, 0 (t43 ^ 2 / 0.2e1 + (t46 ^ 2 / 0.2e1 + t48 ^ 2 / 0.2e1) * t42 ^ 2) * t50, t45 ^ 2 * t63, t45 * t49 * t47, t45 * t53, t38, t47 * t53, qJD(3) ^ 2 / 0.2e1, t23 * qJD(3) - t33 * t55, -t24 * qJD(3) + t33 * t56 (-t23 * t45 + t24 * t47) * qJD(2), t24 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1, t31 ^ 2 / 0.2e1, -t31 * t29, -t31 * t55, t29 ^ 2 / 0.2e1, t29 * t55, t38, -t10 * t55 + t19 * t29, t11 * t55 + t19 * t31, -t10 * t31 - t11 * t29, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t14, -t61, -t13, t54, t60, t34, t12 * t15 - t3 * t36, t12 * t17 + t4 * t36, -t4 * t15 - t3 * t17, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1, t14, -t13, t61, t34, -t60, t54, t1 * t36 + t5 * t15, t1 * t17 - t2 * t15, -t5 * t17 - t2 * t36, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
