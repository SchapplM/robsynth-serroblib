% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPRR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR5_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:20:46
% EndTime: 2019-03-08 22:20:46
% DurationCPUTime: 0.18s
% Computational Cost: add. (693->61), mult. (1607->149), div. (0->0), fcn. (1185->12), ass. (0->51)
t56 = qJD(2) ^ 2;
t68 = t56 / 0.2e1;
t53 = sin(qJ(2));
t48 = sin(pkin(6));
t64 = qJD(1) * t48;
t37 = qJD(2) * pkin(8) + t53 * t64;
t52 = sin(qJ(3));
t54 = cos(qJ(3));
t49 = cos(pkin(6));
t63 = qJD(1) * t49;
t29 = t54 * t37 + t52 * t63;
t27 = qJD(3) * qJ(4) + t29;
t55 = cos(qJ(2));
t59 = t55 * t64;
t30 = -t59 + (-pkin(3) * t54 - qJ(4) * t52 - pkin(2)) * qJD(2);
t47 = sin(pkin(12));
t65 = cos(pkin(12));
t16 = -t47 * t27 + t65 * t30;
t62 = qJD(2) * t52;
t36 = t47 * qJD(3) + t65 * t62;
t61 = t54 * qJD(2);
t10 = -pkin(4) * t61 - t36 * pkin(9) + t16;
t17 = t65 * t27 + t47 * t30;
t34 = -t65 * qJD(3) + t47 * t62;
t15 = -t34 * pkin(9) + t17;
t51 = sin(qJ(5));
t67 = cos(qJ(5));
t6 = t51 * t10 + t67 * t15;
t66 = cos(qJ(6));
t60 = qJD(2) * qJD(3);
t5 = t67 * t10 - t51 * t15;
t58 = qJD(2) * t64;
t42 = -qJD(5) + t61;
t28 = -t52 * t37 + t54 * t63;
t24 = -qJD(3) * pkin(3) + qJD(4) - t28;
t18 = t34 * pkin(4) + t24;
t57 = qJD(1) ^ 2;
t50 = sin(qJ(6));
t44 = t54 ^ 2 * t68;
t39 = -qJD(6) + t42;
t38 = -qJD(2) * pkin(2) - t59;
t22 = -t51 * t34 + t67 * t36;
t20 = t67 * t34 + t51 * t36;
t13 = -t50 * t20 + t66 * t22;
t11 = t66 * t20 + t50 * t22;
t9 = t20 * pkin(5) + t18;
t4 = -t20 * pkin(10) + t6;
t3 = -t42 * pkin(5) - t22 * pkin(10) + t5;
t2 = t50 * t3 + t66 * t4;
t1 = t66 * t3 - t50 * t4;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t57 / 0.2e1, 0, 0, 0, 0, 0, t68, t55 * t58, -t53 * t58, 0 (t49 ^ 2 / 0.2e1 + (t53 ^ 2 / 0.2e1 + t55 ^ 2 / 0.2e1) * t48 ^ 2) * t57, t52 ^ 2 * t68, t52 * t56 * t54, t52 * t60, t44, t54 * t60, qJD(3) ^ 2 / 0.2e1, t28 * qJD(3) - t38 * t61, -t29 * qJD(3) + t38 * t62 (-t28 * t52 + t29 * t54) * qJD(2), t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1, t36 ^ 2 / 0.2e1, -t36 * t34, -t36 * t61, t34 ^ 2 / 0.2e1, t34 * t61, t44, -t16 * t61 + t24 * t34, t17 * t61 + t24 * t36, -t16 * t36 - t17 * t34, t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, -t22 * t42, t20 ^ 2 / 0.2e1, t20 * t42, t42 ^ 2 / 0.2e1, t18 * t20 - t5 * t42, t18 * t22 + t6 * t42, -t6 * t20 - t5 * t22, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t13 * t11, -t13 * t39, t11 ^ 2 / 0.2e1, t11 * t39, t39 ^ 2 / 0.2e1, -t1 * t39 + t9 * t11, t9 * t13 + t2 * t39, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1;];
T_reg  = t7;
