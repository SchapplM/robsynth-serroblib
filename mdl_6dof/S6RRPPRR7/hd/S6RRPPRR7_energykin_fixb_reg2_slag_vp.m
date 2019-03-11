% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPRR7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR7_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:20:57
% EndTime: 2019-03-09 09:20:57
% DurationCPUTime: 0.20s
% Computational Cost: add. (556->66), mult. (1328->140), div. (0->0), fcn. (875->8), ass. (0->52)
t67 = -pkin(2) - pkin(3);
t49 = cos(pkin(6));
t62 = t49 * qJD(1);
t44 = qJD(2) + t62;
t51 = sin(qJ(2));
t48 = sin(pkin(6));
t63 = qJD(1) * t48;
t58 = t51 * t63;
t39 = pkin(8) * t58;
t53 = cos(qJ(2));
t55 = qJD(3) + t39 + (-pkin(1) * t49 * t53 - qJ(4) * t48 * t51) * qJD(1);
t10 = (-pkin(9) + t67) * t44 + t55;
t59 = t53 * t63;
t21 = -pkin(1) * t63 - pkin(2) * t59 - qJ(3) * t58;
t18 = pkin(3) * t59 + qJD(4) - t21;
t11 = (pkin(4) * t51 + pkin(9) * t53) * t63 + t18;
t52 = cos(qJ(5));
t65 = sin(qJ(5));
t6 = t52 * t10 + t65 * t11;
t66 = cos(qJ(6));
t54 = qJD(1) ^ 2;
t64 = t48 ^ 2 * t54;
t60 = pkin(1) * t62;
t28 = pkin(8) * t59 + t51 * t60;
t34 = t44 ^ 2 / 0.2e1;
t61 = t53 * t64;
t20 = t44 * qJ(3) + t28;
t57 = t64 / 0.2e1;
t29 = t44 * t58;
t56 = t44 * t59;
t27 = t53 * t60 - t39;
t5 = -t65 * t10 + t52 * t11;
t17 = qJ(4) * t59 - t20;
t24 = -t52 * t44 + t65 * t59;
t13 = t44 * pkin(4) - t17;
t50 = sin(qJ(6));
t36 = t53 ^ 2 * t57;
t35 = t51 ^ 2 * t57;
t33 = qJD(5) + t58;
t32 = t51 * t61;
t25 = t65 * t44 + t52 * t59;
t22 = -qJD(6) + t24;
t19 = -t44 * pkin(2) + qJD(3) - t27;
t16 = -t66 * t25 + t50 * t33;
t14 = -t50 * t25 - t66 * t33;
t12 = t67 * t44 + t55;
t7 = -t24 * pkin(5) + t25 * pkin(10) + t13;
t4 = t33 * pkin(10) + t6;
t3 = -t33 * pkin(5) - t5;
t2 = t66 * t4 + t50 * t7;
t1 = -t50 * t4 + t66 * t7;
t8 = [0, 0, 0, 0, 0, t54 / 0.2e1, 0, 0, 0, 0, t35, t32, t29, t36, t56, t34, pkin(1) * t61 + t27 * t44, -pkin(1) * t51 * t64 - t28 * t44 (-t27 * t51 + t28 * t53) * t63, t28 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t57, t35, t29, -t32, t34, -t56, t36, -t19 * t44 - t21 * t59 (t19 * t51 + t20 * t53) * t63, t20 * t44 - t21 * t58, t20 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t36, t32, t56, t35, t29, t34, -t17 * t44 + t18 * t58, t12 * t44 - t18 * t59 (-t12 * t51 + t17 * t53) * t63, t12 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t25 ^ 2 / 0.2e1, -t25 * t24, -t25 * t33, t24 ^ 2 / 0.2e1, t24 * t33, t33 ^ 2 / 0.2e1, -t13 * t24 + t5 * t33, -t13 * t25 - t6 * t33, t6 * t24 + t5 * t25, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1, t16 ^ 2 / 0.2e1, -t16 * t14, -t16 * t22, t14 ^ 2 / 0.2e1, t14 * t22, t22 ^ 2 / 0.2e1, -t1 * t22 + t3 * t14, t3 * t16 + t2 * t22, -t1 * t16 - t2 * t14, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
