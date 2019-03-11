% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRP8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP8_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:25:14
% EndTime: 2019-03-09 12:25:14
% DurationCPUTime: 0.22s
% Computational Cost: add. (943->65), mult. (2166->142), div. (0->0), fcn. (1533->8), ass. (0->52)
t55 = qJD(1) ^ 2;
t69 = t55 / 0.2e1;
t51 = sin(qJ(5));
t67 = cos(qJ(5));
t53 = sin(qJ(2));
t54 = cos(qJ(2));
t33 = (-pkin(2) * t54 - qJ(3) * t53 - pkin(1)) * qJD(1);
t61 = t54 * qJD(1);
t40 = pkin(7) * t61 + qJD(2) * qJ(3);
t50 = sin(pkin(10));
t63 = cos(pkin(10));
t27 = t63 * t33 - t50 * t40;
t62 = qJD(1) * t53;
t36 = t50 * qJD(2) + t63 * t62;
t20 = -pkin(3) * t61 - t36 * pkin(8) + t27;
t28 = t50 * t33 + t63 * t40;
t34 = -t63 * qJD(2) + t50 * t62;
t22 = -t34 * pkin(8) + t28;
t52 = sin(qJ(4));
t68 = cos(qJ(4));
t10 = t68 * t20 - t52 * t22;
t26 = -t52 * t34 + t68 * t36;
t43 = -qJD(4) + t61;
t7 = -t43 * pkin(4) - t26 * pkin(9) + t10;
t11 = t52 * t20 + t68 * t22;
t24 = t68 * t34 + t52 * t36;
t9 = -t24 * pkin(9) + t11;
t4 = t51 * t7 + t67 * t9;
t14 = t67 * t24 + t51 * t26;
t16 = -t51 * t24 + t67 * t26;
t66 = t16 * t14;
t41 = -qJD(5) + t43;
t65 = t41 * t14;
t64 = t54 * t55;
t60 = t14 ^ 2 / 0.2e1;
t59 = qJD(1) * qJD(2);
t58 = t53 * t59;
t57 = t54 * t59;
t38 = -qJD(2) * pkin(2) + pkin(7) * t62 + qJD(3);
t3 = -t51 * t9 + t67 * t7;
t29 = t34 * pkin(3) + t38;
t17 = t24 * pkin(4) + t29;
t49 = t54 ^ 2;
t48 = t53 ^ 2;
t45 = t49 * t69;
t39 = t41 ^ 2 / 0.2e1;
t13 = t16 ^ 2 / 0.2e1;
t12 = t16 * t41;
t5 = t14 * pkin(5) - t16 * qJ(6) + t17;
t2 = -t41 * qJ(6) + t4;
t1 = t41 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t69, 0, 0, 0, 0, t48 * t69, t53 * t64, t58, t45, t57, qJD(2) ^ 2 / 0.2e1, pkin(1) * t64 - pkin(7) * t58, -t55 * pkin(1) * t53 - pkin(7) * t57 (t48 + t49) * t55 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t49 / 0.2e1 + t48 / 0.2e1) * pkin(7) ^ 2) * t55, t36 ^ 2 / 0.2e1, -t36 * t34, -t36 * t61, t34 ^ 2 / 0.2e1, t34 * t61, t45, -t27 * t61 + t38 * t34, t28 * t61 + t38 * t36, -t27 * t36 - t28 * t34, t28 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1, t26 ^ 2 / 0.2e1, -t26 * t24, -t26 * t43, t24 ^ 2 / 0.2e1, t24 * t43, t43 ^ 2 / 0.2e1, -t10 * t43 + t29 * t24, t11 * t43 + t29 * t26, -t10 * t26 - t11 * t24, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1, t13, -t66, -t12, t60, t65, t39, t17 * t14 - t3 * t41, t17 * t16 + t4 * t41, -t4 * t14 - t3 * t16, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t13, -t12, t66, t39, -t65, t60, t1 * t41 + t5 * t14, t1 * t16 - t2 * t14, -t5 * t16 - t2 * t41, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
