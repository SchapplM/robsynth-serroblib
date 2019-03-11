% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRRP1
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
% Datum: 2019-03-09 11:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRP1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP1_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:41:13
% EndTime: 2019-03-09 11:41:14
% DurationCPUTime: 0.21s
% Computational Cost: add. (852->64), mult. (2110->140), div. (0->0), fcn. (1538->8), ass. (0->53)
t56 = qJD(1) ^ 2;
t68 = t56 / 0.2e1;
t51 = sin(pkin(10));
t55 = cos(qJ(2));
t61 = qJD(1) * t55;
t54 = sin(qJ(2));
t62 = qJD(1) * t54;
t63 = cos(pkin(10));
t38 = t51 * t62 - t63 * t61;
t40 = (t51 * t55 + t63 * t54) * qJD(1);
t53 = sin(qJ(4));
t67 = cos(qJ(4));
t29 = t67 * t38 + t53 * t40;
t31 = -t53 * t38 + t67 * t40;
t44 = qJD(3) + (-pkin(2) * t55 - pkin(1)) * qJD(1);
t34 = t38 * pkin(3) + t44;
t13 = t29 * pkin(4) - t31 * pkin(9) + t34;
t52 = sin(qJ(5));
t66 = cos(qJ(5));
t64 = pkin(7) + qJ(3);
t42 = qJD(2) * pkin(2) - t64 * t62;
t43 = t64 * t61;
t32 = t63 * t42 - t51 * t43;
t20 = qJD(2) * pkin(3) - t40 * pkin(8) + t32;
t33 = t51 * t42 + t63 * t43;
t21 = -t38 * pkin(8) + t33;
t10 = t53 * t20 + t67 * t21;
t47 = qJD(2) + qJD(4);
t8 = t47 * pkin(9) + t10;
t4 = t52 * t13 + t66 * t8;
t65 = t55 * t56;
t60 = qJD(1) * qJD(2);
t3 = t66 * t13 - t52 * t8;
t59 = t54 * t60;
t58 = t55 * t60;
t9 = t67 * t20 - t53 * t21;
t7 = -t47 * pkin(4) - t9;
t50 = t55 ^ 2;
t49 = t54 ^ 2;
t48 = qJD(2) ^ 2 / 0.2e1;
t28 = qJD(5) + t29;
t27 = t28 ^ 2 / 0.2e1;
t26 = t66 * t31 + t52 * t47;
t24 = t52 * t31 - t66 * t47;
t23 = t26 ^ 2 / 0.2e1;
t22 = t24 ^ 2 / 0.2e1;
t16 = t26 * t28;
t15 = t24 * t28;
t14 = t26 * t24;
t5 = t24 * pkin(5) + qJD(6) + t7;
t2 = -t24 * qJ(6) + t4;
t1 = t28 * pkin(5) - t26 * qJ(6) + t3;
t6 = [0, 0, 0, 0, 0, t68, 0, 0, 0, 0, t49 * t68, t54 * t65, t59, t50 * t68, t58, t48, pkin(1) * t65 - pkin(7) * t59, -t56 * pkin(1) * t54 - pkin(7) * t58 (t49 + t50) * t56 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t50 / 0.2e1 + t49 / 0.2e1) * pkin(7) ^ 2) * t56, t40 ^ 2 / 0.2e1, -t40 * t38, t40 * qJD(2), t38 ^ 2 / 0.2e1, -t38 * qJD(2), t48, t32 * qJD(2) + t44 * t38, -t33 * qJD(2) + t44 * t40, -t32 * t40 - t33 * t38, t33 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1 + t44 ^ 2 / 0.2e1, t31 ^ 2 / 0.2e1, -t31 * t29, t31 * t47, t29 ^ 2 / 0.2e1, -t29 * t47, t47 ^ 2 / 0.2e1, t34 * t29 + t9 * t47, -t10 * t47 + t34 * t31, -t10 * t29 - t9 * t31, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1, t23, -t14, t16, t22, -t15, t27, t7 * t24 + t3 * t28, t7 * t26 - t4 * t28, -t4 * t24 - t3 * t26, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t23, -t14, t16, t22, -t15, t27, t1 * t28 + t5 * t24, -t2 * t28 + t5 * t26, -t1 * t26 - t2 * t24, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t6;
