% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRP1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP1_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:33:00
% EndTime: 2019-03-09 16:33:00
% DurationCPUTime: 0.20s
% Computational Cost: add. (878->64), mult. (2110->140), div. (0->0), fcn. (1538->8), ass. (0->53)
t56 = qJD(1) ^ 2;
t68 = t56 / 0.2e1;
t67 = -pkin(8) - pkin(7);
t53 = sin(qJ(3));
t55 = cos(qJ(2));
t61 = qJD(1) * t55;
t54 = sin(qJ(2));
t62 = qJD(1) * t54;
t66 = cos(qJ(3));
t38 = t53 * t62 - t66 * t61;
t40 = (t53 * t55 + t66 * t54) * qJD(1);
t51 = sin(pkin(10));
t63 = cos(pkin(10));
t29 = t63 * t38 + t51 * t40;
t31 = -t51 * t38 + t63 * t40;
t44 = (-pkin(2) * t55 - pkin(1)) * qJD(1);
t34 = t38 * pkin(3) + qJD(4) + t44;
t13 = t29 * pkin(4) - t31 * pkin(9) + t34;
t52 = sin(qJ(5));
t65 = cos(qJ(5));
t42 = qJD(2) * pkin(2) + t67 * t62;
t43 = t67 * t61;
t32 = t66 * t42 + t53 * t43;
t48 = qJD(2) + qJD(3);
t18 = t48 * pkin(3) - t40 * qJ(4) + t32;
t33 = t53 * t42 - t66 * t43;
t23 = -t38 * qJ(4) + t33;
t10 = t51 * t18 + t63 * t23;
t8 = t48 * pkin(9) + t10;
t4 = t52 * t13 + t65 * t8;
t64 = t55 * t56;
t60 = qJD(1) * qJD(2);
t3 = t65 * t13 - t52 * t8;
t59 = t54 * t60;
t58 = t55 * t60;
t9 = t63 * t18 - t51 * t23;
t7 = -t48 * pkin(4) - t9;
t50 = t55 ^ 2;
t49 = t54 ^ 2;
t47 = t48 ^ 2 / 0.2e1;
t28 = qJD(5) + t29;
t27 = t28 ^ 2 / 0.2e1;
t26 = t65 * t31 + t52 * t48;
t24 = t52 * t31 - t65 * t48;
t22 = t26 ^ 2 / 0.2e1;
t21 = t24 ^ 2 / 0.2e1;
t16 = t26 * t28;
t15 = t24 * t28;
t14 = t26 * t24;
t5 = t24 * pkin(5) + qJD(6) + t7;
t2 = -t24 * qJ(6) + t4;
t1 = t28 * pkin(5) - t26 * qJ(6) + t3;
t6 = [0, 0, 0, 0, 0, t68, 0, 0, 0, 0, t49 * t68, t54 * t64, t59, t50 * t68, t58, qJD(2) ^ 2 / 0.2e1, pkin(1) * t64 - pkin(7) * t59, -t56 * pkin(1) * t54 - pkin(7) * t58 (t49 + t50) * t56 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t50 / 0.2e1 + t49 / 0.2e1) * pkin(7) ^ 2) * t56, t40 ^ 2 / 0.2e1, -t40 * t38, t40 * t48, t38 ^ 2 / 0.2e1, -t38 * t48, t47, t32 * t48 + t44 * t38, -t33 * t48 + t44 * t40, -t32 * t40 - t33 * t38, t33 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1 + t44 ^ 2 / 0.2e1, t31 ^ 2 / 0.2e1, -t31 * t29, t31 * t48, t29 ^ 2 / 0.2e1, -t29 * t48, t47, t34 * t29 + t9 * t48, -t10 * t48 + t34 * t31, -t10 * t29 - t9 * t31, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1, t22, -t14, t16, t21, -t15, t27, t7 * t24 + t3 * t28, t7 * t26 - t4 * t28, -t4 * t24 - t3 * t26, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t22, -t14, t16, t21, -t15, t27, t1 * t28 + t5 * t24, -t2 * t28 + t5 * t26, -t1 * t26 - t2 * t24, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t6;
