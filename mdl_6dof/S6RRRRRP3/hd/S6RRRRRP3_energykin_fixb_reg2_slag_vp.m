% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRP3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:09:29
% EndTime: 2019-03-10 01:09:29
% DurationCPUTime: 0.22s
% Computational Cost: add. (946->64), mult. (2016->143), div. (0->0), fcn. (1447->8), ass. (0->53)
t58 = qJD(1) ^ 2;
t68 = t58 / 0.2e1;
t67 = -pkin(8) - pkin(7);
t53 = sin(qJ(3));
t56 = cos(qJ(3));
t57 = cos(qJ(2));
t63 = qJD(1) * t57;
t54 = sin(qJ(2));
t64 = qJD(1) * t54;
t38 = t53 * t64 - t56 * t63;
t40 = (t53 * t57 + t54 * t56) * qJD(1);
t45 = (-pkin(2) * t57 - pkin(1)) * qJD(1);
t24 = t38 * pkin(3) - t40 * pkin(9) + t45;
t43 = qJD(2) * pkin(2) + t67 * t64;
t44 = t67 * t63;
t29 = t53 * t43 - t56 * t44;
t48 = qJD(2) + qJD(3);
t27 = t48 * pkin(9) + t29;
t52 = sin(qJ(4));
t55 = cos(qJ(4));
t13 = t52 * t24 + t55 * t27;
t31 = t52 * t40 - t55 * t48;
t11 = -t31 * pkin(10) + t13;
t51 = sin(qJ(5));
t66 = cos(qJ(5));
t12 = t55 * t24 - t52 * t27;
t33 = t55 * t40 + t52 * t48;
t37 = qJD(4) + t38;
t7 = t37 * pkin(4) - t33 * pkin(10) + t12;
t4 = t66 * t11 + t51 * t7;
t65 = t57 * t58;
t62 = qJD(1) * qJD(2);
t3 = -t51 * t11 + t66 * t7;
t61 = t54 * t62;
t60 = t57 * t62;
t28 = t56 * t43 + t53 * t44;
t26 = -t48 * pkin(3) - t28;
t16 = t31 * pkin(4) + t26;
t50 = t57 ^ 2;
t49 = t54 ^ 2;
t35 = qJD(5) + t37;
t34 = t35 ^ 2 / 0.2e1;
t21 = -t51 * t31 + t66 * t33;
t19 = t66 * t31 + t51 * t33;
t18 = t21 ^ 2 / 0.2e1;
t17 = t19 ^ 2 / 0.2e1;
t15 = t21 * t35;
t14 = t19 * t35;
t10 = t21 * t19;
t8 = t19 * pkin(5) + qJD(6) + t16;
t2 = -t19 * qJ(6) + t4;
t1 = t35 * pkin(5) - t21 * qJ(6) + t3;
t5 = [0, 0, 0, 0, 0, t68, 0, 0, 0, 0, t49 * t68, t54 * t65, t61, t50 * t68, t60, qJD(2) ^ 2 / 0.2e1, pkin(1) * t65 - pkin(7) * t61, -t58 * pkin(1) * t54 - pkin(7) * t60 (t49 + t50) * t58 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t50 / 0.2e1 + t49 / 0.2e1) * pkin(7) ^ 2) * t58, t40 ^ 2 / 0.2e1, -t40 * t38, t40 * t48, t38 ^ 2 / 0.2e1, -t38 * t48, t48 ^ 2 / 0.2e1, t28 * t48 + t45 * t38, -t29 * t48 + t45 * t40, -t28 * t40 - t29 * t38, t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t45 ^ 2 / 0.2e1, t33 ^ 2 / 0.2e1, -t33 * t31, t33 * t37, t31 ^ 2 / 0.2e1, -t31 * t37, t37 ^ 2 / 0.2e1, t12 * t37 + t26 * t31, -t13 * t37 + t26 * t33, -t12 * t33 - t13 * t31, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1, t18, -t10, t15, t17, -t14, t34, t16 * t19 + t3 * t35, t16 * t21 - t4 * t35, -t4 * t19 - t3 * t21, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1, t18, -t10, t15, t17, -t14, t34, t1 * t35 + t8 * t19, -t2 * t35 + t8 * t21, -t1 * t21 - t2 * t19, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1;];
T_reg  = t5;
