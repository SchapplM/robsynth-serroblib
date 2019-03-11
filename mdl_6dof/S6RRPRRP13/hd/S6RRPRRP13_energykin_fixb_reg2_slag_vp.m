% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRP13_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP13_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:01:31
% EndTime: 2019-03-09 13:01:31
% DurationCPUTime: 0.21s
% Computational Cost: add. (709->67), mult. (1690->136), div. (0->0), fcn. (1197->8), ass. (0->58)
t71 = -pkin(2) - pkin(9);
t55 = sin(qJ(2));
t57 = cos(qJ(2));
t51 = sin(pkin(6));
t68 = qJD(1) * t51;
t64 = t57 * t68;
t52 = cos(pkin(6));
t67 = t52 * qJD(1);
t65 = pkin(1) * t67;
t37 = pkin(8) * t64 + t55 * t65;
t49 = qJD(2) + t67;
t29 = -t49 * qJ(3) - t37;
t26 = pkin(3) * t64 - t29;
t54 = sin(qJ(4));
t56 = cos(qJ(4));
t33 = t54 * t49 + t56 * t64;
t35 = t56 * t49 - t54 * t64;
t13 = t33 * pkin(4) - t35 * pkin(10) + t26;
t53 = sin(qJ(5));
t70 = cos(qJ(5));
t63 = t55 * t68;
t45 = pkin(8) * t63;
t18 = qJD(3) + t45 + t71 * t49 + (-pkin(1) * t52 * t57 + pkin(3) * t51 * t55) * qJD(1);
t61 = -qJ(3) * t55 - pkin(1);
t27 = (t71 * t57 + t61) * t68;
t10 = t54 * t18 + t56 * t27;
t40 = qJD(4) + t63;
t8 = t40 * pkin(10) + t10;
t4 = t53 * t13 + t70 * t8;
t58 = qJD(1) ^ 2;
t69 = t51 ^ 2 * t58;
t66 = t57 * t69;
t62 = t69 / 0.2e1;
t3 = t70 * t13 - t53 * t8;
t9 = t56 * t18 - t54 * t27;
t60 = t49 * t63;
t59 = t49 * t64;
t36 = t57 * t65 - t45;
t7 = -t40 * pkin(4) - t9;
t43 = t57 ^ 2 * t62;
t42 = t55 ^ 2 * t62;
t41 = t49 ^ 2 / 0.2e1;
t39 = t55 * t66;
t32 = qJD(5) + t33;
t31 = (-pkin(2) * t57 + t61) * t68;
t30 = t32 ^ 2 / 0.2e1;
t28 = -t49 * pkin(2) + qJD(3) - t36;
t23 = t70 * t35 + t53 * t40;
t21 = t53 * t35 - t70 * t40;
t20 = t23 ^ 2 / 0.2e1;
t19 = t21 ^ 2 / 0.2e1;
t16 = t23 * t32;
t15 = t21 * t32;
t14 = t23 * t21;
t5 = t21 * pkin(5) + qJD(6) + t7;
t2 = -t21 * qJ(6) + t4;
t1 = t32 * pkin(5) - t23 * qJ(6) + t3;
t6 = [0, 0, 0, 0, 0, t58 / 0.2e1, 0, 0, 0, 0, t42, t39, t60, t43, t59, t41, pkin(1) * t66 + t36 * t49, -pkin(1) * t55 * t69 - t37 * t49 (-t36 * t55 + t37 * t57) * t68, t37 ^ 2 / 0.2e1 + t36 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t62, t41, -t60, -t59, t42, t39, t43 (t28 * t55 - t29 * t57) * t68, t28 * t49 + t31 * t64, -t29 * t49 - t31 * t63, t31 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1, t35 ^ 2 / 0.2e1, -t35 * t33, t35 * t40, t33 ^ 2 / 0.2e1, -t33 * t40, t40 ^ 2 / 0.2e1, t26 * t33 + t9 * t40, -t10 * t40 + t26 * t35, -t10 * t33 - t9 * t35, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1, t20, -t14, t16, t19, -t15, t30, t7 * t21 + t3 * t32, t7 * t23 - t4 * t32, -t4 * t21 - t3 * t23, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t20, -t14, t16, t19, -t15, t30, t1 * t32 + t5 * t21, -t2 * t32 + t5 * t23, -t1 * t23 - t2 * t21, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t6;
