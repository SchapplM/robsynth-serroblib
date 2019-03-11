% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPRR11_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR11_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:42:42
% EndTime: 2019-03-09 09:42:43
% DurationCPUTime: 0.24s
% Computational Cost: add. (944->71), mult. (2265->154), div. (0->0), fcn. (1653->10), ass. (0->58)
t56 = sin(qJ(2));
t51 = sin(pkin(6));
t68 = qJD(1) * t51;
t63 = t56 * t68;
t44 = pkin(8) * t63;
t53 = cos(pkin(6));
t67 = t53 * qJD(1);
t48 = qJD(2) + t67;
t57 = cos(qJ(2));
t69 = -pkin(2) - qJ(4);
t24 = qJD(3) + t44 + t69 * t48 + (-pkin(1) * t53 * t57 + pkin(3) * t51 * t56) * qJD(1);
t61 = -qJ(3) * t56 - pkin(1);
t27 = (t69 * t57 + t61) * t68;
t50 = sin(pkin(11));
t52 = cos(pkin(11));
t13 = t50 * t24 + t52 * t27;
t64 = t57 * t68;
t32 = t50 * t48 + t52 * t64;
t11 = -t32 * pkin(9) + t13;
t55 = sin(qJ(5));
t72 = cos(qJ(5));
t12 = t52 * t24 - t50 * t27;
t34 = t52 * t48 - t50 * t64;
t9 = pkin(4) * t63 - t34 * pkin(9) + t12;
t6 = t72 * t11 + t55 * t9;
t71 = cos(qJ(6));
t58 = qJD(1) ^ 2;
t70 = t51 ^ 2 * t58;
t65 = pkin(1) * t67;
t36 = pkin(8) * t64 + t56 * t65;
t66 = t57 * t70;
t29 = -t48 * qJ(3) - t36;
t62 = t70 / 0.2e1;
t19 = t72 * t32 + t55 * t34;
t60 = t48 * t63;
t59 = t48 * t64;
t25 = pkin(3) * t64 + qJD(4) - t29;
t35 = t57 * t65 - t44;
t5 = -t55 * t11 + t72 * t9;
t17 = t32 * pkin(4) + t25;
t54 = sin(qJ(6));
t42 = t57 ^ 2 * t62;
t41 = t56 ^ 2 * t62;
t40 = t48 ^ 2 / 0.2e1;
t39 = qJD(5) + t63;
t38 = t56 * t66;
t31 = (-pkin(2) * t57 + t61) * t68;
t28 = -t48 * pkin(2) + qJD(3) - t35;
t21 = -t55 * t32 + t72 * t34;
t18 = qJD(6) + t19;
t16 = t71 * t21 + t54 * t39;
t14 = t54 * t21 - t71 * t39;
t7 = t19 * pkin(5) - t21 * pkin(10) + t17;
t4 = t39 * pkin(10) + t6;
t3 = -t39 * pkin(5) - t5;
t2 = t71 * t4 + t54 * t7;
t1 = -t54 * t4 + t71 * t7;
t8 = [0, 0, 0, 0, 0, t58 / 0.2e1, 0, 0, 0, 0, t41, t38, t60, t42, t59, t40, pkin(1) * t66 + t35 * t48, -pkin(1) * t56 * t70 - t36 * t48 (-t35 * t56 + t36 * t57) * t68, t36 ^ 2 / 0.2e1 + t35 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t62, t40, -t60, -t59, t41, t38, t42 (t28 * t56 - t29 * t57) * t68, t28 * t48 + t31 * t64, -t29 * t48 - t31 * t63, t31 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, t34 * t63, t32 ^ 2 / 0.2e1, -t32 * t63, t41, t12 * t63 + t25 * t32, -t13 * t63 + t25 * t34, -t12 * t34 - t13 * t32, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1, t21 ^ 2 / 0.2e1, -t21 * t19, t21 * t39, t19 ^ 2 / 0.2e1, -t19 * t39, t39 ^ 2 / 0.2e1, t17 * t19 + t5 * t39, t17 * t21 - t6 * t39, -t6 * t19 - t5 * t21, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t16 ^ 2 / 0.2e1, -t16 * t14, t16 * t18, t14 ^ 2 / 0.2e1, -t14 * t18, t18 ^ 2 / 0.2e1, t1 * t18 + t3 * t14, t3 * t16 - t2 * t18, -t1 * t16 - t2 * t14, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
