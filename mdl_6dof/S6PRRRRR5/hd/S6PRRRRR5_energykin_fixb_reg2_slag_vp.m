% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRRR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR5_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_energykin_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:10:00
% EndTime: 2019-03-09 01:10:00
% DurationCPUTime: 0.26s
% Computational Cost: add. (983->64), mult. (2342->160), div. (0->0), fcn. (1897->14), ass. (0->57)
t62 = cos(qJ(2));
t52 = sin(pkin(6));
t73 = qJD(1) * t52;
t42 = qJD(2) * pkin(2) + t62 * t73;
t51 = sin(pkin(7));
t53 = cos(pkin(7));
t54 = cos(pkin(6));
t72 = qJD(1) * t54;
t78 = t42 * t53 + t51 * t72;
t59 = sin(qJ(2));
t71 = qJD(2) * t51;
t40 = pkin(9) * t71 + t59 * t73;
t58 = sin(qJ(3));
t61 = cos(qJ(3));
t23 = -t58 * t40 + t78 * t61;
t24 = t61 * t40 + t78 * t58;
t48 = t53 * qJD(2) + qJD(3);
t22 = t48 * pkin(10) + t24;
t47 = t53 * t72;
t27 = t47 + (-t42 + (-pkin(3) * t61 - pkin(10) * t58) * qJD(2)) * t51;
t57 = sin(qJ(4));
t60 = cos(qJ(4));
t14 = t60 * t22 + t57 * t27;
t68 = t61 * t71;
t45 = -qJD(4) + t68;
t10 = -t45 * pkin(11) + t14;
t21 = -t48 * pkin(3) - t23;
t69 = t58 * t71;
t34 = -t60 * t48 + t57 * t69;
t36 = t57 * t48 + t60 * t69;
t15 = t34 * pkin(4) - t36 * pkin(11) + t21;
t56 = sin(qJ(5));
t77 = cos(qJ(5));
t6 = t77 * t10 + t56 * t15;
t76 = cos(qJ(6));
t63 = qJD(2) ^ 2;
t74 = t51 ^ 2 * t63;
t67 = t74 / 0.2e1;
t66 = qJD(2) * t73;
t5 = -t56 * t10 + t77 * t15;
t13 = -t57 * t22 + t60 * t27;
t9 = t45 * pkin(4) - t13;
t33 = qJD(5) + t34;
t64 = qJD(1) ^ 2;
t55 = sin(qJ(6));
t32 = -t51 * t42 + t47;
t31 = qJD(6) + t33;
t30 = t77 * t36 - t56 * t45;
t28 = t56 * t36 + t77 * t45;
t18 = -t55 * t28 + t76 * t30;
t16 = t76 * t28 + t55 * t30;
t7 = t28 * pkin(5) + t9;
t4 = -t28 * pkin(12) + t6;
t3 = t33 * pkin(5) - t30 * pkin(12) + t5;
t2 = t55 * t3 + t76 * t4;
t1 = t76 * t3 - t55 * t4;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t64 / 0.2e1, 0, 0, 0, 0, 0, t63 / 0.2e1, t62 * t66, -t59 * t66, 0 (t54 ^ 2 / 0.2e1 + (t59 ^ 2 / 0.2e1 + t62 ^ 2 / 0.2e1) * t52 ^ 2) * t64, t58 ^ 2 * t67, t58 * t61 * t74, t48 * t69, t61 ^ 2 * t67, t48 * t68, t48 ^ 2 / 0.2e1, t23 * t48 - t32 * t68, -t24 * t48 + t32 * t69 (-t23 * t58 + t24 * t61) * t71, t24 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1, t36 ^ 2 / 0.2e1, -t36 * t34, -t36 * t45, t34 ^ 2 / 0.2e1, t34 * t45, t45 ^ 2 / 0.2e1, -t13 * t45 + t21 * t34, t14 * t45 + t21 * t36, -t13 * t36 - t14 * t34, t14 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t30 ^ 2 / 0.2e1, -t30 * t28, t30 * t33, t28 ^ 2 / 0.2e1, -t28 * t33, t33 ^ 2 / 0.2e1, t9 * t28 + t5 * t33, t9 * t30 - t6 * t33, -t6 * t28 - t5 * t30, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t31, t16 ^ 2 / 0.2e1, -t16 * t31, t31 ^ 2 / 0.2e1, t1 * t31 + t7 * t16, t7 * t18 - t2 * t31, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
