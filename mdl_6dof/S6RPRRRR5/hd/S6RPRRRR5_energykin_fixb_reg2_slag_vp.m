% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR5_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:09:48
% EndTime: 2019-03-09 07:09:48
% DurationCPUTime: 0.22s
% Computational Cost: add. (1081->67), mult. (2702->155), div. (0->0), fcn. (2082->10), ass. (0->50)
t59 = qJD(1) ^ 2;
t67 = t59 / 0.2e1;
t52 = sin(pkin(11));
t62 = qJD(1) * t52;
t63 = pkin(7) + qJ(2);
t43 = t63 * t62;
t53 = cos(pkin(11));
t61 = qJD(1) * t53;
t44 = t63 * t61;
t57 = sin(qJ(3));
t66 = cos(qJ(3));
t33 = -t66 * t43 - t57 * t44;
t42 = (t66 * t52 + t53 * t57) * qJD(1);
t22 = qJD(3) * pkin(3) - t42 * pkin(8) + t33;
t34 = -t57 * t43 + t66 * t44;
t40 = t57 * t62 - t66 * t61;
t24 = -t40 * pkin(8) + t34;
t56 = sin(qJ(4));
t58 = cos(qJ(4));
t14 = t56 * t22 + t58 * t24;
t51 = qJD(3) + qJD(4);
t10 = t51 * pkin(9) + t14;
t30 = t58 * t40 + t56 * t42;
t32 = -t56 * t40 + t58 * t42;
t45 = qJD(2) + (-pkin(2) * t53 - pkin(1)) * qJD(1);
t35 = t40 * pkin(3) + t45;
t15 = t30 * pkin(4) - t32 * pkin(9) + t35;
t55 = sin(qJ(5));
t65 = cos(qJ(5));
t6 = t65 * t10 + t55 * t15;
t64 = cos(qJ(6));
t5 = -t55 * t10 + t65 * t15;
t13 = t58 * t22 - t56 * t24;
t29 = qJD(5) + t30;
t9 = -t51 * pkin(4) - t13;
t54 = sin(qJ(6));
t50 = t53 ^ 2;
t49 = t52 ^ 2;
t48 = -qJD(1) * pkin(1) + qJD(2);
t28 = qJD(6) + t29;
t27 = t65 * t32 + t55 * t51;
t25 = t55 * t32 - t65 * t51;
t18 = -t54 * t25 + t64 * t27;
t16 = t64 * t25 + t54 * t27;
t7 = t25 * pkin(5) + t9;
t4 = -t25 * pkin(10) + t6;
t3 = t29 * pkin(5) - t27 * pkin(10) + t5;
t2 = t54 * t3 + t64 * t4;
t1 = t64 * t3 - t54 * t4;
t8 = [0, 0, 0, 0, 0, t67, 0, 0, 0, 0, t49 * t67, t52 * t59 * t53, 0, t50 * t67, 0, 0, -t48 * t61, t48 * t62 (t49 + t50) * t59 * qJ(2), t48 ^ 2 / 0.2e1 + (t50 / 0.2e1 + t49 / 0.2e1) * qJ(2) ^ 2 * t59, t42 ^ 2 / 0.2e1, -t42 * t40, t42 * qJD(3), t40 ^ 2 / 0.2e1, -t40 * qJD(3), qJD(3) ^ 2 / 0.2e1, t33 * qJD(3) + t45 * t40, -t34 * qJD(3) + t45 * t42, -t33 * t42 - t34 * t40, t34 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1 + t45 ^ 2 / 0.2e1, t32 ^ 2 / 0.2e1, -t32 * t30, t32 * t51, t30 ^ 2 / 0.2e1, -t30 * t51, t51 ^ 2 / 0.2e1, t13 * t51 + t35 * t30, -t14 * t51 + t35 * t32, -t13 * t32 - t14 * t30, t14 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1 + t35 ^ 2 / 0.2e1, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t29, t25 ^ 2 / 0.2e1, -t25 * t29, t29 ^ 2 / 0.2e1, t9 * t25 + t5 * t29, t9 * t27 - t6 * t29, -t6 * t25 - t5 * t27, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t28, t16 ^ 2 / 0.2e1, -t16 * t28, t28 ^ 2 / 0.2e1, t1 * t28 + t7 * t16, t7 * t18 - t2 * t28, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
