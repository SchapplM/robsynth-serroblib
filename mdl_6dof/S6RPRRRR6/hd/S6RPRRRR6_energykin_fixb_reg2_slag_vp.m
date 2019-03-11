% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRRR6
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
% Datum: 2019-03-09 07:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRR6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR6_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:14:32
% EndTime: 2019-03-09 07:14:32
% DurationCPUTime: 0.22s
% Computational Cost: add. (1102->67), mult. (2673->155), div. (0->0), fcn. (2023->10), ass. (0->50)
t59 = qJD(1) ^ 2;
t67 = t59 / 0.2e1;
t57 = sin(qJ(3));
t58 = cos(qJ(3));
t53 = cos(pkin(11));
t61 = qJD(1) * t53;
t52 = sin(pkin(11));
t62 = qJD(1) * t52;
t41 = t57 * t62 - t58 * t61;
t43 = (t52 * t58 + t53 * t57) * qJD(1);
t46 = qJD(2) + (-pkin(2) * t53 - pkin(1)) * qJD(1);
t25 = t41 * pkin(3) - t43 * pkin(8) + t46;
t63 = pkin(7) + qJ(2);
t44 = t63 * t62;
t45 = t63 * t61;
t30 = -t57 * t44 + t58 * t45;
t28 = qJD(3) * pkin(8) + t30;
t56 = sin(qJ(4));
t66 = cos(qJ(4));
t17 = t56 * t25 + t66 * t28;
t32 = -t66 * qJD(3) + t56 * t43;
t15 = -t32 * pkin(9) + t17;
t55 = sin(qJ(5));
t65 = cos(qJ(5));
t16 = t66 * t25 - t56 * t28;
t34 = t56 * qJD(3) + t66 * t43;
t37 = qJD(4) + t41;
t9 = t37 * pkin(4) - t34 * pkin(9) + t16;
t6 = t65 * t15 + t55 * t9;
t64 = cos(qJ(6));
t5 = -t55 * t15 + t65 * t9;
t29 = -t58 * t44 - t57 * t45;
t27 = -qJD(3) * pkin(3) - t29;
t36 = qJD(5) + t37;
t18 = t32 * pkin(4) + t27;
t54 = sin(qJ(6));
t51 = t53 ^ 2;
t50 = t52 ^ 2;
t48 = -qJD(1) * pkin(1) + qJD(2);
t35 = qJD(6) + t36;
t22 = -t55 * t32 + t65 * t34;
t20 = t65 * t32 + t55 * t34;
t13 = t20 * pkin(5) + t18;
t12 = -t54 * t20 + t64 * t22;
t10 = t64 * t20 + t54 * t22;
t4 = -t20 * pkin(10) + t6;
t3 = t36 * pkin(5) - t22 * pkin(10) + t5;
t2 = t54 * t3 + t64 * t4;
t1 = t64 * t3 - t54 * t4;
t7 = [0, 0, 0, 0, 0, t67, 0, 0, 0, 0, t50 * t67, t52 * t59 * t53, 0, t51 * t67, 0, 0, -t48 * t61, t48 * t62 (t50 + t51) * t59 * qJ(2), t48 ^ 2 / 0.2e1 + (t51 / 0.2e1 + t50 / 0.2e1) * qJ(2) ^ 2 * t59, t43 ^ 2 / 0.2e1, -t43 * t41, t43 * qJD(3), t41 ^ 2 / 0.2e1, -t41 * qJD(3), qJD(3) ^ 2 / 0.2e1, t29 * qJD(3) + t46 * t41, -t30 * qJD(3) + t46 * t43, -t29 * t43 - t30 * t41, t30 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1 + t46 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, t34 * t37, t32 ^ 2 / 0.2e1, -t32 * t37, t37 ^ 2 / 0.2e1, t16 * t37 + t27 * t32, -t17 * t37 + t27 * t34, -t16 * t34 - t17 * t32, t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * t36, t20 ^ 2 / 0.2e1, -t20 * t36, t36 ^ 2 / 0.2e1, t18 * t20 + t5 * t36, t18 * t22 - t6 * t36, -t6 * t20 - t5 * t22, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t12 ^ 2 / 0.2e1, -t12 * t10, t12 * t35, t10 ^ 2 / 0.2e1, -t10 * t35, t35 ^ 2 / 0.2e1, t1 * t35 + t13 * t10, t13 * t12 - t2 * t35, -t1 * t12 - t2 * t10, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1;];
T_reg  = t7;
