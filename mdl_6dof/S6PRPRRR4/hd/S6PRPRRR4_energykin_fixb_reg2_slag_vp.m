% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRRR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR4_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:38:55
% EndTime: 2019-03-08 20:38:55
% DurationCPUTime: 0.16s
% Computational Cost: add. (597->60), mult. (1459->146), div. (0->0), fcn. (1112->12), ass. (0->50)
t54 = qJD(2) ^ 2;
t64 = t54 / 0.2e1;
t51 = sin(qJ(2));
t45 = sin(pkin(6));
t61 = qJD(1) * t45;
t38 = qJD(2) * qJ(3) + t51 * t61;
t46 = cos(pkin(12));
t47 = cos(pkin(6));
t60 = qJD(1) * t47;
t40 = t46 * t60;
t44 = sin(pkin(12));
t22 = t40 + (-pkin(8) * qJD(2) - t38) * t44;
t29 = t46 * t38 + t44 * t60;
t58 = qJD(2) * t46;
t23 = pkin(8) * t58 + t29;
t50 = sin(qJ(4));
t52 = cos(qJ(4));
t12 = t50 * t22 + t52 * t23;
t10 = qJD(4) * pkin(9) + t12;
t53 = cos(qJ(2));
t56 = -t53 * t61 + qJD(3);
t31 = (-pkin(3) * t46 - pkin(2)) * qJD(2) + t56;
t59 = qJD(2) * t44;
t33 = t50 * t59 - t52 * t58;
t35 = (t44 * t52 + t46 * t50) * qJD(2);
t18 = t33 * pkin(4) - t35 * pkin(9) + t31;
t49 = sin(qJ(5));
t63 = cos(qJ(5));
t6 = t63 * t10 + t49 * t18;
t62 = cos(qJ(6));
t57 = qJD(2) * t61;
t5 = -t49 * t10 + t63 * t18;
t11 = t52 * t22 - t50 * t23;
t32 = qJD(5) + t33;
t9 = -qJD(4) * pkin(4) - t11;
t55 = qJD(1) ^ 2;
t48 = sin(qJ(6));
t37 = -qJD(2) * pkin(2) + t56;
t30 = qJD(6) + t32;
t28 = -t44 * t38 + t40;
t27 = t49 * qJD(4) + t63 * t35;
t25 = -t63 * qJD(4) + t49 * t35;
t17 = -t48 * t25 + t62 * t27;
t15 = t62 * t25 + t48 * t27;
t7 = t25 * pkin(5) + t9;
t4 = -t25 * pkin(10) + t6;
t3 = t32 * pkin(5) - t27 * pkin(10) + t5;
t2 = t48 * t3 + t62 * t4;
t1 = t62 * t3 - t48 * t4;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t55 / 0.2e1, 0, 0, 0, 0, 0, t64, t53 * t57, -t51 * t57, 0 (t47 ^ 2 / 0.2e1 + (t51 ^ 2 / 0.2e1 + t53 ^ 2 / 0.2e1) * t45 ^ 2) * t55, t44 ^ 2 * t64, t44 * t54 * t46, 0, t46 ^ 2 * t64, 0, 0, -t37 * t58, t37 * t59 (-t28 * t44 + t29 * t46) * qJD(2), t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t37 ^ 2 / 0.2e1, t35 ^ 2 / 0.2e1, -t35 * t33, t35 * qJD(4), t33 ^ 2 / 0.2e1, -t33 * qJD(4), qJD(4) ^ 2 / 0.2e1, t11 * qJD(4) + t31 * t33, -t12 * qJD(4) + t31 * t35, -t11 * t35 - t12 * t33, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t32, t25 ^ 2 / 0.2e1, -t25 * t32, t32 ^ 2 / 0.2e1, t9 * t25 + t5 * t32, t9 * t27 - t6 * t32, -t6 * t25 - t5 * t27, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t17 ^ 2 / 0.2e1, -t17 * t15, t17 * t30, t15 ^ 2 / 0.2e1, -t15 * t30, t30 ^ 2 / 0.2e1, t1 * t30 + t7 * t15, t7 * t17 - t2 * t30, -t1 * t17 - t2 * t15, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
