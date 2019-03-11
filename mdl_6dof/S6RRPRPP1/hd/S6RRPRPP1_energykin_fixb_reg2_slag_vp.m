% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPP1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP1_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:47:50
% EndTime: 2019-03-09 09:47:50
% DurationCPUTime: 0.19s
% Computational Cost: add. (821->62), mult. (1986->138), div. (0->0), fcn. (1417->8), ass. (0->52)
t51 = qJD(1) ^ 2;
t66 = t51 / 0.2e1;
t46 = sin(pkin(10));
t59 = cos(pkin(10));
t47 = sin(pkin(9));
t50 = cos(qJ(2));
t57 = qJD(1) * t50;
t49 = sin(qJ(2));
t58 = qJD(1) * t49;
t60 = cos(pkin(9));
t33 = t47 * t58 - t60 * t57;
t35 = (t47 * t50 + t60 * t49) * qJD(1);
t40 = qJD(3) + (-pkin(2) * t50 - pkin(1)) * qJD(1);
t20 = t33 * pkin(3) - t35 * pkin(8) + t40;
t61 = pkin(7) + qJ(3);
t38 = qJD(2) * pkin(2) - t61 * t58;
t39 = t61 * t57;
t25 = t47 * t38 + t60 * t39;
t23 = qJD(2) * pkin(8) + t25;
t48 = sin(qJ(4));
t65 = cos(qJ(4));
t10 = t65 * t20 - t48 * t23;
t29 = t48 * qJD(2) + t65 * t35;
t32 = qJD(4) + t33;
t7 = t32 * pkin(4) - t29 * qJ(5) + t10;
t11 = t48 * t20 + t65 * t23;
t27 = -t65 * qJD(2) + t48 * t35;
t9 = -t27 * qJ(5) + t11;
t4 = t46 * t7 + t59 * t9;
t15 = t59 * t27 + t46 * t29;
t17 = -t46 * t27 + t59 * t29;
t64 = t17 * t15;
t63 = t32 * t15;
t62 = t50 * t51;
t56 = t15 ^ 2 / 0.2e1;
t55 = qJD(1) * qJD(2);
t54 = t49 * t55;
t53 = t50 * t55;
t24 = t60 * t38 - t47 * t39;
t3 = -t46 * t9 + t59 * t7;
t22 = -qJD(2) * pkin(3) - t24;
t13 = t27 * pkin(4) + qJD(5) + t22;
t45 = t50 ^ 2;
t44 = t49 ^ 2;
t43 = qJD(2) ^ 2 / 0.2e1;
t30 = t32 ^ 2 / 0.2e1;
t14 = t17 ^ 2 / 0.2e1;
t12 = t17 * t32;
t5 = t15 * pkin(5) - t17 * qJ(6) + t13;
t2 = t32 * qJ(6) + t4;
t1 = -t32 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t66, 0, 0, 0, 0, t44 * t66, t49 * t62, t54, t45 * t66, t53, t43, pkin(1) * t62 - pkin(7) * t54, -t51 * pkin(1) * t49 - pkin(7) * t53 (t44 + t45) * t51 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t45 / 0.2e1 + t44 / 0.2e1) * pkin(7) ^ 2) * t51, t35 ^ 2 / 0.2e1, -t35 * t33, t35 * qJD(2), t33 ^ 2 / 0.2e1, -t33 * qJD(2), t43, t24 * qJD(2) + t40 * t33, -t25 * qJD(2) + t40 * t35, -t24 * t35 - t25 * t33, t25 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1 + t40 ^ 2 / 0.2e1, t29 ^ 2 / 0.2e1, -t29 * t27, t29 * t32, t27 ^ 2 / 0.2e1, -t27 * t32, t30, t10 * t32 + t22 * t27, -t11 * t32 + t22 * t29, -t10 * t29 - t11 * t27, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t14, -t64, t12, t56, -t63, t30, t13 * t15 + t3 * t32, t13 * t17 - t4 * t32, -t4 * t15 - t3 * t17, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1, t14, t12, t64, t30, t63, t56, -t1 * t32 + t5 * t15, t1 * t17 - t2 * t15, -t5 * t17 + t2 * t32, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
