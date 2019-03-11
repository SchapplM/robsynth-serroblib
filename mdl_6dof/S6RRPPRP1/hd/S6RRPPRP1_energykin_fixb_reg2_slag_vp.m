% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPRP1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP1_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:27:46
% EndTime: 2019-03-09 08:27:46
% DurationCPUTime: 0.21s
% Computational Cost: add. (807->62), mult. (1986->137), div. (0->0), fcn. (1417->8), ass. (0->53)
t51 = qJD(1) ^ 2;
t67 = t51 / 0.2e1;
t48 = sin(qJ(5));
t66 = cos(qJ(5));
t47 = sin(pkin(9));
t50 = cos(qJ(2));
t58 = qJD(1) * t50;
t49 = sin(qJ(2));
t59 = qJD(1) * t49;
t61 = cos(pkin(9));
t33 = t47 * t59 - t61 * t58;
t35 = (t47 * t50 + t61 * t49) * qJD(1);
t40 = qJD(3) + (-pkin(2) * t50 - pkin(1)) * qJD(1);
t20 = t33 * pkin(3) - t35 * qJ(4) + t40;
t62 = pkin(7) + qJ(3);
t38 = qJD(2) * pkin(2) - t62 * t59;
t39 = t62 * t58;
t25 = t47 * t38 + t61 * t39;
t23 = qJD(2) * qJ(4) + t25;
t46 = sin(pkin(10));
t60 = cos(pkin(10));
t10 = t60 * t20 - t46 * t23;
t29 = t46 * qJD(2) + t60 * t35;
t7 = t33 * pkin(4) - t29 * pkin(8) + t10;
t11 = t46 * t20 + t60 * t23;
t27 = -t60 * qJD(2) + t46 * t35;
t9 = -t27 * pkin(8) + t11;
t4 = t48 * t7 + t66 * t9;
t15 = t66 * t27 + t48 * t29;
t17 = -t48 * t27 + t66 * t29;
t65 = t17 * t15;
t32 = qJD(5) + t33;
t64 = t32 * t15;
t63 = t50 * t51;
t57 = t15 ^ 2 / 0.2e1;
t56 = t33 ^ 2 / 0.2e1;
t55 = qJD(1) * qJD(2);
t54 = t49 * t55;
t53 = t50 * t55;
t24 = t61 * t38 - t47 * t39;
t3 = -t48 * t9 + t66 * t7;
t22 = -qJD(2) * pkin(3) + qJD(4) - t24;
t13 = t27 * pkin(4) + t22;
t45 = t50 ^ 2;
t44 = t49 ^ 2;
t43 = qJD(2) ^ 2 / 0.2e1;
t30 = t32 ^ 2 / 0.2e1;
t14 = t17 ^ 2 / 0.2e1;
t12 = t17 * t32;
t5 = t15 * pkin(5) - t17 * qJ(6) + t13;
t2 = t32 * qJ(6) + t4;
t1 = -t32 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t67, 0, 0, 0, 0, t44 * t67, t49 * t63, t54, t45 * t67, t53, t43, pkin(1) * t63 - pkin(7) * t54, -t51 * pkin(1) * t49 - pkin(7) * t53 (t44 + t45) * t51 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t45 / 0.2e1 + t44 / 0.2e1) * pkin(7) ^ 2) * t51, t35 ^ 2 / 0.2e1, -t35 * t33, t35 * qJD(2), t56, -t33 * qJD(2), t43, t24 * qJD(2) + t40 * t33, -t25 * qJD(2) + t40 * t35, -t24 * t35 - t25 * t33, t25 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1 + t40 ^ 2 / 0.2e1, t29 ^ 2 / 0.2e1, -t29 * t27, t29 * t33, t27 ^ 2 / 0.2e1, -t27 * t33, t56, t10 * t33 + t22 * t27, -t11 * t33 + t22 * t29, -t10 * t29 - t11 * t27, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t14, -t65, t12, t57, -t64, t30, t13 * t15 + t3 * t32, t13 * t17 - t4 * t32, -t4 * t15 - t3 * t17, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1, t14, t12, t65, t30, t64, t57, -t1 * t32 + t5 * t15, t1 * t17 - t2 * t15, -t5 * t17 + t2 * t32, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
