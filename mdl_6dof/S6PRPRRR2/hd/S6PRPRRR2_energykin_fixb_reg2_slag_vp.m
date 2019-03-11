% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRPRRR2
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
% Datum: 2019-03-08 20:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRRR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR2_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:29:54
% EndTime: 2019-03-08 20:29:54
% DurationCPUTime: 0.14s
% Computational Cost: add. (424->52), mult. (974->132), div. (0->0), fcn. (695->12), ass. (0->46)
t52 = qJD(2) ^ 2;
t41 = t52 / 0.2e1;
t51 = cos(qJ(2));
t43 = sin(pkin(6));
t59 = qJD(1) * t43;
t30 = qJD(2) * pkin(2) + t51 * t59;
t42 = sin(pkin(12));
t44 = cos(pkin(12));
t49 = sin(qJ(2));
t55 = t49 * t59;
t24 = t42 * t30 + t44 * t55;
t22 = qJD(2) * pkin(8) + t24;
t45 = cos(pkin(6));
t36 = t45 * qJD(1) + qJD(3);
t48 = sin(qJ(4));
t50 = cos(qJ(4));
t12 = t50 * t22 + t48 * t36;
t10 = qJD(4) * pkin(9) + t12;
t23 = t44 * t30 - t42 * t55;
t15 = (-pkin(4) * t50 - pkin(9) * t48 - pkin(3)) * qJD(2) - t23;
t47 = sin(qJ(5));
t61 = cos(qJ(5));
t6 = t61 * t10 + t47 * t15;
t60 = cos(qJ(6));
t58 = qJD(2) * t48;
t57 = t50 * qJD(2);
t56 = qJD(2) * qJD(4);
t54 = qJD(2) * t59;
t5 = -t47 * t10 + t61 * t15;
t11 = -t48 * t22 + t50 * t36;
t37 = -qJD(5) + t57;
t9 = -qJD(4) * pkin(4) - t11;
t53 = qJD(1) ^ 2;
t46 = sin(qJ(6));
t34 = -qJD(6) + t37;
t29 = t47 * qJD(4) + t61 * t58;
t27 = -t61 * qJD(4) + t47 * t58;
t21 = -qJD(2) * pkin(3) - t23;
t18 = -t46 * t27 + t60 * t29;
t16 = t60 * t27 + t46 * t29;
t7 = t27 * pkin(5) + t9;
t4 = -t27 * pkin(10) + t6;
t3 = -t37 * pkin(5) - t29 * pkin(10) + t5;
t2 = t46 * t3 + t60 * t4;
t1 = t60 * t3 - t46 * t4;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t53 / 0.2e1, 0, 0, 0, 0, 0, t41, t51 * t54, -t49 * t54, 0 (t45 ^ 2 / 0.2e1 + (t49 ^ 2 / 0.2e1 + t51 ^ 2 / 0.2e1) * t43 ^ 2) * t53, 0, 0, 0, 0, 0, t41, t23 * qJD(2), -t24 * qJD(2), 0, t24 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + t36 ^ 2 / 0.2e1, t48 ^ 2 * t41, t48 * t52 * t50, t48 * t56, t50 ^ 2 * t41, t50 * t56, qJD(4) ^ 2 / 0.2e1, t11 * qJD(4) - t21 * t57, -t12 * qJD(4) + t21 * t58 (-t11 * t48 + t12 * t50) * qJD(2), t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t29 ^ 2 / 0.2e1, -t29 * t27, -t29 * t37, t27 ^ 2 / 0.2e1, t27 * t37, t37 ^ 2 / 0.2e1, t9 * t27 - t5 * t37, t9 * t29 + t6 * t37, -t6 * t27 - t5 * t29, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, -t18 * t34, t16 ^ 2 / 0.2e1, t16 * t34, t34 ^ 2 / 0.2e1, -t1 * t34 + t7 * t16, t7 * t18 + t2 * t34, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
