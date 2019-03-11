% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRRP4
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
% Datum: 2019-03-10 01:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRP4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP4_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:15:39
% EndTime: 2019-03-10 01:15:39
% DurationCPUTime: 0.21s
% Computational Cost: add. (934->62), mult. (1986->143), div. (0->0), fcn. (1417->8), ass. (0->53)
t53 = qJD(1) ^ 2;
t67 = t53 / 0.2e1;
t66 = -pkin(8) - pkin(7);
t47 = sin(qJ(5));
t64 = cos(qJ(5));
t49 = sin(qJ(3));
t51 = cos(qJ(3));
t52 = cos(qJ(2));
t59 = qJD(1) * t52;
t50 = sin(qJ(2));
t60 = qJD(1) * t50;
t34 = t49 * t60 - t51 * t59;
t36 = (t49 * t52 + t50 * t51) * qJD(1);
t41 = (-pkin(2) * t52 - pkin(1)) * qJD(1);
t20 = t34 * pkin(3) - t36 * pkin(9) + t41;
t39 = qJD(2) * pkin(2) + t66 * t60;
t40 = t66 * t59;
t25 = t49 * t39 - t51 * t40;
t44 = qJD(2) + qJD(3);
t23 = t44 * pkin(9) + t25;
t48 = sin(qJ(4));
t65 = cos(qJ(4));
t10 = t65 * t20 - t48 * t23;
t29 = t65 * t36 + t48 * t44;
t33 = qJD(4) + t34;
t7 = t33 * pkin(4) - t29 * pkin(10) + t10;
t11 = t48 * t20 + t65 * t23;
t27 = t48 * t36 - t65 * t44;
t9 = -t27 * pkin(10) + t11;
t4 = t47 * t7 + t64 * t9;
t15 = t64 * t27 + t47 * t29;
t17 = -t47 * t27 + t64 * t29;
t63 = t15 * t17;
t31 = qJD(5) + t33;
t62 = t15 * t31;
t61 = t52 * t53;
t58 = t15 ^ 2 / 0.2e1;
t57 = qJD(1) * qJD(2);
t56 = t50 * t57;
t55 = t52 * t57;
t24 = t51 * t39 + t49 * t40;
t22 = -t44 * pkin(3) - t24;
t3 = -t47 * t9 + t64 * t7;
t13 = t27 * pkin(4) + t22;
t46 = t52 ^ 2;
t45 = t50 ^ 2;
t30 = t31 ^ 2 / 0.2e1;
t14 = t17 ^ 2 / 0.2e1;
t12 = t17 * t31;
t5 = t15 * pkin(5) - t17 * qJ(6) + t13;
t2 = t31 * qJ(6) + t4;
t1 = -t31 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t67, 0, 0, 0, 0, t45 * t67, t50 * t61, t56, t46 * t67, t55, qJD(2) ^ 2 / 0.2e1, pkin(1) * t61 - pkin(7) * t56, -t53 * pkin(1) * t50 - pkin(7) * t55 (t45 + t46) * t53 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t46 / 0.2e1 + t45 / 0.2e1) * pkin(7) ^ 2) * t53, t36 ^ 2 / 0.2e1, -t36 * t34, t36 * t44, t34 ^ 2 / 0.2e1, -t34 * t44, t44 ^ 2 / 0.2e1, t24 * t44 + t41 * t34, -t25 * t44 + t41 * t36, -t24 * t36 - t25 * t34, t25 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1 + t41 ^ 2 / 0.2e1, t29 ^ 2 / 0.2e1, -t29 * t27, t29 * t33, t27 ^ 2 / 0.2e1, -t27 * t33, t33 ^ 2 / 0.2e1, t10 * t33 + t22 * t27, -t11 * t33 + t22 * t29, -t10 * t29 - t11 * t27, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t14, -t63, t12, t58, -t62, t30, t13 * t15 + t3 * t31, t13 * t17 - t4 * t31, -t4 * t15 - t3 * t17, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1, t14, t12, t63, t30, t62, t58, -t1 * t31 + t5 * t15, t1 * t17 - t2 * t15, -t5 * t17 + t2 * t31, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
