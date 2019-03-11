% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRP4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP4_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:55:38
% EndTime: 2019-03-09 11:55:38
% DurationCPUTime: 0.21s
% Computational Cost: add. (831->62), mult. (1986->140), div. (0->0), fcn. (1417->8), ass. (0->53)
t53 = qJD(1) ^ 2;
t67 = t53 / 0.2e1;
t49 = sin(qJ(5));
t65 = cos(qJ(5));
t47 = sin(pkin(10));
t48 = cos(pkin(10));
t52 = cos(qJ(2));
t59 = qJD(1) * t52;
t51 = sin(qJ(2));
t60 = qJD(1) * t51;
t34 = t47 * t60 - t48 * t59;
t36 = (t47 * t52 + t48 * t51) * qJD(1);
t41 = qJD(3) + (-pkin(2) * t52 - pkin(1)) * qJD(1);
t20 = t34 * pkin(3) - t36 * pkin(8) + t41;
t61 = pkin(7) + qJ(3);
t39 = qJD(2) * pkin(2) - t61 * t60;
t40 = t61 * t59;
t25 = t47 * t39 + t48 * t40;
t23 = qJD(2) * pkin(8) + t25;
t50 = sin(qJ(4));
t66 = cos(qJ(4));
t10 = t66 * t20 - t50 * t23;
t29 = t50 * qJD(2) + t66 * t36;
t33 = qJD(4) + t34;
t7 = t33 * pkin(4) - t29 * pkin(9) + t10;
t11 = t50 * t20 + t66 * t23;
t27 = -t66 * qJD(2) + t50 * t36;
t9 = -t27 * pkin(9) + t11;
t4 = t49 * t7 + t65 * t9;
t15 = t65 * t27 + t49 * t29;
t17 = -t49 * t27 + t65 * t29;
t64 = t17 * t15;
t31 = qJD(5) + t33;
t63 = t31 * t15;
t62 = t52 * t53;
t58 = t15 ^ 2 / 0.2e1;
t57 = qJD(1) * qJD(2);
t56 = t51 * t57;
t55 = t52 * t57;
t24 = t48 * t39 - t47 * t40;
t3 = -t49 * t9 + t65 * t7;
t22 = -qJD(2) * pkin(3) - t24;
t13 = t27 * pkin(4) + t22;
t46 = t52 ^ 2;
t45 = t51 ^ 2;
t44 = qJD(2) ^ 2 / 0.2e1;
t30 = t31 ^ 2 / 0.2e1;
t14 = t17 ^ 2 / 0.2e1;
t12 = t17 * t31;
t5 = t15 * pkin(5) - t17 * qJ(6) + t13;
t2 = t31 * qJ(6) + t4;
t1 = -t31 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t67, 0, 0, 0, 0, t45 * t67, t51 * t62, t56, t46 * t67, t55, t44, pkin(1) * t62 - pkin(7) * t56, -t53 * pkin(1) * t51 - pkin(7) * t55 (t45 + t46) * t53 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t46 / 0.2e1 + t45 / 0.2e1) * pkin(7) ^ 2) * t53, t36 ^ 2 / 0.2e1, -t36 * t34, t36 * qJD(2), t34 ^ 2 / 0.2e1, -t34 * qJD(2), t44, t24 * qJD(2) + t41 * t34, -t25 * qJD(2) + t41 * t36, -t24 * t36 - t25 * t34, t25 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1 + t41 ^ 2 / 0.2e1, t29 ^ 2 / 0.2e1, -t29 * t27, t29 * t33, t27 ^ 2 / 0.2e1, -t27 * t33, t33 ^ 2 / 0.2e1, t10 * t33 + t22 * t27, -t11 * t33 + t22 * t29, -t10 * t29 - t11 * t27, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t14, -t64, t12, t58, -t63, t30, t13 * t15 + t3 * t31, t13 * t17 - t4 * t31, -t4 * t15 - t3 * t17, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1, t14, t12, t64, t30, t63, t58, -t1 * t31 + t5 * t15, t1 * t17 - t2 * t15, -t5 * t17 + t2 * t31, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
