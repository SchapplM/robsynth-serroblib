% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR15_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR15_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:43:07
% EndTime: 2019-12-31 20:43:07
% DurationCPUTime: 0.17s
% Computational Cost: add. (296->50), mult. (667->115), div. (0->0), fcn. (357->6), ass. (0->43)
t41 = qJD(1) ^ 2;
t53 = t41 / 0.2e1;
t52 = -pkin(2) - pkin(7);
t51 = cos(qJ(5));
t40 = cos(qJ(2));
t50 = t40 * t41;
t38 = sin(qJ(2));
t43 = -qJ(3) * t38 - pkin(1);
t14 = (t52 * t40 + t43) * qJD(1);
t48 = t38 * qJD(1);
t47 = pkin(6) * t48 + qJD(3);
t15 = pkin(3) * t48 + t52 * qJD(2) + t47;
t37 = sin(qJ(4));
t39 = cos(qJ(4));
t6 = t39 * t14 + t37 * t15;
t49 = qJD(1) * t40;
t23 = -pkin(6) * t49 - qJD(2) * qJ(3);
t46 = qJD(1) * qJD(2);
t17 = pkin(3) * t49 - t23;
t45 = t38 * t46;
t44 = t40 * t46;
t5 = -t37 * t14 + t39 * t15;
t26 = qJD(4) + t48;
t36 = sin(qJ(5));
t35 = t40 ^ 2;
t34 = t38 ^ 2;
t32 = qJD(2) ^ 2 / 0.2e1;
t28 = t35 * t53;
t27 = t34 * t53;
t25 = t38 * t50;
t24 = qJD(5) + t26;
t22 = -qJD(2) * pkin(2) + t47;
t21 = t39 * qJD(2) - t37 * t49;
t19 = t37 * qJD(2) + t39 * t49;
t18 = (-pkin(2) * t40 + t43) * qJD(1);
t10 = t19 * pkin(4) + t17;
t9 = -t36 * t19 + t51 * t21;
t7 = t51 * t19 + t36 * t21;
t4 = -t19 * pkin(8) + t6;
t3 = t26 * pkin(4) - t21 * pkin(8) + t5;
t2 = t36 * t3 + t51 * t4;
t1 = t51 * t3 - t36 * t4;
t8 = [0, 0, 0, 0, 0, t53, 0, 0, 0, 0, t27, t25, t45, t28, t44, t32, pkin(1) * t50 - pkin(6) * t45, -t41 * pkin(1) * t38 - pkin(6) * t44, (t34 + t35) * t41 * pkin(6), (pkin(1) ^ 2 / 0.2e1 + (t35 / 0.2e1 + t34 / 0.2e1) * pkin(6) ^ 2) * t41, t32, -t45, -t44, t27, t25, t28, (t22 * t38 - t23 * t40) * qJD(1), t22 * qJD(2) + t18 * t49, -t23 * qJD(2) - t18 * t48, t18 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t21 ^ 2 / 0.2e1, -t21 * t19, t21 * t26, t19 ^ 2 / 0.2e1, -t19 * t26, t26 ^ 2 / 0.2e1, t17 * t19 + t5 * t26, t17 * t21 - t6 * t26, -t6 * t19 - t5 * t21, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t9 ^ 2 / 0.2e1, -t9 * t7, t9 * t24, t7 ^ 2 / 0.2e1, -t7 * t24, t24 ^ 2 / 0.2e1, t1 * t24 + t10 * t7, t10 * t9 - t2 * t24, -t1 * t9 - t2 * t7, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg = t8;
