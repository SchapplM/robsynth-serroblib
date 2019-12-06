% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRR9_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR9_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:21:09
% EndTime: 2019-12-05 17:21:09
% DurationCPUTime: 0.17s
% Computational Cost: add. (280->45), mult. (680->119), div. (0->0), fcn. (466->10), ass. (0->41)
t42 = qJD(2) ^ 2;
t53 = t42 / 0.2e1;
t52 = cos(qJ(4));
t51 = cos(qJ(5));
t39 = sin(qJ(2));
t34 = sin(pkin(5));
t50 = qJD(1) * t34;
t25 = qJD(2) * pkin(7) + t39 * t50;
t38 = sin(qJ(3));
t40 = cos(qJ(3));
t35 = cos(pkin(5));
t49 = qJD(1) * t35;
t17 = t40 * t25 + t38 * t49;
t13 = qJD(3) * pkin(8) + t17;
t41 = cos(qJ(2));
t45 = t41 * t50;
t18 = -t45 + (-pkin(3) * t40 - pkin(8) * t38 - pkin(2)) * qJD(2);
t37 = sin(qJ(4));
t6 = t52 * t13 + t37 * t18;
t48 = qJD(2) * t38;
t47 = t40 * qJD(2);
t46 = qJD(2) * qJD(3);
t44 = qJD(2) * t50;
t5 = -t37 * t13 + t52 * t18;
t30 = -qJD(4) + t47;
t16 = -t38 * t25 + t40 * t49;
t12 = -qJD(3) * pkin(3) - t16;
t43 = qJD(1) ^ 2;
t36 = sin(qJ(5));
t27 = -qJD(5) + t30;
t26 = -qJD(2) * pkin(2) - t45;
t24 = t37 * qJD(3) + t52 * t48;
t22 = -t52 * qJD(3) + t37 * t48;
t10 = -t36 * t22 + t51 * t24;
t8 = t51 * t22 + t36 * t24;
t7 = t22 * pkin(4) + t12;
t4 = -t22 * pkin(9) + t6;
t3 = -t30 * pkin(4) - t24 * pkin(9) + t5;
t2 = t36 * t3 + t51 * t4;
t1 = t51 * t3 - t36 * t4;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t43 / 0.2e1, 0, 0, 0, 0, 0, t53, t41 * t44, -t39 * t44, 0, (t35 ^ 2 / 0.2e1 + (t39 ^ 2 / 0.2e1 + t41 ^ 2 / 0.2e1) * t34 ^ 2) * t43, t38 ^ 2 * t53, t38 * t42 * t40, t38 * t46, t40 ^ 2 * t53, t40 * t46, qJD(3) ^ 2 / 0.2e1, t16 * qJD(3) - t26 * t47, -t17 * qJD(3) + t26 * t48, (-t16 * t38 + t17 * t40) * qJD(2), t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1, t24 ^ 2 / 0.2e1, -t24 * t22, -t24 * t30, t22 ^ 2 / 0.2e1, t22 * t30, t30 ^ 2 / 0.2e1, t12 * t22 - t5 * t30, t12 * t24 + t6 * t30, -t6 * t22 - t5 * t24, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, -t10 * t27, t8 ^ 2 / 0.2e1, t8 * t27, t27 ^ 2 / 0.2e1, -t1 * t27 + t7 * t8, t7 * t10 + t2 * t27, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg = t9;
