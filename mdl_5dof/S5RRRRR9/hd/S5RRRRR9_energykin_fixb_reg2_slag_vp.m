% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRR9_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR9_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:29:51
% EndTime: 2019-12-31 22:29:51
% DurationCPUTime: 0.18s
% Computational Cost: add. (559->53), mult. (1269->131), div. (0->0), fcn. (856->8), ass. (0->43)
t48 = qJD(1) ^ 2;
t59 = t48 / 0.2e1;
t58 = cos(qJ(3));
t57 = cos(qJ(4));
t56 = cos(qJ(5));
t47 = cos(qJ(2));
t55 = t47 * t48;
t46 = sin(qJ(2));
t26 = (-pkin(2) * t47 - pkin(7) * t46 - pkin(1)) * qJD(1);
t53 = t47 * qJD(1);
t34 = pkin(6) * t53 + qJD(2) * pkin(7);
t45 = sin(qJ(3));
t20 = t58 * t26 - t45 * t34;
t54 = qJD(1) * t46;
t29 = t45 * qJD(2) + t58 * t54;
t37 = -qJD(3) + t53;
t13 = -t37 * pkin(3) - t29 * pkin(8) + t20;
t21 = t45 * t26 + t58 * t34;
t27 = -t58 * qJD(2) + t45 * t54;
t15 = -t27 * pkin(8) + t21;
t44 = sin(qJ(4));
t6 = t44 * t13 + t57 * t15;
t52 = qJD(1) * qJD(2);
t51 = t46 * t52;
t50 = t47 * t52;
t5 = t57 * t13 - t44 * t15;
t33 = -qJD(2) * pkin(2) + pkin(6) * t54;
t35 = -qJD(4) + t37;
t22 = t27 * pkin(3) + t33;
t43 = sin(qJ(5));
t42 = t47 ^ 2;
t41 = t46 ^ 2;
t31 = -qJD(5) + t35;
t19 = -t44 * t27 + t57 * t29;
t17 = t57 * t27 + t44 * t29;
t10 = t17 * pkin(4) + t22;
t9 = -t43 * t17 + t56 * t19;
t7 = t56 * t17 + t43 * t19;
t4 = -t17 * pkin(9) + t6;
t3 = -t35 * pkin(4) - t19 * pkin(9) + t5;
t2 = t43 * t3 + t56 * t4;
t1 = t56 * t3 - t43 * t4;
t8 = [0, 0, 0, 0, 0, t59, 0, 0, 0, 0, t41 * t59, t46 * t55, t51, t42 * t59, t50, qJD(2) ^ 2 / 0.2e1, pkin(1) * t55 - pkin(6) * t51, -t48 * pkin(1) * t46 - pkin(6) * t50, (t41 + t42) * t48 * pkin(6), (pkin(1) ^ 2 / 0.2e1 + (t42 / 0.2e1 + t41 / 0.2e1) * pkin(6) ^ 2) * t48, t29 ^ 2 / 0.2e1, -t29 * t27, -t29 * t37, t27 ^ 2 / 0.2e1, t27 * t37, t37 ^ 2 / 0.2e1, -t20 * t37 + t33 * t27, t21 * t37 + t33 * t29, -t20 * t29 - t21 * t27, t21 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1, t19 ^ 2 / 0.2e1, -t19 * t17, -t35 * t19, t17 ^ 2 / 0.2e1, t35 * t17, t35 ^ 2 / 0.2e1, t22 * t17 - t5 * t35, t22 * t19 + t6 * t35, -t6 * t17 - t5 * t19, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t9 ^ 2 / 0.2e1, -t9 * t7, -t9 * t31, t7 ^ 2 / 0.2e1, t7 * t31, t31 ^ 2 / 0.2e1, -t1 * t31 + t10 * t7, t10 * t9 + t2 * t31, -t1 * t9 - t2 * t7, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg = t8;
