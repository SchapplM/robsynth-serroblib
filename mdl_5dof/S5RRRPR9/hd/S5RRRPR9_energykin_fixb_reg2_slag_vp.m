% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPR9_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR9_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:24:37
% EndTime: 2019-12-31 21:24:37
% DurationCPUTime: 0.18s
% Computational Cost: add. (549->53), mult. (1269->128), div. (0->0), fcn. (856->8), ass. (0->43)
t47 = qJD(1) ^ 2;
t58 = t47 / 0.2e1;
t57 = cos(qJ(3));
t56 = cos(qJ(5));
t46 = cos(qJ(2));
t55 = t46 * t47;
t45 = sin(qJ(2));
t26 = (-pkin(2) * t46 - pkin(7) * t45 - pkin(1)) * qJD(1);
t52 = t46 * qJD(1);
t32 = pkin(6) * t52 + qJD(2) * pkin(7);
t44 = sin(qJ(3));
t20 = t57 * t26 - t44 * t32;
t53 = qJD(1) * t45;
t29 = t44 * qJD(2) + t57 * t53;
t36 = -qJD(3) + t52;
t13 = -t36 * pkin(3) - t29 * qJ(4) + t20;
t21 = t44 * t26 + t57 * t32;
t27 = -t57 * qJD(2) + t44 * t53;
t15 = -t27 * qJ(4) + t21;
t42 = sin(pkin(9));
t54 = cos(pkin(9));
t6 = t42 * t13 + t54 * t15;
t51 = qJD(1) * qJD(2);
t50 = t45 * t51;
t49 = t46 * t51;
t5 = t54 * t13 - t42 * t15;
t31 = -qJD(2) * pkin(2) + pkin(6) * t53;
t22 = t27 * pkin(3) + qJD(4) + t31;
t43 = sin(qJ(5));
t41 = t46 ^ 2;
t40 = t45 ^ 2;
t34 = -qJD(5) + t36;
t33 = t36 ^ 2 / 0.2e1;
t19 = -t42 * t27 + t54 * t29;
t17 = t54 * t27 + t42 * t29;
t10 = t17 * pkin(4) + t22;
t9 = -t43 * t17 + t56 * t19;
t7 = t56 * t17 + t43 * t19;
t4 = -t17 * pkin(8) + t6;
t3 = -t36 * pkin(4) - t19 * pkin(8) + t5;
t2 = t43 * t3 + t56 * t4;
t1 = t56 * t3 - t43 * t4;
t8 = [0, 0, 0, 0, 0, t58, 0, 0, 0, 0, t40 * t58, t45 * t55, t50, t41 * t58, t49, qJD(2) ^ 2 / 0.2e1, pkin(1) * t55 - pkin(6) * t50, -t47 * pkin(1) * t45 - pkin(6) * t49, (t40 + t41) * t47 * pkin(6), (pkin(1) ^ 2 / 0.2e1 + (t41 / 0.2e1 + t40 / 0.2e1) * pkin(6) ^ 2) * t47, t29 ^ 2 / 0.2e1, -t29 * t27, -t29 * t36, t27 ^ 2 / 0.2e1, t27 * t36, t33, -t20 * t36 + t31 * t27, t21 * t36 + t31 * t29, -t20 * t29 - t21 * t27, t21 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1, t19 ^ 2 / 0.2e1, -t19 * t17, -t36 * t19, t17 ^ 2 / 0.2e1, t36 * t17, t33, t22 * t17 - t5 * t36, t22 * t19 + t6 * t36, -t6 * t17 - t5 * t19, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t9 ^ 2 / 0.2e1, -t9 * t7, -t9 * t34, t7 ^ 2 / 0.2e1, t7 * t34, t34 ^ 2 / 0.2e1, -t1 * t34 + t10 * t7, t10 * t9 + t2 * t34, -t1 * t9 - t2 * t7, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg = t8;
