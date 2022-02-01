% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:25:52
% EndTime: 2022-01-23 09:25:52
% DurationCPUTime: 0.18s
% Computational Cost: add. (402->51), mult. (1069->124), div. (0->0), fcn. (707->8), ass. (0->46)
t45 = qJD(1) ^ 2;
t58 = t45 / 0.2e1;
t57 = cos(qJ(5));
t39 = sin(pkin(8));
t36 = t39 ^ 2;
t56 = t36 * t45;
t41 = cos(pkin(8));
t24 = qJD(2) + (-pkin(2) * t41 - pkin(6) * t39 - pkin(1)) * qJD(1);
t44 = cos(qJ(3));
t23 = t44 * t24;
t53 = t41 * qJD(1);
t32 = -qJD(3) + t53;
t43 = sin(qJ(3));
t12 = -t32 * pkin(3) + t23 + (-qJ(2) * t41 * t43 - qJ(4) * t39 * t44) * qJD(1);
t52 = qJ(2) * qJD(1);
t47 = t41 * t52;
t17 = t43 * t24 + t44 * t47;
t54 = qJD(1) * t39;
t50 = t43 * t54;
t15 = -qJ(4) * t50 + t17;
t38 = sin(pkin(9));
t40 = cos(pkin(9));
t6 = t38 * t12 + t40 * t15;
t55 = qJ(2) * t45;
t51 = t36 * t55;
t49 = t56 / 0.2e1;
t48 = qJ(2) ^ 2 * t58;
t25 = pkin(3) * t50 + t39 * t52 + qJD(4);
t5 = t40 * t12 - t38 * t15;
t42 = sin(qJ(5));
t37 = t41 ^ 2;
t35 = -qJD(1) * pkin(1) + qJD(2);
t31 = t36 * t48;
t28 = -qJD(5) + t32;
t26 = t32 ^ 2 / 0.2e1;
t21 = (-t38 * t43 + t40 * t44) * t54;
t19 = (-t38 * t44 - t40 * t43) * t54;
t16 = -t43 * t47 + t23;
t14 = -t19 * pkin(4) + t25;
t9 = t42 * t19 + t57 * t21;
t7 = -t57 * t19 + t42 * t21;
t4 = t19 * pkin(7) + t6;
t3 = -t32 * pkin(4) - t21 * pkin(7) + t5;
t2 = t42 * t3 + t57 * t4;
t1 = t57 * t3 - t42 * t4;
t8 = [0, 0, 0, 0, 0, t58, 0, 0, 0, 0, t49, t39 * t45 * t41, 0, t37 * t58, 0, 0, -t35 * t53, t35 * t54, (t36 + t37) * t55, t37 * t48 + t31 + t35 ^ 2 / 0.2e1, t44 ^ 2 * t49, -t44 * t43 * t56, -t44 * t32 * t54, t43 ^ 2 * t49, t32 * t50, t26, -t16 * t32 + t43 * t51, t17 * t32 + t44 * t51, (-t16 * t44 - t17 * t43) * t54, t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t31, t21 ^ 2 / 0.2e1, t21 * t19, -t21 * t32, t19 ^ 2 / 0.2e1, -t19 * t32, t26, -t25 * t19 - t5 * t32, t25 * t21 + t6 * t32, t6 * t19 - t5 * t21, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1, t9 ^ 2 / 0.2e1, -t9 * t7, -t9 * t28, t7 ^ 2 / 0.2e1, t7 * t28, t28 ^ 2 / 0.2e1, -t1 * t28 + t14 * t7, t14 * t9 + t2 * t28, -t1 * t9 - t2 * t7, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1;];
T_reg = t8;
