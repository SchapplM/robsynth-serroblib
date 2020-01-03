% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPR11_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR11_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:35:12
% EndTime: 2019-12-31 21:35:13
% DurationCPUTime: 0.17s
% Computational Cost: add. (312->50), mult. (708->111), div. (0->0), fcn. (415->6), ass. (0->45)
t39 = qJD(1) ^ 2;
t55 = t39 / 0.2e1;
t54 = pkin(3) + pkin(4);
t53 = cos(qJ(3));
t52 = cos(qJ(5));
t36 = sin(qJ(3));
t37 = sin(qJ(2));
t48 = qJD(1) * t37;
t19 = -t53 * qJD(2) + t36 * t48;
t21 = t36 * qJD(2) + t53 * t48;
t51 = t21 * t19;
t38 = cos(qJ(2));
t47 = t38 * qJD(1);
t29 = -qJD(3) + t47;
t50 = t29 * t19;
t49 = t38 * t39;
t17 = (-pkin(2) * t38 - pkin(7) * t37 - pkin(1)) * qJD(1);
t25 = pkin(6) * t47 + qJD(2) * pkin(7);
t13 = t36 * t17 + t53 * t25;
t46 = t19 ^ 2 / 0.2e1;
t45 = qJD(1) * qJD(2);
t7 = -t29 * qJ(4) + t13;
t44 = t37 * t45;
t43 = t38 * t45;
t24 = -qJD(2) * pkin(2) + pkin(6) * t48;
t12 = t53 * t17 - t36 * t25;
t42 = qJD(4) - t12;
t41 = t21 * qJ(4) - t24;
t35 = sin(qJ(5));
t34 = t38 ^ 2;
t33 = t37 ^ 2;
t28 = qJD(5) + t29;
t26 = t29 ^ 2 / 0.2e1;
t18 = t21 ^ 2 / 0.2e1;
t14 = t21 * t29;
t11 = t35 * t19 + t52 * t21;
t9 = -t52 * t19 + t35 * t21;
t8 = t19 * pkin(3) - t41;
t6 = t29 * pkin(3) + t42;
t5 = -t54 * t19 + t41;
t4 = t19 * pkin(8) + t7;
t3 = -t21 * pkin(8) + t54 * t29 + t42;
t2 = t35 * t3 + t52 * t4;
t1 = t52 * t3 - t35 * t4;
t10 = [0, 0, 0, 0, 0, t55, 0, 0, 0, 0, t33 * t55, t37 * t49, t44, t34 * t55, t43, qJD(2) ^ 2 / 0.2e1, pkin(1) * t49 - pkin(6) * t44, -t39 * pkin(1) * t37 - pkin(6) * t43, (t33 + t34) * t39 * pkin(6), (pkin(1) ^ 2 / 0.2e1 + (t34 / 0.2e1 + t33 / 0.2e1) * pkin(6) ^ 2) * t39, t18, -t51, -t14, t46, t50, t26, -t12 * t29 + t24 * t19, t13 * t29 + t24 * t21, -t12 * t21 - t13 * t19, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1, t18, -t14, t51, t26, -t50, t46, t8 * t19 + t6 * t29, -t7 * t19 + t6 * t21, -t8 * t21 - t7 * t29, t7 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t11 ^ 2 / 0.2e1, -t11 * t9, t11 * t28, t9 ^ 2 / 0.2e1, -t9 * t28, t28 ^ 2 / 0.2e1, t1 * t28 + t5 * t9, t5 * t11 - t2 * t28, -t1 * t11 - t2 * t9, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t10;
