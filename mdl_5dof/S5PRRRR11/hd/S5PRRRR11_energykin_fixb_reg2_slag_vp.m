% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRR11_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR11_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:48
% EndTime: 2024-09-27 21:45:48
% DurationCPUTime: 0.09s
% Computational Cost: add. (328->44), mult. (998->108), div. (0->0), fcn. (750->8), ass. (0->40)
t36 = sin(pkin(5));
t37 = cos(pkin(5));
t47 = t37 * qJD(2);
t54 = pkin(2) * t47 + qJD(1) * t36;
t53 = cos(qJ(4));
t52 = cos(qJ(5));
t42 = qJD(2) ^ 2;
t51 = t36 ^ 2 * t42;
t32 = qJD(3) + t47;
t40 = sin(qJ(3));
t48 = qJD(2) * t36;
t45 = t40 * t48;
t41 = cos(qJ(3));
t50 = t54 * t41;
t13 = t32 * pkin(3) + (-pkin(7) - pkin(8)) * t45 + t50;
t44 = t41 * t48;
t17 = pkin(7) * t44 + t54 * t40;
t15 = pkin(8) * t44 + t17;
t39 = sin(qJ(4));
t6 = t39 * t13 + t53 * t15;
t43 = t51 / 0.2e1;
t5 = t53 * t13 - t39 * t15;
t26 = qJD(4) + t32;
t33 = t37 * qJD(1);
t22 = t33 + (-pkin(3) * t41 - pkin(2)) * t48;
t38 = sin(qJ(5));
t35 = qJD(1) ^ 2 / 0.2e1;
t25 = qJD(5) + t26;
t23 = -pkin(2) * t48 + t33;
t21 = (t39 * t41 + t53 * t40) * t48;
t19 = t39 * t45 - t53 * t44;
t16 = -pkin(7) * t45 + t50;
t12 = t19 * pkin(4) + t22;
t9 = -t38 * t19 + t52 * t21;
t7 = t52 * t19 + t38 * t21;
t4 = -t19 * pkin(9) + t6;
t3 = t26 * pkin(4) - t21 * pkin(9) + t5;
t2 = t38 * t3 + t52 * t4;
t1 = t52 * t3 - t38 * t4;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, t42 / 0.2e1, 0, 0, 0, t35, t40 ^ 2 * t43, t40 * t41 * t51, t32 * t45, t41 ^ 2 * t43, t32 * t44, t32 ^ 2 / 0.2e1, t16 * t32 - t23 * t44, -t17 * t32 + t23 * t45, (-t16 * t40 + t17 * t41) * t48, t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1, t21 ^ 2 / 0.2e1, -t21 * t19, t21 * t26, t19 ^ 2 / 0.2e1, -t19 * t26, t26 ^ 2 / 0.2e1, t22 * t19 + t5 * t26, t22 * t21 - t6 * t26, -t6 * t19 - t5 * t21, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t9 ^ 2 / 0.2e1, -t9 * t7, t9 * t25, t7 ^ 2 / 0.2e1, -t7 * t25, t25 ^ 2 / 0.2e1, t1 * t25 + t12 * t7, t12 * t9 - t2 * t25, -t1 * t9 - t2 * t7, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1;];
T_reg = t8;
