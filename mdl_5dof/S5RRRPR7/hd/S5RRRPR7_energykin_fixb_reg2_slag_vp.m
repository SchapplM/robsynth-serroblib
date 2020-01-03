% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRPR7
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
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPR7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR7_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:17:15
% EndTime: 2019-12-31 21:17:15
% DurationCPUTime: 0.17s
% Computational Cost: add. (513->52), mult. (1178->127), div. (0->0), fcn. (801->8), ass. (0->44)
t43 = qJD(1) ^ 2;
t56 = t43 / 0.2e1;
t55 = -pkin(7) - pkin(6);
t54 = cos(qJ(3));
t53 = cos(qJ(5));
t42 = cos(qJ(2));
t52 = t42 * t43;
t40 = sin(qJ(3));
t49 = qJD(1) * t42;
t41 = sin(qJ(2));
t50 = qJD(1) * t41;
t25 = t40 * t50 - t54 * t49;
t27 = (t40 * t42 + t54 * t41) * qJD(1);
t32 = (-pkin(2) * t42 - pkin(1)) * qJD(1);
t13 = t25 * pkin(3) - t27 * qJ(4) + t32;
t30 = qJD(2) * pkin(2) + t55 * t50;
t31 = t55 * t49;
t18 = t40 * t30 - t54 * t31;
t35 = qJD(2) + qJD(3);
t16 = t35 * qJ(4) + t18;
t38 = sin(pkin(9));
t51 = cos(pkin(9));
t6 = t38 * t13 + t51 * t16;
t48 = t25 ^ 2 / 0.2e1;
t47 = qJD(1) * qJD(2);
t46 = t41 * t47;
t45 = t42 * t47;
t5 = t51 * t13 - t38 * t16;
t17 = t54 * t30 + t40 * t31;
t15 = -t35 * pkin(3) + qJD(4) - t17;
t39 = sin(qJ(5));
t37 = t42 ^ 2;
t36 = t41 ^ 2;
t24 = qJD(5) + t25;
t22 = t51 * t27 + t38 * t35;
t20 = t38 * t27 - t51 * t35;
t10 = -t39 * t20 + t53 * t22;
t8 = t53 * t20 + t39 * t22;
t7 = t20 * pkin(4) + t15;
t4 = -t20 * pkin(8) + t6;
t3 = t25 * pkin(4) - t22 * pkin(8) + t5;
t2 = t39 * t3 + t53 * t4;
t1 = t53 * t3 - t39 * t4;
t9 = [0, 0, 0, 0, 0, t56, 0, 0, 0, 0, t36 * t56, t41 * t52, t46, t37 * t56, t45, qJD(2) ^ 2 / 0.2e1, pkin(1) * t52 - pkin(6) * t46, -t43 * pkin(1) * t41 - pkin(6) * t45, (t36 + t37) * t43 * pkin(6), (pkin(1) ^ 2 / 0.2e1 + (t37 / 0.2e1 + t36 / 0.2e1) * pkin(6) ^ 2) * t43, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t35, t48, -t25 * t35, t35 ^ 2 / 0.2e1, t17 * t35 + t32 * t25, -t18 * t35 + t32 * t27, -t17 * t27 - t18 * t25, t18 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * t25, t20 ^ 2 / 0.2e1, -t20 * t25, t48, t15 * t20 + t5 * t25, t15 * t22 - t6 * t25, -t6 * t20 - t5 * t22, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, t10 * t24, t8 ^ 2 / 0.2e1, -t8 * t24, t24 ^ 2 / 0.2e1, t1 * t24 + t7 * t8, t7 * t10 - t2 * t24, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg = t9;
