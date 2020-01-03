% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRPR12
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
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR12_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR12_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:30:14
% EndTime: 2019-12-31 18:30:14
% DurationCPUTime: 0.16s
% Computational Cost: add. (432->51), mult. (1135->119), div. (0->0), fcn. (797->8), ass. (0->40)
t43 = qJD(1) ^ 2;
t52 = t43 / 0.2e1;
t51 = cos(qJ(3));
t50 = cos(qJ(5));
t49 = pkin(6) + qJ(2);
t42 = sin(qJ(3));
t40 = cos(pkin(8));
t46 = qJD(1) * t40;
t39 = sin(pkin(8));
t47 = qJD(1) * t39;
t27 = t42 * t47 - t51 * t46;
t29 = (t51 * t39 + t40 * t42) * qJD(1);
t32 = qJD(2) + (-pkin(2) * t40 - pkin(1)) * qJD(1);
t13 = t27 * pkin(3) - t29 * qJ(4) + t32;
t30 = t49 * t47;
t31 = t49 * t46;
t18 = -t42 * t30 + t51 * t31;
t16 = qJD(3) * qJ(4) + t18;
t38 = sin(pkin(9));
t48 = cos(pkin(9));
t6 = t38 * t13 + t48 * t16;
t45 = t27 ^ 2 / 0.2e1;
t5 = t48 * t13 - t38 * t16;
t17 = -t51 * t30 - t42 * t31;
t15 = -qJD(3) * pkin(3) + qJD(4) - t17;
t41 = sin(qJ(5));
t37 = t40 ^ 2;
t36 = t39 ^ 2;
t34 = -qJD(1) * pkin(1) + qJD(2);
t23 = qJD(5) + t27;
t22 = t38 * qJD(3) + t48 * t29;
t20 = -t48 * qJD(3) + t38 * t29;
t10 = -t41 * t20 + t50 * t22;
t8 = t50 * t20 + t41 * t22;
t7 = t20 * pkin(4) + t15;
t4 = -t20 * pkin(7) + t6;
t3 = t27 * pkin(4) - t22 * pkin(7) + t5;
t2 = t41 * t3 + t50 * t4;
t1 = t50 * t3 - t41 * t4;
t9 = [0, 0, 0, 0, 0, t52, 0, 0, 0, 0, t36 * t52, t39 * t43 * t40, 0, t37 * t52, 0, 0, -t34 * t46, t34 * t47, (t36 + t37) * t43 * qJ(2), t34 ^ 2 / 0.2e1 + (t37 / 0.2e1 + t36 / 0.2e1) * qJ(2) ^ 2 * t43, t29 ^ 2 / 0.2e1, -t29 * t27, t29 * qJD(3), t45, -t27 * qJD(3), qJD(3) ^ 2 / 0.2e1, t17 * qJD(3) + t32 * t27, -t18 * qJD(3) + t32 * t29, -t17 * t29 - t18 * t27, t18 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * t27, t20 ^ 2 / 0.2e1, -t20 * t27, t45, t15 * t20 + t5 * t27, t15 * t22 - t6 * t27, -t6 * t20 - t5 * t22, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, t10 * t23, t8 ^ 2 / 0.2e1, -t8 * t23, t23 ^ 2 / 0.2e1, t1 * t23 + t7 * t8, t7 * t10 - t2 * t23, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg = t9;
