% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRR10_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR10_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:10:43
% EndTime: 2019-12-31 19:10:43
% DurationCPUTime: 0.16s
% Computational Cost: add. (442->51), mult. (1135->122), div. (0->0), fcn. (797->8), ass. (0->40)
t45 = qJD(1) ^ 2;
t52 = t45 / 0.2e1;
t51 = cos(qJ(4));
t50 = cos(qJ(5));
t49 = pkin(6) + qJ(2);
t43 = sin(qJ(3));
t44 = cos(qJ(3));
t40 = cos(pkin(9));
t47 = qJD(1) * t40;
t39 = sin(pkin(9));
t48 = qJD(1) * t39;
t28 = t43 * t48 - t44 * t47;
t30 = (t39 * t44 + t40 * t43) * qJD(1);
t33 = qJD(2) + (-pkin(2) * t40 - pkin(1)) * qJD(1);
t13 = t28 * pkin(3) - t30 * pkin(7) + t33;
t31 = t49 * t48;
t32 = t49 * t47;
t18 = -t43 * t31 + t44 * t32;
t16 = qJD(3) * pkin(7) + t18;
t42 = sin(qJ(4));
t6 = t42 * t13 + t51 * t16;
t5 = t51 * t13 - t42 * t16;
t17 = -t44 * t31 - t43 * t32;
t24 = qJD(4) + t28;
t15 = -qJD(3) * pkin(3) - t17;
t41 = sin(qJ(5));
t38 = t40 ^ 2;
t37 = t39 ^ 2;
t35 = -qJD(1) * pkin(1) + qJD(2);
t23 = qJD(5) + t24;
t22 = t42 * qJD(3) + t51 * t30;
t20 = -t51 * qJD(3) + t42 * t30;
t10 = -t41 * t20 + t50 * t22;
t8 = t50 * t20 + t41 * t22;
t7 = t20 * pkin(4) + t15;
t4 = -t20 * pkin(8) + t6;
t3 = t24 * pkin(4) - t22 * pkin(8) + t5;
t2 = t41 * t3 + t50 * t4;
t1 = t50 * t3 - t41 * t4;
t9 = [0, 0, 0, 0, 0, t52, 0, 0, 0, 0, t37 * t52, t39 * t45 * t40, 0, t38 * t52, 0, 0, -t35 * t47, t35 * t48, (t37 + t38) * t45 * qJ(2), t35 ^ 2 / 0.2e1 + (t38 / 0.2e1 + t37 / 0.2e1) * qJ(2) ^ 2 * t45, t30 ^ 2 / 0.2e1, -t30 * t28, t30 * qJD(3), t28 ^ 2 / 0.2e1, -t28 * qJD(3), qJD(3) ^ 2 / 0.2e1, t17 * qJD(3) + t33 * t28, -t18 * qJD(3) + t33 * t30, -t17 * t30 - t18 * t28, t18 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * t24, t20 ^ 2 / 0.2e1, -t20 * t24, t24 ^ 2 / 0.2e1, t15 * t20 + t5 * t24, t15 * t22 - t6 * t24, -t6 * t20 - t5 * t22, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, t10 * t23, t8 ^ 2 / 0.2e1, -t8 * t23, t23 ^ 2 / 0.2e1, t1 * t23 + t7 * t8, t7 * t10 - t2 * t23, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg = t9;
