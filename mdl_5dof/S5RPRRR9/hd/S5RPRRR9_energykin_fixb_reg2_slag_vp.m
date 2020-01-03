% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRRR9
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
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRR9_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR9_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:07:58
% EndTime: 2019-12-31 19:07:58
% DurationCPUTime: 0.16s
% Computational Cost: add. (449->51), mult. (1204->122), div. (0->0), fcn. (866->8), ass. (0->40)
t43 = qJD(1) ^ 2;
t51 = t43 / 0.2e1;
t38 = sin(pkin(9));
t46 = qJD(1) * t38;
t47 = pkin(6) + qJ(2);
t29 = t47 * t46;
t39 = cos(pkin(9));
t45 = qJD(1) * t39;
t30 = t47 * t45;
t42 = sin(qJ(3));
t50 = cos(qJ(3));
t19 = -t50 * t29 - t42 * t30;
t28 = (t50 * t38 + t39 * t42) * qJD(1);
t10 = qJD(3) * pkin(3) - t28 * pkin(7) + t19;
t20 = -t42 * t29 + t50 * t30;
t26 = t42 * t46 - t50 * t45;
t11 = -t26 * pkin(7) + t20;
t41 = sin(qJ(4));
t49 = cos(qJ(4));
t6 = t41 * t10 + t49 * t11;
t48 = cos(qJ(5));
t16 = t49 * t26 + t41 * t28;
t5 = t49 * t10 - t41 * t11;
t31 = qJD(2) + (-pkin(2) * t39 - pkin(1)) * qJD(1);
t21 = t26 * pkin(3) + t31;
t40 = sin(qJ(5));
t37 = qJD(3) + qJD(4);
t36 = t39 ^ 2;
t35 = t38 ^ 2;
t34 = -qJD(1) * pkin(1) + qJD(2);
t18 = -t41 * t26 + t49 * t28;
t15 = qJD(5) + t16;
t14 = t48 * t18 + t40 * t37;
t12 = t40 * t18 - t48 * t37;
t7 = t16 * pkin(4) - t18 * pkin(8) + t21;
t4 = t37 * pkin(8) + t6;
t3 = -t37 * pkin(4) - t5;
t2 = t48 * t4 + t40 * t7;
t1 = -t40 * t4 + t48 * t7;
t8 = [0, 0, 0, 0, 0, t51, 0, 0, 0, 0, t35 * t51, t38 * t43 * t39, 0, t36 * t51, 0, 0, -t34 * t45, t34 * t46, (t35 + t36) * t43 * qJ(2), t34 ^ 2 / 0.2e1 + (t36 / 0.2e1 + t35 / 0.2e1) * qJ(2) ^ 2 * t43, t28 ^ 2 / 0.2e1, -t28 * t26, t28 * qJD(3), t26 ^ 2 / 0.2e1, -t26 * qJD(3), qJD(3) ^ 2 / 0.2e1, t19 * qJD(3) + t31 * t26, -t20 * qJD(3) + t31 * t28, -t19 * t28 - t20 * t26, t20 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t37, t16 ^ 2 / 0.2e1, -t16 * t37, t37 ^ 2 / 0.2e1, t21 * t16 + t5 * t37, t21 * t18 - t6 * t37, -t6 * t16 - t5 * t18, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t14 ^ 2 / 0.2e1, -t14 * t12, t14 * t15, t12 ^ 2 / 0.2e1, -t12 * t15, t15 ^ 2 / 0.2e1, t1 * t15 + t3 * t12, t3 * t14 - t2 * t15, -t1 * t14 - t2 * t12, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t8;
