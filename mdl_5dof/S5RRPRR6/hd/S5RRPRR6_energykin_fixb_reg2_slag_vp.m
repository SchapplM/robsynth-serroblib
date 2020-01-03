% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:05:53
% EndTime: 2020-01-03 12:05:53
% DurationCPUTime: 0.16s
% Computational Cost: add. (355->39), mult. (538->106), div. (0->0), fcn. (294->8), ass. (0->45)
t28 = cos(pkin(9));
t25 = t28 ^ 2;
t51 = t25 / 0.2e1;
t30 = sin(qJ(4));
t33 = cos(qJ(4));
t26 = qJD(1) + qJD(2);
t31 = sin(qJ(2));
t44 = pkin(1) * qJD(1);
t40 = t31 * t44;
t18 = t26 * qJ(3) + t40;
t49 = t18 * t28;
t27 = sin(pkin(9));
t34 = cos(qJ(2));
t39 = t34 * t44;
t37 = qJD(3) - t39;
t9 = (-pkin(3) * t28 - pkin(7) * t27 - pkin(2)) * t26 + t37;
t6 = t30 * t9 + t33 * t49;
t50 = t18 * t26;
t23 = t26 ^ 2;
t24 = t27 ^ 2;
t48 = t24 * t23;
t47 = t26 * t27;
t46 = t26 * t30;
t45 = t28 * t26;
t43 = t24 * t50;
t42 = t27 * t46;
t41 = t33 * t47;
t38 = t48 / 0.2e1;
t21 = -qJD(4) + t45;
t5 = -t30 * t49 + t33 * t9;
t35 = qJD(1) ^ 2;
t32 = cos(qJ(5));
t29 = sin(qJ(5));
t19 = -qJD(5) + t21;
t17 = t18 ^ 2;
t16 = -t26 * pkin(2) + t37;
t15 = t24 * t17 / 0.2e1;
t13 = (-t29 * t30 + t32 * t33) * t47;
t11 = (-t29 * t33 - t30 * t32) * t47;
t10 = (pkin(4) * t46 + t18) * t27;
t4 = -pkin(8) * t42 + t6;
t3 = -t21 * pkin(4) - pkin(8) * t41 + t5;
t2 = t29 * t3 + t32 * t4;
t1 = -t29 * t4 + t32 * t3;
t7 = [0, 0, 0, 0, 0, t35 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 / 0.2e1, t26 * t39, -t26 * t40, 0, (t31 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t35, t38, t27 * t23 * t28, 0, t23 * t51, 0, 0, -t16 * t45, t16 * t47, (t24 + t25) * t50, t17 * t51 + t15 + t16 ^ 2 / 0.2e1, t33 ^ 2 * t38, -t33 * t30 * t48, -t21 * t41, t30 ^ 2 * t38, t21 * t42, t21 ^ 2 / 0.2e1, -t5 * t21 + t30 * t43, t6 * t21 + t33 * t43, (-t30 * t6 - t33 * t5) * t47, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t15, t13 ^ 2 / 0.2e1, t13 * t11, -t13 * t19, t11 ^ 2 / 0.2e1, -t11 * t19, t19 ^ 2 / 0.2e1, -t1 * t19 - t10 * t11, t10 * t13 + t2 * t19, -t1 * t13 + t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg = t7;
