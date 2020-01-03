% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP10_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP10_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:52:03
% EndTime: 2019-12-31 18:52:03
% DurationCPUTime: 0.15s
% Computational Cost: add. (311->47), mult. (830->102), div. (0->0), fcn. (555->6), ass. (0->39)
t42 = qJD(1) ^ 2;
t49 = t42 / 0.2e1;
t38 = sin(pkin(8));
t45 = qJD(1) * t38;
t46 = pkin(6) + qJ(2);
t30 = t46 * t45;
t39 = cos(pkin(8));
t44 = qJD(1) * t39;
t31 = t46 * t44;
t41 = sin(qJ(3));
t48 = cos(qJ(3));
t16 = -t41 * t30 + t48 * t31;
t14 = qJD(3) * pkin(7) + t16;
t40 = sin(qJ(4));
t47 = cos(qJ(4));
t27 = t41 * t45 - t48 * t44;
t29 = (t48 * t38 + t39 * t41) * qJD(1);
t32 = qJD(2) + (-pkin(2) * t39 - pkin(1)) * qJD(1);
t9 = t27 * pkin(3) - t29 * pkin(7) + t32;
t4 = t47 * t14 + t40 * t9;
t3 = -t40 * t14 + t47 * t9;
t15 = -t48 * t30 - t41 * t31;
t13 = -qJD(3) * pkin(3) - t15;
t37 = t39 ^ 2;
t36 = t38 ^ 2;
t34 = -qJD(1) * pkin(1) + qJD(2);
t23 = qJD(4) + t27;
t22 = t23 ^ 2 / 0.2e1;
t21 = t40 * qJD(3) + t47 * t29;
t19 = -t47 * qJD(3) + t40 * t29;
t18 = t21 ^ 2 / 0.2e1;
t17 = t19 ^ 2 / 0.2e1;
t11 = t21 * t23;
t10 = t19 * t23;
t6 = t21 * t19;
t5 = t19 * pkin(4) + qJD(5) + t13;
t2 = -t19 * qJ(5) + t4;
t1 = t23 * pkin(4) - t21 * qJ(5) + t3;
t7 = [0, 0, 0, 0, 0, t49, 0, 0, 0, 0, t36 * t49, t38 * t42 * t39, 0, t37 * t49, 0, 0, -t34 * t44, t34 * t45, (t36 + t37) * t42 * qJ(2), t34 ^ 2 / 0.2e1 + (t37 / 0.2e1 + t36 / 0.2e1) * qJ(2) ^ 2 * t42, t29 ^ 2 / 0.2e1, -t29 * t27, t29 * qJD(3), t27 ^ 2 / 0.2e1, -t27 * qJD(3), qJD(3) ^ 2 / 0.2e1, t15 * qJD(3) + t32 * t27, -t16 * qJD(3) + t32 * t29, -t15 * t29 - t16 * t27, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1, t18, -t6, t11, t17, -t10, t22, t13 * t19 + t3 * t23, t13 * t21 - t4 * t23, -t4 * t19 - t3 * t21, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1, t18, -t6, t11, t17, -t10, t22, t1 * t23 + t5 * t19, -t2 * t23 + t5 * t21, -t1 * t21 - t2 * t19, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t7;
