% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRPR4
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
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR4_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:11:38
% EndTime: 2019-12-31 21:11:38
% DurationCPUTime: 0.15s
% Computational Cost: add. (234->41), mult. (364->98), div. (0->0), fcn. (158->6), ass. (0->43)
t29 = sin(qJ(3));
t26 = t29 ^ 2;
t49 = t26 / 0.2e1;
t32 = cos(qJ(3));
t27 = t32 ^ 2;
t48 = t27 / 0.2e1;
t47 = pkin(3) + pkin(4);
t23 = qJD(1) + qJD(2);
t46 = t23 * t29;
t45 = t23 * t32;
t30 = sin(qJ(2));
t44 = pkin(1) * qJD(1);
t39 = t30 * t44;
t15 = t23 * pkin(7) + t39;
t11 = qJD(3) * qJ(4) + t32 * t15;
t43 = qJD(3) * t29;
t42 = qJD(3) * t32;
t41 = t29 * t15 + qJD(4);
t20 = t23 ^ 2;
t40 = t29 * t20 * t32;
t33 = cos(qJ(2));
t38 = t33 * t44;
t37 = t23 * t42;
t36 = qJ(4) * t29 + pkin(2);
t34 = qJD(1) ^ 2;
t31 = cos(qJ(5));
t28 = sin(qJ(5));
t24 = qJD(3) ^ 2 / 0.2e1;
t21 = qJD(3) - qJD(5);
t19 = t23 * t43;
t18 = t20 * t48;
t17 = t20 * t49;
t16 = -t23 * pkin(2) - t38;
t10 = (-t28 * t32 + t29 * t31) * t23;
t8 = (-t28 * t29 - t31 * t32) * t23;
t7 = -qJD(3) * pkin(3) + t41;
t6 = -pkin(8) * t45 + t11;
t5 = -t38 + (-pkin(3) * t32 - t36) * t23;
t4 = -pkin(8) * t46 - t47 * qJD(3) + t41;
t3 = t38 + (t47 * t32 + t36) * t23;
t2 = t28 * t4 + t31 * t6;
t1 = -t28 * t6 + t31 * t4;
t9 = [0, 0, 0, 0, 0, t34 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20 / 0.2e1, t23 * t38, -t23 * t39, 0, (t30 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t34, t17, t40, t19, t18, t37, t24, -t15 * t43 - t16 * t45, -t15 * t42 + t16 * t46, (t26 + t27) * t23 * t15, t16 ^ 2 / 0.2e1 + (t48 + t49) * t15 ^ 2, t17, t19, -t40, t24, -t37, t18, -t7 * qJD(3) - t5 * t45, (t11 * t32 + t29 * t7) * t23, t11 * qJD(3) - t5 * t46, t11 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, t10 * t8, -t10 * t21, t8 ^ 2 / 0.2e1, -t8 * t21, t21 ^ 2 / 0.2e1, -t1 * t21 - t3 * t8, t3 * t10 + t2 * t21, -t1 * t10 + t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t9;
