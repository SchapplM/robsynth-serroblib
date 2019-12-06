% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRPR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR5_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:27:52
% EndTime: 2019-12-05 16:27:53
% DurationCPUTime: 0.18s
% Computational Cost: add. (259->46), mult. (677->118), div. (0->0), fcn. (472->10), ass. (0->43)
t39 = qJD(2) ^ 2;
t51 = t39 / 0.2e1;
t36 = sin(qJ(2));
t32 = sin(pkin(5));
t48 = qJD(1) * t32;
t23 = qJD(2) * pkin(7) + t36 * t48;
t37 = cos(qJ(3));
t33 = cos(pkin(5));
t47 = qJD(1) * t33;
t27 = t37 * t47;
t35 = sin(qJ(3));
t44 = qJ(4) * qJD(2);
t10 = qJD(3) * pkin(3) + t27 + (-t23 - t44) * t35;
t16 = t37 * t23 + t35 * t47;
t11 = t37 * t44 + t16;
t31 = sin(pkin(10));
t49 = cos(pkin(10));
t6 = t31 * t10 + t49 * t11;
t50 = cos(qJ(5));
t46 = qJD(2) * t35;
t45 = qJD(2) * t37;
t43 = qJD(2) * qJD(3);
t38 = cos(qJ(2));
t42 = t38 * t48;
t41 = qJD(2) * t48;
t19 = t31 * t46 - t49 * t45;
t5 = t49 * t10 - t31 * t11;
t17 = -t42 + qJD(4) + (-pkin(3) * t37 - pkin(2)) * qJD(2);
t40 = qJD(1) ^ 2;
t34 = sin(qJ(5));
t30 = qJD(3) ^ 2 / 0.2e1;
t24 = -qJD(2) * pkin(2) - t42;
t21 = (t31 * t37 + t49 * t35) * qJD(2);
t18 = qJD(5) + t19;
t15 = -t35 * t23 + t27;
t14 = t34 * qJD(3) + t50 * t21;
t12 = -t50 * qJD(3) + t34 * t21;
t7 = t19 * pkin(4) - t21 * pkin(8) + t17;
t4 = qJD(3) * pkin(8) + t6;
t3 = -qJD(3) * pkin(4) - t5;
t2 = t34 * t7 + t50 * t4;
t1 = -t34 * t4 + t50 * t7;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t40 / 0.2e1, 0, 0, 0, 0, 0, t51, t38 * t41, -t36 * t41, 0, (t33 ^ 2 / 0.2e1 + (t36 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1) * t32 ^ 2) * t40, t35 ^ 2 * t51, t35 * t39 * t37, t35 * t43, t37 ^ 2 * t51, t37 * t43, t30, t15 * qJD(3) - t24 * t45, -t16 * qJD(3) + t24 * t46, (-t15 * t35 + t16 * t37) * qJD(2), t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1, t21 ^ 2 / 0.2e1, -t21 * t19, t21 * qJD(3), t19 ^ 2 / 0.2e1, -t19 * qJD(3), t30, t5 * qJD(3) + t17 * t19, -t6 * qJD(3) + t17 * t21, -t6 * t19 - t5 * t21, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t14 ^ 2 / 0.2e1, -t14 * t12, t14 * t18, t12 ^ 2 / 0.2e1, -t12 * t18, t18 ^ 2 / 0.2e1, t1 * t18 + t3 * t12, t3 * t14 - t2 * t18, -t1 * t14 - t2 * t12, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t8;
