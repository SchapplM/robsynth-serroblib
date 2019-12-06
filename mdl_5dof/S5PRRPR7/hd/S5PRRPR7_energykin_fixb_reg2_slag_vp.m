% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRRPR7
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
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRPR7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR7_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:37:32
% EndTime: 2019-12-05 16:37:32
% DurationCPUTime: 0.19s
% Computational Cost: add. (255->45), mult. (638->116), div. (0->0), fcn. (432->10), ass. (0->41)
t39 = qJD(2) ^ 2;
t49 = t39 / 0.2e1;
t35 = sin(qJ(2));
t31 = sin(pkin(5));
t47 = qJD(1) * t31;
t24 = qJD(2) * pkin(7) + t35 * t47;
t34 = sin(qJ(3));
t37 = cos(qJ(3));
t32 = cos(pkin(5));
t46 = qJD(1) * t32;
t16 = t37 * t24 + t34 * t46;
t11 = qJD(3) * qJ(4) + t16;
t38 = cos(qJ(2));
t42 = t38 * t47;
t17 = -t42 + (-pkin(3) * t37 - qJ(4) * t34 - pkin(2)) * qJD(2);
t30 = sin(pkin(10));
t48 = cos(pkin(10));
t7 = t48 * t11 + t30 * t17;
t45 = qJD(2) * t34;
t44 = qJD(2) * t37;
t43 = qJD(2) * qJD(3);
t41 = qJD(2) * t47;
t15 = -t34 * t24 + t37 * t46;
t21 = -t48 * qJD(3) + t30 * t45;
t6 = -t30 * t11 + t48 * t17;
t9 = -qJD(3) * pkin(3) + qJD(4) - t15;
t40 = qJD(1) ^ 2;
t36 = cos(qJ(5));
t33 = sin(qJ(5));
t27 = t37 ^ 2 * t49;
t25 = -qJD(2) * pkin(2) - t42;
t23 = t30 * qJD(3) + t48 * t45;
t18 = qJD(5) + t21;
t14 = t36 * t23 - t33 * t44;
t12 = t33 * t23 + t36 * t44;
t5 = t21 * pkin(4) - t23 * pkin(8) + t9;
t4 = -pkin(8) * t44 + t7;
t3 = pkin(4) * t44 - t6;
t2 = t33 * t5 + t36 * t4;
t1 = -t33 * t4 + t36 * t5;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t40 / 0.2e1, 0, 0, 0, 0, 0, t49, t38 * t41, -t35 * t41, 0, (t32 ^ 2 / 0.2e1 + (t35 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1) * t31 ^ 2) * t40, t34 ^ 2 * t49, t34 * t39 * t37, t34 * t43, t27, t37 * t43, qJD(3) ^ 2 / 0.2e1, t15 * qJD(3) - t25 * t44, -t16 * qJD(3) + t25 * t45, (-t15 * t34 + t16 * t37) * qJD(2), t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1, t23 ^ 2 / 0.2e1, -t23 * t21, -t23 * t44, t21 ^ 2 / 0.2e1, t21 * t44, t27, t9 * t21 - t6 * t44, t9 * t23 + t7 * t44, -t7 * t21 - t6 * t23, t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t14 ^ 2 / 0.2e1, -t14 * t12, t14 * t18, t12 ^ 2 / 0.2e1, -t12 * t18, t18 ^ 2 / 0.2e1, t1 * t18 + t3 * t12, t3 * t14 - t2 * t18, -t1 * t14 - t2 * t12, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t8;
