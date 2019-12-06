% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRR8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR8_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:16:39
% EndTime: 2019-12-05 17:16:39
% DurationCPUTime: 0.18s
% Computational Cost: add. (274->46), mult. (677->121), div. (0->0), fcn. (472->10), ass. (0->42)
t39 = qJD(2) ^ 2;
t50 = t39 / 0.2e1;
t36 = sin(qJ(2));
t31 = sin(pkin(5));
t47 = qJD(1) * t31;
t23 = qJD(2) * pkin(7) + t36 * t47;
t37 = cos(qJ(3));
t32 = cos(pkin(5));
t46 = qJD(1) * t32;
t26 = t37 * t46;
t35 = sin(qJ(3));
t10 = qJD(3) * pkin(3) + t26 + (-pkin(8) * qJD(2) - t23) * t35;
t16 = t37 * t23 + t35 * t46;
t44 = qJD(2) * t37;
t11 = pkin(8) * t44 + t16;
t34 = sin(qJ(4));
t49 = cos(qJ(4));
t6 = t34 * t10 + t49 * t11;
t48 = cos(qJ(5));
t45 = qJD(2) * t35;
t43 = qJD(2) * qJD(3);
t38 = cos(qJ(2));
t42 = t38 * t47;
t41 = qJD(2) * t47;
t19 = t34 * t45 - t49 * t44;
t5 = t49 * t10 - t34 * t11;
t17 = -t42 + (-pkin(3) * t37 - pkin(2)) * qJD(2);
t40 = qJD(1) ^ 2;
t33 = sin(qJ(5));
t30 = qJD(3) + qJD(4);
t24 = -qJD(2) * pkin(2) - t42;
t21 = (t34 * t37 + t49 * t35) * qJD(2);
t18 = qJD(5) + t19;
t15 = -t35 * t23 + t26;
t14 = t48 * t21 + t33 * t30;
t12 = t33 * t21 - t48 * t30;
t7 = t19 * pkin(4) - t21 * pkin(9) + t17;
t4 = t30 * pkin(9) + t6;
t3 = -t30 * pkin(4) - t5;
t2 = t33 * t7 + t48 * t4;
t1 = -t33 * t4 + t48 * t7;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t40 / 0.2e1, 0, 0, 0, 0, 0, t50, t38 * t41, -t36 * t41, 0, (t32 ^ 2 / 0.2e1 + (t36 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1) * t31 ^ 2) * t40, t35 ^ 2 * t50, t35 * t39 * t37, t35 * t43, t37 ^ 2 * t50, t37 * t43, qJD(3) ^ 2 / 0.2e1, t15 * qJD(3) - t24 * t44, -t16 * qJD(3) + t24 * t45, (-t15 * t35 + t16 * t37) * qJD(2), t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1, t21 ^ 2 / 0.2e1, -t21 * t19, t21 * t30, t19 ^ 2 / 0.2e1, -t19 * t30, t30 ^ 2 / 0.2e1, t17 * t19 + t5 * t30, t17 * t21 - t6 * t30, -t6 * t19 - t5 * t21, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t14 ^ 2 / 0.2e1, -t14 * t12, t14 * t18, t12 ^ 2 / 0.2e1, -t12 * t18, t18 ^ 2 / 0.2e1, t1 * t18 + t3 * t12, t3 * t14 - t2 * t18, -t1 * t14 - t2 * t12, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t8;
