% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_energykin_fixb_reg2_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:03:21
% EndTime: 2019-12-05 17:03:22
% DurationCPUTime: 0.14s
% Computational Cost: add. (117->33), mult. (375->96), div. (0->0), fcn. (267->8), ass. (0->38)
t28 = qJD(2) ^ 2;
t44 = t28 / 0.2e1;
t29 = qJD(1) ^ 2;
t43 = t29 / 0.2e1;
t42 = cos(qJ(4));
t41 = cos(qJ(5));
t25 = sin(qJ(2));
t40 = t25 ^ 2 * t29;
t23 = sin(qJ(4));
t26 = cos(qJ(3));
t39 = t23 * t26;
t24 = sin(qJ(3));
t37 = qJD(1) * t25;
t30 = qJD(3) * pkin(2) - t24 * t37;
t6 = -t42 * t30 + t37 * t39;
t38 = t6 ^ 2 / 0.2e1;
t27 = cos(qJ(2));
t36 = qJD(2) * t27;
t35 = qJD(3) * t25;
t34 = qJD(1) * qJD(2);
t33 = qJD(2) * qJD(3);
t32 = t26 * t42;
t31 = t25 * t34;
t10 = (t23 * t24 - t32) * qJD(2);
t22 = sin(qJ(5));
t21 = t26 ^ 2;
t19 = t24 ^ 2;
t18 = qJD(3) + qJD(4);
t17 = t27 ^ 2 * t43;
t14 = -t26 * qJD(2) * pkin(2) - t27 * qJD(1);
t12 = (t42 * t24 + t39) * qJD(2);
t9 = qJD(5) + t10;
t8 = t23 * t30 + t32 * t37;
t5 = t41 * t12 + t22 * t18;
t3 = t22 * t12 - t41 * t18;
t2 = t22 * t14 + t41 * t8;
t1 = t41 * t14 - t22 * t8;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, 0, 0, 0, 0, t44, t27 * t34, -t31, 0, t40 / 0.2e1 + t17, t19 * t44, t24 * t28 * t26, t24 * t33, t21 * t44, t26 * t33, qJD(3) ^ 2 / 0.2e1, (-t24 * t35 + t26 * t36) * qJD(1), (-t24 * t36 - t26 * t35) * qJD(1), (t19 + t21) * t31, t17 + (t21 / 0.2e1 + t19 / 0.2e1) * t40, t12 ^ 2 / 0.2e1, -t12 * t10, t12 * t18, t10 ^ 2 / 0.2e1, -t10 * t18, t18 ^ 2 / 0.2e1, t14 * t10 - t6 * t18, t14 * t12 - t8 * t18, -t8 * t10 + t6 * t12, t8 ^ 2 / 0.2e1 + t38 + t14 ^ 2 / 0.2e1, t5 ^ 2 / 0.2e1, -t5 * t3, t5 * t9, t3 ^ 2 / 0.2e1, -t3 * t9, t9 ^ 2 / 0.2e1, t1 * t9 + t6 * t3, -t2 * t9 + t6 * t5, -t1 * t5 - t2 * t3, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t38;];
T_reg = t4;
