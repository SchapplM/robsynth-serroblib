% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPRR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR4_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:51:34
% EndTime: 2019-12-05 15:51:34
% DurationCPUTime: 0.16s
% Computational Cost: add. (169->36), mult. (436->99), div. (0->0), fcn. (292->10), ass. (0->36)
t36 = qJD(2) ^ 2;
t26 = t36 / 0.2e1;
t35 = cos(qJ(2));
t28 = sin(pkin(5));
t43 = qJD(1) * t28;
t17 = qJD(2) * pkin(2) + t35 * t43;
t27 = sin(pkin(10));
t29 = cos(pkin(10));
t33 = sin(qJ(2));
t39 = t33 * t43;
t12 = t27 * t17 + t29 * t39;
t10 = qJD(2) * pkin(7) + t12;
t30 = cos(pkin(5));
t21 = t30 * qJD(1) + qJD(3);
t32 = sin(qJ(4));
t34 = cos(qJ(4));
t6 = t34 * t10 + t32 * t21;
t44 = cos(qJ(5));
t42 = qJD(2) * t32;
t41 = t34 * qJD(2);
t40 = qJD(2) * qJD(4);
t38 = qJD(2) * t43;
t11 = t29 * t17 - t27 * t39;
t5 = -t32 * t10 + t34 * t21;
t37 = qJD(1) ^ 2;
t31 = sin(qJ(5));
t22 = -qJD(5) + t41;
t16 = t31 * qJD(4) + t44 * t42;
t14 = -t44 * qJD(4) + t31 * t42;
t9 = -qJD(2) * pkin(3) - t11;
t7 = (-pkin(4) * t34 - pkin(8) * t32 - pkin(3)) * qJD(2) - t11;
t4 = qJD(4) * pkin(8) + t6;
t3 = -qJD(4) * pkin(4) - t5;
t2 = t31 * t7 + t44 * t4;
t1 = -t31 * t4 + t44 * t7;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t37 / 0.2e1, 0, 0, 0, 0, 0, t26, t35 * t38, -t33 * t38, 0, (t30 ^ 2 / 0.2e1 + (t33 ^ 2 / 0.2e1 + t35 ^ 2 / 0.2e1) * t28 ^ 2) * t37, 0, 0, 0, 0, 0, t26, t11 * qJD(2), -t12 * qJD(2), 0, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t32 ^ 2 * t26, t32 * t36 * t34, t32 * t40, t34 ^ 2 * t26, t34 * t40, qJD(4) ^ 2 / 0.2e1, t5 * qJD(4) - t9 * t41, -t6 * qJD(4) + t9 * t42, (-t32 * t5 + t34 * t6) * qJD(2), t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t16 ^ 2 / 0.2e1, -t16 * t14, -t16 * t22, t14 ^ 2 / 0.2e1, t14 * t22, t22 ^ 2 / 0.2e1, -t1 * t22 + t3 * t14, t3 * t16 + t2 * t22, -t1 * t16 - t2 * t14, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t8;
