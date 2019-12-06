% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRP1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:36:24
% EndTime: 2019-12-05 17:36:24
% DurationCPUTime: 0.16s
% Computational Cost: add. (159->38), mult. (409->87), div. (0->0), fcn. (210->6), ass. (0->39)
t30 = qJD(1) ^ 2;
t23 = t30 / 0.2e1;
t25 = sin(pkin(7));
t15 = (pkin(1) * t25 + qJ(3)) * qJD(1);
t24 = sin(pkin(8));
t26 = cos(pkin(8));
t12 = t24 * qJD(2) + t26 * t15;
t28 = sin(qJ(4));
t29 = cos(qJ(4));
t27 = cos(pkin(7));
t35 = -pkin(1) * t27 - pkin(2);
t8 = qJD(3) + (-pkin(3) * t26 - pkin(6) * t24 + t35) * qJD(1);
t4 = t29 * t12 + t28 * t8;
t45 = pkin(1) * t30;
t44 = t24 ^ 2 * t30;
t43 = qJD(1) * t24;
t42 = qJD(1) * t28;
t41 = t26 * qJD(1);
t21 = t26 * qJD(2);
t10 = t24 * t15 - t21;
t40 = t10 ^ 2 / 0.2e1;
t39 = t24 * t42;
t38 = t29 * t43;
t37 = t10 * t43;
t36 = t44 / 0.2e1;
t34 = qJ(5) * t43;
t3 = -t28 * t12 + t29 * t8;
t33 = t29 * t28 * t44;
t19 = -qJD(4) + t41;
t32 = t19 * t39;
t18 = t29 ^ 2 * t36;
t17 = t28 ^ 2 * t36;
t16 = t19 ^ 2 / 0.2e1;
t14 = t35 * qJD(1) + qJD(3);
t13 = t19 * t38;
t5 = qJD(5) - t21 + (pkin(4) * t42 + t15) * t24;
t2 = -t28 * t34 + t4;
t1 = -t19 * pkin(4) - t29 * t34 + t3;
t6 = [0, 0, 0, 0, 0, t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t27 * t45, -t25 * t45, 0, qJD(2) ^ 2 / 0.2e1 + (t25 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t30, t36, t24 * t30 * t26, 0, t26 ^ 2 * t23, 0, 0, -t14 * t41, t14 * t43, (t10 * t24 + t12 * t26) * qJD(1), t12 ^ 2 / 0.2e1 + t40 + t14 ^ 2 / 0.2e1, t18, -t33, -t13, t17, t32, t16, -t3 * t19 + t28 * t37, t4 * t19 + t29 * t37, (-t28 * t4 - t29 * t3) * t43, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t40, t18, -t33, -t13, t17, t32, t16, -t1 * t19 + t5 * t39, t2 * t19 + t5 * t38, (-t1 * t29 - t2 * t28) * t43, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t6;
