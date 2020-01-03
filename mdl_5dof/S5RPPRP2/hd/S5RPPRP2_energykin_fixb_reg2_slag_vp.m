% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPPRP2
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
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRP2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP2_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:23
% EndTime: 2019-12-31 17:49:23
% DurationCPUTime: 0.14s
% Computational Cost: add. (174->38), mult. (458->89), div. (0->0), fcn. (263->6), ass. (0->34)
t31 = qJD(1) ^ 2;
t25 = t31 / 0.2e1;
t30 = sin(qJ(4));
t39 = cos(qJ(4));
t27 = sin(pkin(7));
t20 = (pkin(1) * t27 + qJ(3)) * qJD(1);
t28 = cos(pkin(8));
t23 = t28 * qJD(2);
t26 = sin(pkin(8));
t8 = t23 + (-pkin(6) * qJD(1) - t20) * t26;
t11 = t26 * qJD(2) + t28 * t20;
t36 = qJD(1) * t28;
t9 = pkin(6) * t36 + t11;
t4 = t30 * t8 + t39 * t9;
t40 = pkin(1) * t31;
t37 = qJD(1) * t26;
t15 = t30 * t37 - t39 * t36;
t17 = (t39 * t26 + t28 * t30) * qJD(1);
t38 = t17 * t15;
t35 = qJD(4) * t15;
t34 = t15 ^ 2 / 0.2e1;
t29 = cos(pkin(7));
t33 = -pkin(1) * t29 - pkin(2);
t3 = -t30 * t9 + t39 * t8;
t14 = qJD(3) + (-pkin(3) * t28 + t33) * qJD(1);
t24 = qJD(4) ^ 2 / 0.2e1;
t19 = t33 * qJD(1) + qJD(3);
t13 = t17 * qJD(4);
t12 = t17 ^ 2 / 0.2e1;
t10 = -t26 * t20 + t23;
t5 = t15 * pkin(4) - t17 * qJ(5) + t14;
t2 = qJD(4) * qJ(5) + t4;
t1 = -qJD(4) * pkin(4) + qJD(5) - t3;
t6 = [0, 0, 0, 0, 0, t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t29 * t40, -t27 * t40, 0, qJD(2) ^ 2 / 0.2e1 + (t27 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t31, t26 ^ 2 * t25, t26 * t31 * t28, 0, t28 ^ 2 * t25, 0, 0, -t19 * t36, t19 * t37, (-t10 * t26 + t11 * t28) * qJD(1), t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t12, -t38, t13, t34, -t35, t24, t3 * qJD(4) + t14 * t15, -t4 * qJD(4) + t14 * t17, -t4 * t15 - t3 * t17, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t12, t13, t38, t24, t35, t34, -t1 * qJD(4) + t5 * t15, t1 * t17 - t2 * t15, t2 * qJD(4) - t5 * t17, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg = t6;
