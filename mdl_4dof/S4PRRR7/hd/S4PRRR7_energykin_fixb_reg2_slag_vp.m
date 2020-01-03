% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PRRR7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR7_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:47
% EndTime: 2019-12-31 16:36:47
% DurationCPUTime: 0.11s
% Computational Cost: add. (94->29), mult. (277->86), div. (0->0), fcn. (171->8), ass. (0->31)
t26 = qJD(2) ^ 2;
t36 = t26 / 0.2e1;
t23 = sin(qJ(2));
t19 = sin(pkin(4));
t34 = qJD(1) * t19;
t12 = qJD(2) * pkin(6) + t23 * t34;
t22 = sin(qJ(3));
t24 = cos(qJ(3));
t20 = cos(pkin(4));
t33 = qJD(1) * t20;
t6 = t24 * t12 + t22 * t33;
t35 = cos(qJ(4));
t32 = qJD(2) * t22;
t31 = t24 * qJD(2);
t30 = qJD(2) * qJD(3);
t25 = cos(qJ(2));
t29 = t25 * t34;
t28 = qJD(2) * t34;
t5 = -t22 * t12 + t24 * t33;
t27 = qJD(1) ^ 2;
t21 = sin(qJ(4));
t15 = -qJD(4) + t31;
t13 = -qJD(2) * pkin(2) - t29;
t11 = t21 * qJD(3) + t35 * t32;
t9 = -t35 * qJD(3) + t21 * t32;
t7 = -t29 + (-pkin(3) * t24 - pkin(7) * t22 - pkin(2)) * qJD(2);
t4 = qJD(3) * pkin(7) + t6;
t3 = -qJD(3) * pkin(3) - t5;
t2 = t21 * t7 + t35 * t4;
t1 = -t21 * t4 + t35 * t7;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t27 / 0.2e1, 0, 0, 0, 0, 0, t36, t25 * t28, -t23 * t28, 0, (t20 ^ 2 / 0.2e1 + (t23 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1) * t19 ^ 2) * t27, t22 ^ 2 * t36, t22 * t26 * t24, t22 * t30, t24 ^ 2 * t36, t24 * t30, qJD(3) ^ 2 / 0.2e1, t5 * qJD(3) - t13 * t31, -t6 * qJD(3) + t13 * t32, (-t22 * t5 + t24 * t6) * qJD(2), t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1, t11 ^ 2 / 0.2e1, -t11 * t9, -t11 * t15, t9 ^ 2 / 0.2e1, t9 * t15, t15 ^ 2 / 0.2e1, -t1 * t15 + t3 * t9, t3 * t11 + t2 * t15, -t1 * t11 - t2 * t9, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t8;
