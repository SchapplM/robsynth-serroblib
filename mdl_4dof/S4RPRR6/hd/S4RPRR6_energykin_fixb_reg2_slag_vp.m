% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPRR6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR6_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:40
% EndTime: 2019-12-31 16:52:40
% DurationCPUTime: 0.12s
% Computational Cost: add. (170->35), mult. (515->89), div. (0->0), fcn. (336->6), ass. (0->30)
t30 = qJD(1) ^ 2;
t37 = t30 / 0.2e1;
t36 = cos(qJ(3));
t35 = cos(qJ(4));
t34 = pkin(5) + qJ(2);
t26 = sin(pkin(7));
t33 = qJD(1) * t26;
t18 = t34 * t33;
t27 = cos(pkin(7));
t32 = qJD(1) * t27;
t19 = t34 * t32;
t29 = sin(qJ(3));
t9 = -t29 * t18 + t36 * t19;
t8 = -t36 * t18 - t29 * t19;
t20 = qJD(2) + (-pkin(2) * t27 - pkin(1)) * qJD(1);
t28 = sin(qJ(4));
t25 = qJD(3) + qJD(4);
t24 = t27 ^ 2;
t23 = t26 ^ 2;
t22 = -qJD(1) * pkin(1) + qJD(2);
t17 = (t36 * t26 + t27 * t29) * qJD(1);
t15 = t29 * t33 - t36 * t32;
t10 = t15 * pkin(3) + t20;
t7 = -t28 * t15 + t35 * t17;
t5 = t35 * t15 + t28 * t17;
t4 = -t15 * pkin(6) + t9;
t3 = qJD(3) * pkin(3) - t17 * pkin(6) + t8;
t2 = t28 * t3 + t35 * t4;
t1 = -t28 * t4 + t35 * t3;
t6 = [0, 0, 0, 0, 0, t37, 0, 0, 0, 0, t23 * t37, t26 * t30 * t27, 0, t24 * t37, 0, 0, -t22 * t32, t22 * t33, (t23 + t24) * t30 * qJ(2), t22 ^ 2 / 0.2e1 + (t24 / 0.2e1 + t23 / 0.2e1) * qJ(2) ^ 2 * t30, t17 ^ 2 / 0.2e1, -t17 * t15, t17 * qJD(3), t15 ^ 2 / 0.2e1, -t15 * qJD(3), qJD(3) ^ 2 / 0.2e1, t8 * qJD(3) + t20 * t15, -t9 * qJD(3) + t20 * t17, -t9 * t15 - t8 * t17, t9 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1, t7 ^ 2 / 0.2e1, -t7 * t5, t7 * t25, t5 ^ 2 / 0.2e1, -t5 * t25, t25 ^ 2 / 0.2e1, t1 * t25 + t10 * t5, t10 * t7 - t2 * t25, -t1 * t7 - t2 * t5, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg = t6;
