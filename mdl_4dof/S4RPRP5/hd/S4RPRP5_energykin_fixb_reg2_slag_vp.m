% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPRP5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP5_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:00
% EndTime: 2019-12-31 16:45:01
% DurationCPUTime: 0.10s
% Computational Cost: add. (106->29), mult. (338->69), div. (0->0), fcn. (192->4), ass. (0->29)
t24 = qJD(1) ^ 2;
t33 = t24 / 0.2e1;
t21 = sin(pkin(6));
t29 = qJD(1) * t21;
t30 = pkin(5) + qJ(2);
t13 = t30 * t29;
t22 = cos(pkin(6));
t28 = qJD(1) * t22;
t14 = t30 * t28;
t23 = sin(qJ(3));
t32 = cos(qJ(3));
t5 = -t23 * t13 + t32 * t14;
t10 = t23 * t29 - t32 * t28;
t12 = (t32 * t21 + t22 * t23) * qJD(1);
t31 = t12 * t10;
t27 = qJD(3) * t10;
t26 = t10 ^ 2 / 0.2e1;
t4 = -t32 * t13 - t23 * t14;
t15 = qJD(2) + (-pkin(2) * t22 - pkin(1)) * qJD(1);
t20 = qJD(3) ^ 2 / 0.2e1;
t19 = t22 ^ 2;
t18 = t21 ^ 2;
t17 = -qJD(1) * pkin(1) + qJD(2);
t7 = t12 * qJD(3);
t6 = t12 ^ 2 / 0.2e1;
t3 = qJD(3) * qJ(4) + t5;
t2 = -qJD(3) * pkin(3) + qJD(4) - t4;
t1 = t10 * pkin(3) - t12 * qJ(4) + t15;
t8 = [0, 0, 0, 0, 0, t33, 0, 0, 0, 0, t18 * t33, t21 * t24 * t22, 0, t19 * t33, 0, 0, -t17 * t28, t17 * t29, (t18 + t19) * t24 * qJ(2), t17 ^ 2 / 0.2e1 + (t19 / 0.2e1 + t18 / 0.2e1) * qJ(2) ^ 2 * t24, t6, -t31, t7, t26, -t27, t20, t4 * qJD(3) + t15 * t10, -t5 * qJD(3) + t15 * t12, -t5 * t10 - t4 * t12, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t6, t7, t31, t20, t27, t26, -t2 * qJD(3) + t1 * t10, -t3 * t10 + t2 * t12, t3 * qJD(3) - t1 * t12, t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1;];
T_reg = t8;
