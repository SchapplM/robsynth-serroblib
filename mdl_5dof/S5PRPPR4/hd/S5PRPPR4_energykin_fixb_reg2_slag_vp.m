% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPPR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR4_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:54
% EndTime: 2019-12-31 17:36:54
% DurationCPUTime: 0.10s
% Computational Cost: add. (95->31), mult. (252->68), div. (0->0), fcn. (135->4), ass. (0->29)
t26 = qJD(2) ^ 2;
t32 = t26 / 0.2e1;
t22 = sin(pkin(8));
t23 = cos(pkin(8));
t29 = qJ(3) * qJD(2);
t14 = t22 * qJD(1) + t23 * t29;
t31 = qJD(2) * t22;
t30 = qJD(2) * t23;
t28 = t22 * t26 * t23;
t27 = qJ(4) * t22 + pkin(2);
t13 = t23 * qJD(1) - t22 * t29;
t12 = qJD(4) - t13;
t25 = cos(qJ(5));
t24 = sin(qJ(5));
t21 = qJD(1) ^ 2 / 0.2e1;
t19 = -qJD(2) * pkin(2) + qJD(3);
t16 = t23 ^ 2 * t32;
t15 = t22 ^ 2 * t32;
t11 = t14 ^ 2 / 0.2e1;
t10 = (t22 * t25 - t23 * t24) * qJD(2);
t8 = (-t22 * t24 - t23 * t25) * qJD(2);
t7 = t14 * t30;
t6 = qJD(3) + (-pkin(3) * t23 - t27) * qJD(2);
t5 = -pkin(6) * t30 + t14;
t4 = -pkin(6) * t31 + t12;
t3 = -qJD(3) + ((pkin(3) + pkin(4)) * t23 + t27) * qJD(2);
t2 = t24 * t4 + t25 * t5;
t1 = -t24 * t5 + t25 * t4;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, 0, 0, 0, t32, 0, 0, 0, t21, t15, t28, 0, t16, 0, 0, -t19 * t30, t19 * t31, -t13 * t31 + t7, t11 + t13 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t15, 0, -t28, 0, 0, t16, -t6 * t30, t12 * t31 + t7, -t6 * t31, t11 + t6 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, t10 * t8, t10 * qJD(5), t8 ^ 2 / 0.2e1, t8 * qJD(5), qJD(5) ^ 2 / 0.2e1, t1 * qJD(5) - t3 * t8, -t2 * qJD(5) + t3 * t10, -t1 * t10 + t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t9;
