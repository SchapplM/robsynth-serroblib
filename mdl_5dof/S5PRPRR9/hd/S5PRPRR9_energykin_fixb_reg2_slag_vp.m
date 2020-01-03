% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPRR9_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR9_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:44
% EndTime: 2019-12-31 17:39:44
% DurationCPUTime: 0.09s
% Computational Cost: add. (88->21), mult. (150->55), div. (0->0), fcn. (48->4), ass. (0->22)
t12 = -qJD(2) + qJD(4);
t11 = t12 ^ 2;
t23 = t11 / 0.2e1;
t19 = qJD(2) ^ 2;
t13 = t19 / 0.2e1;
t16 = sin(qJ(4));
t18 = cos(qJ(4));
t20 = qJ(3) * qJD(2);
t8 = qJD(3) + (-pkin(2) - pkin(3)) * qJD(2);
t6 = t16 * t8 + t18 * t20;
t5 = -t16 * t20 + t18 * t8;
t3 = -t12 * pkin(4) - t5;
t22 = t12 * t3;
t21 = qJD(5) * t12;
t17 = cos(qJ(5));
t15 = sin(qJ(5));
t14 = qJD(1) ^ 2 / 0.2e1;
t10 = -qJD(2) * pkin(2) + qJD(3);
t4 = t12 * pkin(7) + t6;
t2 = -t15 * qJD(1) + t17 * t4;
t1 = -t17 * qJD(1) - t15 * t4;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, 0, 0, 0, t13, 0, 0, 0, t14, 0, 0, 0, t13, 0, 0, -t10 * qJD(2), 0, t19 * qJ(3), qJ(3) ^ 2 * t13 + t14 + t10 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t23, t5 * t12, -t6 * t12, 0, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t14, t15 ^ 2 * t23, t15 * t11 * t17, t15 * t21, t17 ^ 2 * t23, t17 * t21, qJD(5) ^ 2 / 0.2e1, t1 * qJD(5) - t17 * t22, -t2 * qJD(5) + t15 * t22, (-t1 * t15 + t17 * t2) * t12, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t7;
