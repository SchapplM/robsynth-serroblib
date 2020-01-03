% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRRR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR2_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:17
% EndTime: 2019-12-31 17:23:18
% DurationCPUTime: 0.10s
% Computational Cost: add. (149->26), mult. (256->80), div. (0->0), fcn. (121->6), ass. (0->32)
t19 = sin(qJ(3));
t16 = t19 ^ 2;
t35 = t16 / 0.2e1;
t21 = cos(qJ(3));
t17 = t21 ^ 2;
t34 = t17 / 0.2e1;
t33 = cos(qJ(4));
t15 = qJD(1) + qJD(2);
t32 = t15 * t19;
t31 = t15 * t21;
t30 = pkin(1) * qJD(1);
t29 = qJD(3) * t19;
t28 = qJD(3) * t21;
t20 = sin(qJ(2));
t27 = t20 * t30;
t22 = cos(qJ(2));
t26 = t22 * t30;
t10 = t15 * pkin(6) + t27;
t25 = pkin(7) * t15 + t10;
t23 = qJD(1) ^ 2;
t18 = sin(qJ(4));
t14 = qJD(3) + qJD(4);
t13 = t15 ^ 2;
t11 = -t15 * pkin(2) - t26;
t8 = -t26 + (-pkin(3) * t21 - pkin(2)) * t15;
t7 = (t18 * t21 + t33 * t19) * t15;
t5 = t18 * t32 - t33 * t31;
t4 = t25 * t21;
t3 = qJD(3) * pkin(3) - t25 * t19;
t2 = t18 * t3 + t33 * t4;
t1 = -t18 * t4 + t33 * t3;
t6 = [0, 0, 0, 0, 0, t23 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13 / 0.2e1, t15 * t26, -t15 * t27, 0, (t20 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t23, t13 * t35, t19 * t13 * t21, t15 * t29, t13 * t34, t15 * t28, qJD(3) ^ 2 / 0.2e1, -t10 * t29 - t11 * t31, -t10 * t28 + t11 * t32, (t16 + t17) * t15 * t10, t11 ^ 2 / 0.2e1 + (t34 + t35) * t10 ^ 2, t7 ^ 2 / 0.2e1, -t7 * t5, t7 * t14, t5 ^ 2 / 0.2e1, -t5 * t14, t14 ^ 2 / 0.2e1, t1 * t14 + t8 * t5, -t2 * t14 + t8 * t7, -t1 * t7 - t2 * t5, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1;];
T_reg = t6;
