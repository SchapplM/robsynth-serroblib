% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR7_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:39
% EndTime: 2019-12-31 20:15:39
% DurationCPUTime: 0.11s
% Computational Cost: add. (220->30), mult. (294->83), div. (0->0), fcn. (125->6), ass. (0->35)
t21 = sin(qJ(4));
t18 = t21 ^ 2;
t38 = t18 / 0.2e1;
t24 = cos(qJ(4));
t19 = t24 ^ 2;
t37 = t19 / 0.2e1;
t17 = qJD(1) + qJD(2);
t22 = sin(qJ(2));
t35 = pkin(1) * qJD(1);
t31 = t22 * t35;
t12 = t17 * qJ(3) + t31;
t36 = t12 * t17;
t34 = qJD(4) * t21;
t33 = qJD(4) * t24;
t32 = t12 ^ 2 / 0.2e1;
t25 = cos(qJ(2));
t30 = t25 * t35;
t28 = qJD(3) - t30;
t10 = (-pkin(2) - pkin(7)) * t17 + t28;
t29 = -pkin(8) * t17 + t10;
t26 = qJD(1) ^ 2;
t23 = cos(qJ(5));
t20 = sin(qJ(5));
t16 = qJD(4) + qJD(5);
t15 = t17 ^ 2;
t14 = t15 / 0.2e1;
t11 = -t17 * pkin(2) + t28;
t8 = t31 + (pkin(4) * t21 + qJ(3)) * t17;
t7 = (-t20 * t21 + t23 * t24) * t17;
t5 = (-t20 * t24 - t21 * t23) * t17;
t4 = t29 * t21;
t3 = qJD(4) * pkin(4) + t29 * t24;
t2 = t20 * t3 + t23 * t4;
t1 = -t20 * t4 + t23 * t3;
t6 = [0, 0, 0, 0, 0, t26 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t17 * t30, -t17 * t31, 0, (t22 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t26, t14, 0, 0, 0, 0, 0, 0, t11 * t17, t36, t32 + t11 ^ 2 / 0.2e1, t15 * t37, -t24 * t15 * t21, t17 * t33, t15 * t38, -t17 * t34, qJD(4) ^ 2 / 0.2e1, t10 * t33 + t21 * t36, -t10 * t34 + t24 * t36, (-t18 - t19) * t17 * t10, t32 + (t38 + t37) * t10 ^ 2, t7 ^ 2 / 0.2e1, t7 * t5, t7 * t16, t5 ^ 2 / 0.2e1, t5 * t16, t16 ^ 2 / 0.2e1, t1 * t16 - t8 * t5, -t2 * t16 + t8 * t7, -t1 * t7 + t2 * t5, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1;];
T_reg = t6;
