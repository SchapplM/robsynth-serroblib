% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRRP6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP6_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:11:16
% EndTime: 2019-03-09 02:11:17
% DurationCPUTime: 0.13s
% Computational Cost: add. (217->41), mult. (404->90), div. (0->0), fcn. (174->4), ass. (0->37)
t32 = qJD(1) ^ 2;
t25 = t32 / 0.2e1;
t22 = qJ(2) * qJD(1) + qJD(3);
t18 = -pkin(7) * qJD(1) + t22;
t30 = sin(qJ(4));
t10 = qJD(4) * pkin(8) + t30 * t18;
t29 = sin(qJ(5));
t43 = cos(qJ(5));
t31 = cos(qJ(4));
t40 = -pkin(1) - qJ(3);
t8 = -qJD(2) + (pkin(4) * t30 - pkin(8) * t31 - t40) * qJD(1);
t5 = t43 * t10 + t29 * t8;
t39 = qJD(1) * t31;
t13 = -t43 * qJD(4) + t29 * t39;
t15 = t29 * qJD(4) + t43 * t39;
t42 = t15 * t13;
t36 = t30 * qJD(1);
t21 = qJD(5) + t36;
t41 = t21 * t13;
t38 = qJD(4) * t18;
t19 = -t40 * qJD(1) - qJD(2);
t37 = t19 * qJD(1);
t35 = t13 ^ 2 / 0.2e1;
t34 = t19 ^ 2 / 0.2e1;
t33 = qJD(1) * qJD(4);
t11 = -qJD(4) * pkin(4) - t31 * t18;
t4 = -t29 * t10 + t43 * t8;
t28 = t31 ^ 2;
t27 = t30 ^ 2;
t23 = -qJD(1) * pkin(1) + qJD(2);
t17 = t21 ^ 2 / 0.2e1;
t12 = t15 ^ 2 / 0.2e1;
t7 = t15 * t21;
t3 = t13 * pkin(5) - t15 * qJ(6) + t11;
t2 = t21 * qJ(6) + t5;
t1 = -t21 * pkin(5) + qJD(6) - t4;
t6 = [0, 0, 0, 0, 0, t25, 0, 0, 0, 0, t25, 0, 0, 0, 0, 0, 0, t23 * qJD(1), t32 * qJ(2), qJ(2) ^ 2 * t25 + t23 ^ 2 / 0.2e1, t25, 0, 0, 0, 0, 0, 0, t22 * qJD(1), t37, t34 + t22 ^ 2 / 0.2e1, t28 * t25, -t31 * t32 * t30, t31 * t33, t27 * t25, -t30 * t33, qJD(4) ^ 2 / 0.2e1, t19 * t36 + t31 * t38, -t30 * t38 + t31 * t37 (-t27 - t28) * t18 * qJD(1), t34 + (t27 / 0.2e1 + t28 / 0.2e1) * t18 ^ 2, t12, -t42, t7, t35, -t41, t17, t11 * t13 + t4 * t21, t11 * t15 - t5 * t21, -t5 * t13 - t4 * t15, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1, t12, t7, t42, t17, t41, t35, -t1 * t21 + t3 * t13, t1 * t15 - t2 * t13, -t3 * t15 + t2 * t21, t2 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
