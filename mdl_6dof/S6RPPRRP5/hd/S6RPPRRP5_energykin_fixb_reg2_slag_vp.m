% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPPRRP5
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
% Datum: 2019-03-09 02:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRRP5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP5_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:08:46
% EndTime: 2019-03-09 02:08:46
% DurationCPUTime: 0.13s
% Computational Cost: add. (217->43), mult. (410->90), div. (0->0), fcn. (180->4), ass. (0->37)
t36 = qJD(1) ^ 2;
t29 = t36 / 0.2e1;
t34 = sin(qJ(4));
t35 = cos(qJ(4));
t43 = -pkin(1) - qJ(3);
t11 = -qJD(2) + (pkin(4) * t34 - pkin(8) * t35 - t43) * qJD(1);
t26 = qJ(2) * qJD(1) + qJD(3);
t22 = -pkin(7) * qJD(1) + t26;
t13 = qJD(4) * pkin(8) + t34 * t22;
t33 = sin(qJ(5));
t44 = cos(qJ(5));
t4 = t33 * t11 + t44 * t13;
t42 = qJD(1) * t35;
t41 = qJD(4) * t22;
t23 = -qJD(1) * t43 - qJD(2);
t40 = t23 * qJD(1);
t39 = t34 * qJD(1);
t38 = t23 ^ 2 / 0.2e1;
t37 = qJD(1) * qJD(4);
t3 = t44 * t11 - t33 * t13;
t14 = -qJD(4) * pkin(4) - t35 * t22;
t32 = t35 ^ 2;
t31 = t34 ^ 2;
t27 = -qJD(1) * pkin(1) + qJD(2);
t25 = qJD(5) + t39;
t21 = t25 ^ 2 / 0.2e1;
t19 = t33 * qJD(4) + t42 * t44;
t17 = -qJD(4) * t44 + t33 * t42;
t16 = t19 ^ 2 / 0.2e1;
t15 = t17 ^ 2 / 0.2e1;
t10 = t19 * t25;
t9 = t17 * t25;
t6 = t19 * t17;
t5 = t17 * pkin(5) + qJD(6) + t14;
t2 = -t17 * qJ(6) + t4;
t1 = t25 * pkin(5) - t19 * qJ(6) + t3;
t7 = [0, 0, 0, 0, 0, t29, 0, 0, 0, 0, t29, 0, 0, 0, 0, 0, 0, t27 * qJD(1), t36 * qJ(2), qJ(2) ^ 2 * t29 + t27 ^ 2 / 0.2e1, t29, 0, 0, 0, 0, 0, 0, t26 * qJD(1), t40, t38 + t26 ^ 2 / 0.2e1, t32 * t29, -t35 * t36 * t34, t35 * t37, t31 * t29, -t34 * t37, qJD(4) ^ 2 / 0.2e1, t23 * t39 + t35 * t41, -t34 * t41 + t35 * t40 (-t31 - t32) * t22 * qJD(1), t38 + (t31 / 0.2e1 + t32 / 0.2e1) * t22 ^ 2, t16, -t6, t10, t15, -t9, t21, t14 * t17 + t3 * t25, t14 * t19 - t4 * t25, -t4 * t17 - t3 * t19, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t16, -t6, t10, t15, -t9, t21, t1 * t25 + t5 * t17, t5 * t19 - t2 * t25, -t1 * t19 - t2 * t17, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t7;
