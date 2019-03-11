% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRRP3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:03:47
% EndTime: 2019-03-09 02:03:47
% DurationCPUTime: 0.14s
% Computational Cost: add. (251->46), mult. (506->98), div. (0->0), fcn. (247->6), ass. (0->40)
t35 = qJD(1) ^ 2;
t29 = t35 / 0.2e1;
t33 = sin(qJ(4));
t34 = cos(qJ(4));
t30 = sin(pkin(9));
t37 = -pkin(1) * t30 - qJ(3);
t12 = (pkin(4) * t33 - pkin(8) * t34 - t37) * qJD(1);
t32 = sin(qJ(5));
t47 = cos(qJ(5));
t31 = cos(pkin(9));
t38 = -pkin(1) * t31 - pkin(2);
t16 = qJD(3) + (-pkin(7) + t38) * qJD(1);
t11 = t34 * qJD(2) + t33 * t16;
t8 = qJD(4) * pkin(8) + t11;
t5 = t32 * t12 + t47 * t8;
t48 = pkin(1) * t35;
t44 = qJD(1) * t34;
t17 = -t47 * qJD(4) + t32 * t44;
t19 = t32 * qJD(4) + t47 * t44;
t46 = t19 * t17;
t42 = t33 * qJD(1);
t24 = qJD(5) + t42;
t45 = t24 * t17;
t21 = t37 * qJD(1);
t43 = t21 * qJD(1);
t41 = t17 ^ 2 / 0.2e1;
t40 = t21 ^ 2 / 0.2e1;
t39 = qJD(1) * qJD(4);
t10 = -t33 * qJD(2) + t34 * t16;
t7 = -qJD(4) * pkin(4) - t10;
t4 = t47 * t12 - t32 * t8;
t28 = qJD(2) ^ 2 / 0.2e1;
t23 = t24 ^ 2 / 0.2e1;
t20 = t38 * qJD(1) + qJD(3);
t15 = t19 ^ 2 / 0.2e1;
t13 = t19 * t24;
t3 = t17 * pkin(5) - t19 * qJ(6) + t7;
t2 = t24 * qJ(6) + t5;
t1 = -t24 * pkin(5) + qJD(6) - t4;
t6 = [0, 0, 0, 0, 0, t29, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t31 * t48, -t30 * t48, 0, t28 + (t30 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t35, t29, 0, 0, 0, 0, 0, 0, t20 * qJD(1), -t43, t28 + t40 + t20 ^ 2 / 0.2e1, t34 ^ 2 * t29, -t34 * t35 * t33, t34 * t39, t33 ^ 2 * t29, -t33 * t39, qJD(4) ^ 2 / 0.2e1, t10 * qJD(4) - t21 * t42, -t11 * qJD(4) - t34 * t43 (-t10 * t34 - t11 * t33) * qJD(1), t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t40, t15, -t46, t13, t41, -t45, t23, t7 * t17 + t4 * t24, t7 * t19 - t5 * t24, -t5 * t17 - t4 * t19, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t15, t13, t46, t23, t45, t41, -t1 * t24 + t3 * t17, t1 * t19 - t2 * t17, -t3 * t19 + t2 * t24, t2 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
