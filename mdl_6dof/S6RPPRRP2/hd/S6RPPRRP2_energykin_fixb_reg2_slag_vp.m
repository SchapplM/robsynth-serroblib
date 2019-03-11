% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRRP2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP2_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:01:16
% EndTime: 2019-03-09 02:01:16
% DurationCPUTime: 0.15s
% Computational Cost: add. (436->54), mult. (1034->122), div. (0->0), fcn. (677->8), ass. (0->44)
t45 = qJD(1) ^ 2;
t38 = t45 / 0.2e1;
t55 = pkin(1) * t45;
t41 = cos(pkin(10));
t42 = cos(pkin(9));
t47 = -pkin(1) * t42 - pkin(2);
t27 = qJD(3) + (-pkin(3) * t41 + t47) * qJD(1);
t44 = sin(qJ(4));
t49 = qJD(1) * t41;
t39 = sin(pkin(10));
t50 = qJD(1) * t39;
t54 = cos(qJ(4));
t28 = t44 * t50 - t54 * t49;
t30 = (t54 * t39 + t41 * t44) * qJD(1);
t12 = t28 * pkin(4) - t30 * pkin(8) + t27;
t43 = sin(qJ(5));
t53 = cos(qJ(5));
t40 = sin(pkin(9));
t33 = (pkin(1) * t40 + qJ(3)) * qJD(1);
t36 = t41 * qJD(2);
t18 = t36 + (-pkin(7) * qJD(1) - t33) * t39;
t25 = t39 * qJD(2) + t41 * t33;
t19 = pkin(7) * t49 + t25;
t10 = t44 * t18 + t54 * t19;
t8 = qJD(4) * pkin(8) + t10;
t4 = t43 * t12 + t53 * t8;
t20 = -t53 * qJD(4) + t43 * t30;
t22 = t43 * qJD(4) + t53 * t30;
t52 = t22 * t20;
t26 = qJD(5) + t28;
t51 = t26 * t20;
t48 = t20 ^ 2 / 0.2e1;
t9 = t54 * t18 - t44 * t19;
t3 = t53 * t12 - t43 * t8;
t7 = -qJD(4) * pkin(4) - t9;
t32 = t47 * qJD(1) + qJD(3);
t24 = -t39 * t33 + t36;
t23 = t26 ^ 2 / 0.2e1;
t17 = t22 ^ 2 / 0.2e1;
t13 = t22 * t26;
t5 = t20 * pkin(5) - t22 * qJ(6) + t7;
t2 = t26 * qJ(6) + t4;
t1 = -t26 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t42 * t55, -t40 * t55, 0, qJD(2) ^ 2 / 0.2e1 + (t40 ^ 2 / 0.2e1 + t42 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t45, t39 ^ 2 * t38, t39 * t45 * t41, 0, t41 ^ 2 * t38, 0, 0, -t32 * t49, t32 * t50 (-t24 * t39 + t25 * t41) * qJD(1), t25 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1, t30 ^ 2 / 0.2e1, -t30 * t28, t30 * qJD(4), t28 ^ 2 / 0.2e1, -t28 * qJD(4), qJD(4) ^ 2 / 0.2e1, t9 * qJD(4) + t27 * t28, -t10 * qJD(4) + t27 * t30, -t10 * t28 - t9 * t30, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1, t17, -t52, t13, t48, -t51, t23, t7 * t20 + t3 * t26, t7 * t22 - t4 * t26, -t4 * t20 - t3 * t22, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t17, t13, t52, t23, t51, t48, -t1 * t26 + t5 * t20, t1 * t22 - t2 * t20, t2 * t26 - t5 * t22, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
