% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPPRRP1
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
% Datum: 2019-03-09 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRRP1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP1_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:58:42
% EndTime: 2019-03-09 01:58:42
% DurationCPUTime: 0.15s
% Computational Cost: add. (436->56), mult. (1040->122), div. (0->0), fcn. (683->8), ass. (0->44)
t49 = qJD(1) ^ 2;
t42 = t49 / 0.2e1;
t56 = pkin(1) * t49;
t45 = cos(pkin(10));
t46 = cos(pkin(9));
t51 = -pkin(1) * t46 - pkin(2);
t31 = qJD(3) + (-pkin(3) * t45 + t51) * qJD(1);
t48 = sin(qJ(4));
t52 = qJD(1) * t45;
t43 = sin(pkin(10));
t53 = qJD(1) * t43;
t55 = cos(qJ(4));
t32 = t48 * t53 - t55 * t52;
t34 = (t55 * t43 + t45 * t48) * qJD(1);
t13 = t32 * pkin(4) - t34 * pkin(8) + t31;
t47 = sin(qJ(5));
t54 = cos(qJ(5));
t44 = sin(pkin(9));
t37 = (pkin(1) * t44 + qJ(3)) * qJD(1);
t40 = t45 * qJD(2);
t22 = t40 + (-pkin(7) * qJD(1) - t37) * t43;
t29 = t43 * qJD(2) + t45 * t37;
t23 = pkin(7) * t52 + t29;
t10 = t48 * t22 + t55 * t23;
t8 = qJD(4) * pkin(8) + t10;
t4 = t47 * t13 + t54 * t8;
t3 = t54 * t13 - t47 * t8;
t9 = t55 * t22 - t48 * t23;
t7 = -qJD(4) * pkin(4) - t9;
t36 = t51 * qJD(1) + qJD(3);
t30 = qJD(5) + t32;
t28 = -t43 * t37 + t40;
t27 = t30 ^ 2 / 0.2e1;
t26 = t47 * qJD(4) + t54 * t34;
t24 = -t54 * qJD(4) + t47 * t34;
t21 = t26 ^ 2 / 0.2e1;
t20 = t24 ^ 2 / 0.2e1;
t16 = t26 * t30;
t15 = t24 * t30;
t14 = t26 * t24;
t5 = t24 * pkin(5) + qJD(6) + t7;
t2 = -t24 * qJ(6) + t4;
t1 = t30 * pkin(5) - t26 * qJ(6) + t3;
t6 = [0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t46 * t56, -t44 * t56, 0, qJD(2) ^ 2 / 0.2e1 + (t44 ^ 2 / 0.2e1 + t46 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t49, t43 ^ 2 * t42, t43 * t49 * t45, 0, t45 ^ 2 * t42, 0, 0, -t36 * t52, t36 * t53 (-t28 * t43 + t29 * t45) * qJD(1), t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t36 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, t34 * qJD(4), t32 ^ 2 / 0.2e1, -t32 * qJD(4), qJD(4) ^ 2 / 0.2e1, t9 * qJD(4) + t31 * t32, -t10 * qJD(4) + t31 * t34, -t10 * t32 - t9 * t34, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1, t21, -t14, t16, t20, -t15, t27, t7 * t24 + t3 * t30, t7 * t26 - t4 * t30, -t4 * t24 - t3 * t26, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t21, -t14, t16, t20, -t15, t27, t1 * t30 + t5 * t24, -t2 * t30 + t5 * t26, -t1 * t26 - t2 * t24, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t6;
