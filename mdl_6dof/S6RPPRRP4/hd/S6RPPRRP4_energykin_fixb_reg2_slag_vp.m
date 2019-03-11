% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRRP4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP4_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:06:18
% EndTime: 2019-03-09 02:06:18
% DurationCPUTime: 0.12s
% Computational Cost: add. (327->49), mult. (598->102), div. (0->0), fcn. (286->6), ass. (0->38)
t41 = qJD(1) ^ 2;
t34 = t41 / 0.2e1;
t27 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t35 = sin(pkin(9));
t36 = cos(pkin(9));
t43 = qJ(2) * qJD(1);
t18 = t36 * t27 - t35 * t43;
t15 = qJD(1) * pkin(3) - t18;
t38 = sin(qJ(4));
t40 = cos(qJ(4));
t10 = (pkin(4) * t40 + pkin(8) * t38) * qJD(1) + t15;
t39 = cos(qJ(5));
t49 = sin(qJ(5));
t19 = t35 * t27 + t36 * t43;
t16 = -qJD(1) * pkin(7) + t19;
t12 = t38 * qJD(3) + t40 * t16;
t9 = qJD(4) * pkin(8) + t12;
t4 = t49 * t10 + t39 * t9;
t46 = qJD(1) * t38;
t21 = t39 * qJD(4) + t49 * t46;
t45 = t40 * qJD(1);
t28 = qJD(5) + t45;
t48 = t21 * t28;
t23 = -t49 * qJD(4) + t39 * t46;
t47 = t23 * t21;
t44 = t21 ^ 2 / 0.2e1;
t42 = qJD(1) * qJD(4);
t11 = t40 * qJD(3) - t38 * t16;
t3 = t39 * t10 - t49 * t9;
t8 = -qJD(4) * pkin(4) - t11;
t31 = -pkin(1) * qJD(1) + qJD(2);
t26 = t28 ^ 2 / 0.2e1;
t20 = t23 ^ 2 / 0.2e1;
t17 = t23 * t28;
t5 = -t21 * pkin(5) + t23 * qJ(6) + t8;
t2 = t28 * qJ(6) + t4;
t1 = -t28 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t34, 0, 0, 0, 0, 0, 0, 0, t34, 0, 0, -t31 * qJD(1), 0, t41 * qJ(2), qJ(2) ^ 2 * t34 + t31 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t34, -t18 * qJD(1), t19 * qJD(1), 0, t19 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t38 ^ 2 * t34, t38 * t41 * t40, -t38 * t42, t40 ^ 2 * t34, -t40 * t42, qJD(4) ^ 2 / 0.2e1, t11 * qJD(4) + t15 * t45, -t12 * qJD(4) - t15 * t46 (t11 * t38 - t12 * t40) * qJD(1), t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t20, -t47, -t17, t44, t48, t26, -t8 * t21 + t3 * t28, -t8 * t23 - t4 * t28, t4 * t21 + t3 * t23, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1, t20, -t17, t47, t26, -t48, t44, -t1 * t28 - t5 * t21, -t1 * t23 + t2 * t21, t2 * t28 + t5 * t23, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
