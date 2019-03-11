% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPPRR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR2_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:32:01
% EndTime: 2019-03-09 01:32:01
% DurationCPUTime: 0.16s
% Computational Cost: add. (319->48), mult. (662->114), div. (0->0), fcn. (379->8), ass. (0->37)
t35 = sin(pkin(10));
t37 = cos(pkin(10));
t40 = sin(qJ(5));
t41 = cos(qJ(5));
t21 = (t35 * t41 + t37 * t40) * qJD(1);
t42 = qJD(1) ^ 2;
t33 = t42 / 0.2e1;
t38 = cos(pkin(9));
t45 = -pkin(1) * t38 - pkin(2);
t24 = qJD(3) + (-qJ(4) + t45) * qJD(1);
t15 = -t35 * qJD(2) + t37 * t24;
t46 = qJD(1) * t37;
t10 = -pkin(7) * t46 + t15;
t16 = t37 * qJD(2) + t35 * t24;
t47 = qJD(1) * t35;
t11 = -pkin(7) * t47 + t16;
t6 = t40 * t10 + t41 * t11;
t49 = pkin(1) * t42;
t48 = cos(qJ(6));
t36 = sin(pkin(9));
t27 = (-pkin(1) * t36 - qJ(3)) * qJD(1);
t25 = qJD(4) - t27;
t20 = pkin(4) * t47 + t25;
t5 = t41 * t10 - t40 * t11;
t39 = sin(qJ(6));
t32 = qJD(2) ^ 2 / 0.2e1;
t26 = t45 * qJD(1) + qJD(3);
t23 = (-t35 * t40 + t37 * t41) * qJD(1);
t19 = qJD(6) + t21;
t14 = t39 * qJD(5) + t48 * t23;
t12 = -t48 * qJD(5) + t39 * t23;
t7 = t21 * pkin(5) - t23 * pkin(8) + t20;
t4 = qJD(5) * pkin(8) + t6;
t3 = -qJD(5) * pkin(5) - t5;
t2 = t39 * t7 + t48 * t4;
t1 = -t39 * t4 + t48 * t7;
t8 = [0, 0, 0, 0, 0, t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t38 * t49, -t36 * t49, 0, t32 + (t36 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t42, t33, 0, 0, 0, 0, 0, 0, t26 * qJD(1), -t27 * qJD(1), t32 + t27 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1, t37 ^ 2 * t33, -t37 * t42 * t35, 0, t35 ^ 2 * t33, 0, 0, t25 * t47, t25 * t46 (-t15 * t37 - t16 * t35) * qJD(1), t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1, t23 ^ 2 / 0.2e1, -t23 * t21, t23 * qJD(5), t21 ^ 2 / 0.2e1, -t21 * qJD(5), qJD(5) ^ 2 / 0.2e1, t5 * qJD(5) + t20 * t21, -t6 * qJD(5) + t20 * t23, -t6 * t21 - t5 * t23, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1, t14 ^ 2 / 0.2e1, -t14 * t12, t14 * t19, t12 ^ 2 / 0.2e1, -t12 * t19, t19 ^ 2 / 0.2e1, t1 * t19 + t3 * t12, t3 * t14 - t2 * t19, -t1 * t14 - t2 * t12, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
