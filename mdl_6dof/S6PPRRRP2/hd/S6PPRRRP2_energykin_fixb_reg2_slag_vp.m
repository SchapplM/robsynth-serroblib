% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PPRRRP2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP2_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:58:07
% EndTime: 2019-03-08 18:58:08
% DurationCPUTime: 0.18s
% Computational Cost: add. (374->48), mult. (940->112), div. (0->0), fcn. (760->12), ass. (0->47)
t31 = cos(pkin(6)) * qJD(1) + qJD(2);
t37 = sin(pkin(7));
t40 = cos(pkin(7));
t39 = cos(pkin(12));
t38 = sin(pkin(6));
t57 = qJD(1) * t38;
t51 = t39 * t57;
t63 = t31 * t37 + t40 * t51;
t43 = sin(qJ(3));
t45 = cos(qJ(3));
t36 = sin(pkin(12));
t52 = t36 * t57;
t17 = -t43 * t52 + t45 * t63;
t46 = qJD(3) ^ 2;
t62 = t46 / 0.2e1;
t42 = sin(qJ(4));
t44 = cos(qJ(4));
t12 = (-pkin(4) * t44 - pkin(10) * t42 - pkin(3)) * qJD(3) - t17;
t41 = sin(qJ(5));
t61 = cos(qJ(5));
t18 = t63 * t43 + t45 * t52;
t16 = qJD(3) * pkin(9) + t18;
t20 = t40 * t31 - t37 * t51;
t10 = t44 * t16 + t42 * t20;
t8 = qJD(4) * pkin(10) + t10;
t4 = t41 * t12 + t61 * t8;
t56 = qJD(3) * t42;
t25 = -t61 * qJD(4) + t41 * t56;
t27 = t41 * qJD(4) + t61 * t56;
t60 = t27 * t25;
t55 = t44 * qJD(3);
t32 = -qJD(5) + t55;
t58 = t32 * t25;
t54 = t25 ^ 2 / 0.2e1;
t53 = qJD(3) * qJD(4);
t9 = -t42 * t16 + t44 * t20;
t7 = -qJD(4) * pkin(4) - t9;
t3 = t61 * t12 - t41 * t8;
t47 = qJD(1) ^ 2;
t30 = t32 ^ 2 / 0.2e1;
t22 = t27 ^ 2 / 0.2e1;
t21 = t27 * t32;
t15 = -qJD(3) * pkin(3) - t17;
t5 = t25 * pkin(5) - t27 * qJ(6) + t7;
t2 = -t32 * qJ(6) + t4;
t1 = t32 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t47 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31 ^ 2 / 0.2e1 + (t36 ^ 2 / 0.2e1 + t39 ^ 2 / 0.2e1) * t47 * t38 ^ 2, 0, 0, 0, 0, 0, t62, t17 * qJD(3), -t18 * qJD(3), 0, t18 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1, t42 ^ 2 * t62, t42 * t46 * t44, t42 * t53, t44 ^ 2 * t62, t44 * t53, qJD(4) ^ 2 / 0.2e1, t9 * qJD(4) - t15 * t55, -t10 * qJD(4) + t15 * t56 (t10 * t44 - t42 * t9) * qJD(3), t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t22, -t60, -t21, t54, t58, t30, t7 * t25 - t3 * t32, t7 * t27 + t4 * t32, -t4 * t25 - t3 * t27, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t22, -t21, t60, t30, -t58, t54, t1 * t32 + t5 * t25, t1 * t27 - t2 * t25, -t2 * t32 - t5 * t27, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
