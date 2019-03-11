% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PPRRRR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR1_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_energykin_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:01:28
% EndTime: 2019-03-08 19:01:29
% DurationCPUTime: 0.20s
% Computational Cost: add. (507->53), mult. (1251->134), div. (0->0), fcn. (1041->14), ass. (0->49)
t35 = cos(pkin(6)) * qJD(1) + qJD(2);
t40 = sin(pkin(7));
t43 = cos(pkin(7));
t42 = cos(pkin(13));
t41 = sin(pkin(6));
t60 = qJD(1) * t41;
t55 = t42 * t60;
t65 = t35 * t40 + t43 * t55;
t47 = sin(qJ(3));
t49 = cos(qJ(3));
t39 = sin(pkin(13));
t56 = t39 * t60;
t18 = -t47 * t56 + t65 * t49;
t50 = qJD(3) ^ 2;
t64 = t50 / 0.2e1;
t19 = t65 * t47 + t49 * t56;
t17 = qJD(3) * pkin(9) + t19;
t25 = t43 * t35 - t40 * t55;
t48 = cos(qJ(4));
t24 = t48 * t25;
t46 = sin(qJ(4));
t10 = qJD(4) * pkin(4) + t24 + (-pkin(10) * qJD(3) - t17) * t46;
t13 = t48 * t17 + t46 * t25;
t58 = qJD(3) * t48;
t11 = pkin(10) * t58 + t13;
t45 = sin(qJ(5));
t63 = cos(qJ(5));
t6 = t45 * t10 + t63 * t11;
t62 = cos(qJ(6));
t59 = qJD(3) * t46;
t57 = qJD(3) * qJD(4);
t27 = t45 * t59 - t63 * t58;
t5 = t63 * t10 - t45 * t11;
t14 = (-pkin(4) * t48 - pkin(3)) * qJD(3) - t18;
t51 = qJD(1) ^ 2;
t44 = sin(qJ(6));
t38 = qJD(4) + qJD(5);
t29 = (t45 * t48 + t63 * t46) * qJD(3);
t26 = qJD(6) + t27;
t22 = t62 * t29 + t44 * t38;
t20 = t44 * t29 - t62 * t38;
t16 = -qJD(3) * pkin(3) - t18;
t12 = -t46 * t17 + t24;
t8 = t27 * pkin(5) - t29 * pkin(11) + t14;
t4 = t38 * pkin(11) + t6;
t3 = -t38 * pkin(5) - t5;
t2 = t62 * t4 + t44 * t8;
t1 = -t44 * t4 + t62 * t8;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t51 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35 ^ 2 / 0.2e1 + (t39 ^ 2 / 0.2e1 + t42 ^ 2 / 0.2e1) * t51 * t41 ^ 2, 0, 0, 0, 0, 0, t64, t18 * qJD(3), -t19 * qJD(3), 0, t19 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1, t46 ^ 2 * t64, t46 * t50 * t48, t46 * t57, t48 ^ 2 * t64, t48 * t57, qJD(4) ^ 2 / 0.2e1, t12 * qJD(4) - t16 * t58, -t13 * qJD(4) + t16 * t59 (-t12 * t46 + t13 * t48) * qJD(3), t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1, t29 ^ 2 / 0.2e1, -t29 * t27, t29 * t38, t27 ^ 2 / 0.2e1, -t27 * t38, t38 ^ 2 / 0.2e1, t14 * t27 + t5 * t38, t14 * t29 - t6 * t38, -t6 * t27 - t5 * t29, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * t26, t20 ^ 2 / 0.2e1, -t20 * t26, t26 ^ 2 / 0.2e1, t1 * t26 + t3 * t20, -t2 * t26 + t3 * t22, -t1 * t22 - t2 * t20, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t7;
