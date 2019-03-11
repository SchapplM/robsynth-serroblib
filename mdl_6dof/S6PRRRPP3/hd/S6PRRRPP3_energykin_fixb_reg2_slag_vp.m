% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRPP3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:58:28
% EndTime: 2019-03-08 22:58:28
% DurationCPUTime: 0.14s
% Computational Cost: add. (317->54), mult. (737->115), div. (0->0), fcn. (481->8), ass. (0->46)
t42 = qJD(2) ^ 2;
t56 = t42 / 0.2e1;
t55 = cos(qJ(4));
t37 = sin(qJ(4));
t38 = sin(qJ(3));
t50 = qJD(2) * t38;
t23 = -t55 * qJD(3) + t37 * t50;
t25 = t37 * qJD(3) + t55 * t50;
t54 = t23 * t25;
t40 = cos(qJ(3));
t49 = t40 * qJD(2);
t31 = -qJD(4) + t49;
t17 = t25 * t31;
t18 = t31 * t23;
t53 = pkin(4) + qJ(6);
t39 = sin(qJ(2));
t35 = sin(pkin(6));
t52 = qJD(1) * t35;
t27 = qJD(2) * pkin(8) + t39 * t52;
t36 = cos(pkin(6));
t51 = qJD(1) * t36;
t15 = t40 * t27 + t38 * t51;
t12 = qJD(3) * pkin(9) + t15;
t41 = cos(qJ(2));
t47 = t41 * t52;
t16 = -t47 + (-pkin(3) * t40 - pkin(9) * t38 - pkin(2)) * qJD(2);
t7 = t55 * t12 + t37 * t16;
t19 = t23 ^ 2 / 0.2e1;
t20 = t25 ^ 2 / 0.2e1;
t48 = qJD(2) * qJD(3);
t46 = qJD(2) * t52;
t5 = t31 * qJ(5) - t7;
t6 = -t37 * t12 + t55 * t16;
t14 = -t38 * t27 + t40 * t51;
t45 = qJD(5) - t6;
t11 = -qJD(3) * pkin(3) - t14;
t44 = -t25 * qJ(5) + t11;
t43 = qJD(1) ^ 2;
t29 = t31 ^ 2 / 0.2e1;
t28 = -qJD(2) * pkin(2) - t47;
t8 = t23 * pkin(4) + t44;
t4 = t31 * pkin(4) + t45;
t3 = t53 * t23 + t44;
t2 = -t23 * pkin(5) + qJD(6) - t5;
t1 = t25 * pkin(5) + t53 * t31 + t45;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t43 / 0.2e1, 0, 0, 0, 0, 0, t56, t41 * t46, -t39 * t46, 0 (t36 ^ 2 / 0.2e1 + (t39 ^ 2 / 0.2e1 + t41 ^ 2 / 0.2e1) * t35 ^ 2) * t43, t38 ^ 2 * t56, t38 * t42 * t40, t38 * t48, t40 ^ 2 * t56, t40 * t48, qJD(3) ^ 2 / 0.2e1, t14 * qJD(3) - t28 * t49, -t15 * qJD(3) + t28 * t50 (-t14 * t38 + t15 * t40) * qJD(2), t15 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1, t20, -t54, -t17, t19, t18, t29, t11 * t23 - t6 * t31, t11 * t25 + t7 * t31, -t7 * t23 - t6 * t25, t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1, t29, t17, -t18, t20, -t54, t19, t5 * t23 + t4 * t25, -t8 * t23 - t4 * t31, -t8 * t25 + t5 * t31, t8 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1, t29, -t18, -t17, t19, t54, t20, t1 * t25 - t2 * t23, -t2 * t31 - t3 * t25, t1 * t31 + t3 * t23, t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1;];
T_reg  = t9;
