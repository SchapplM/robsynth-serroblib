% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PPRRPR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_energykin_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:47:27
% EndTime: 2019-03-08 18:47:27
% DurationCPUTime: 0.20s
% Computational Cost: add. (502->52), mult. (1255->129), div. (0->0), fcn. (1035->14), ass. (0->48)
t35 = cos(pkin(6)) * qJD(1) + qJD(2);
t43 = sin(pkin(7));
t46 = cos(pkin(7));
t45 = cos(pkin(12));
t44 = sin(pkin(6));
t62 = qJD(1) * t44;
t57 = t45 * t62;
t67 = t35 * t43 + t46 * t57;
t49 = sin(qJ(3));
t51 = cos(qJ(3));
t42 = sin(pkin(12));
t58 = t42 * t62;
t20 = -t49 * t58 + t67 * t51;
t52 = qJD(3) ^ 2;
t66 = t52 / 0.2e1;
t21 = t67 * t49 + t51 * t58;
t19 = qJD(3) * pkin(9) + t21;
t26 = t46 * t35 - t43 * t57;
t48 = sin(qJ(4));
t50 = cos(qJ(4));
t12 = t50 * t19 + t48 * t26;
t10 = qJD(4) * qJ(5) + t12;
t15 = (-pkin(4) * t50 - qJ(5) * t48 - pkin(3)) * qJD(3) - t20;
t41 = sin(pkin(13));
t63 = cos(pkin(13));
t6 = t63 * t10 + t41 * t15;
t65 = cos(qJ(6));
t61 = qJD(3) * t48;
t60 = t50 * qJD(3);
t59 = qJD(3) * qJD(4);
t5 = -t41 * t10 + t63 * t15;
t11 = -t48 * t19 + t50 * t26;
t9 = -qJD(4) * pkin(4) + qJD(5) - t11;
t53 = qJD(1) ^ 2;
t47 = sin(qJ(6));
t38 = t50 ^ 2 * t66;
t36 = -qJD(6) + t60;
t31 = t41 * qJD(4) + t63 * t61;
t29 = -t63 * qJD(4) + t41 * t61;
t24 = -t47 * t29 + t65 * t31;
t22 = t65 * t29 + t47 * t31;
t18 = -qJD(3) * pkin(3) - t20;
t7 = t29 * pkin(5) + t9;
t4 = -t29 * pkin(10) + t6;
t3 = -pkin(5) * t60 - t31 * pkin(10) + t5;
t2 = t47 * t3 + t65 * t4;
t1 = t65 * t3 - t47 * t4;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t53 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35 ^ 2 / 0.2e1 + (t42 ^ 2 / 0.2e1 + t45 ^ 2 / 0.2e1) * t53 * t44 ^ 2, 0, 0, 0, 0, 0, t66, t20 * qJD(3), -t21 * qJD(3), 0, t21 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1, t48 ^ 2 * t66, t50 * t52 * t48, t48 * t59, t38, t50 * t59, qJD(4) ^ 2 / 0.2e1, t11 * qJD(4) - t18 * t60, -t12 * qJD(4) + t18 * t61 (-t11 * t48 + t12 * t50) * qJD(3), t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t31 ^ 2 / 0.2e1, -t31 * t29, -t31 * t60, t29 ^ 2 / 0.2e1, t29 * t60, t38, t9 * t29 - t5 * t60, t9 * t31 + t6 * t60, -t6 * t29 - t5 * t31, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t24 ^ 2 / 0.2e1, -t24 * t22, -t24 * t36, t22 ^ 2 / 0.2e1, t22 * t36, t36 ^ 2 / 0.2e1, -t1 * t36 + t7 * t22, t2 * t36 + t7 * t24, -t1 * t24 - t2 * t22, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
