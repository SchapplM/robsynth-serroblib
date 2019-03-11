% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRRP5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:16:42
% EndTime: 2019-03-08 20:16:42
% DurationCPUTime: 0.13s
% Computational Cost: add. (251->48), mult. (540->102), div. (0->0), fcn. (326->8), ass. (0->46)
t42 = qJD(2) ^ 2;
t34 = t42 / 0.2e1;
t43 = qJD(1) ^ 2;
t55 = t43 / 0.2e1;
t38 = sin(qJ(4));
t40 = cos(qJ(4));
t39 = sin(qJ(2));
t35 = sin(pkin(6));
t53 = qJD(1) * t35;
t46 = t39 * t53;
t14 = t46 + (pkin(4) * t38 - pkin(9) * t40 + qJ(3)) * qJD(2);
t37 = sin(qJ(5));
t54 = cos(qJ(5));
t41 = cos(qJ(2));
t44 = -t41 * t53 + qJD(3);
t18 = (-pkin(2) - pkin(8)) * qJD(2) + t44;
t36 = cos(pkin(6));
t52 = qJD(1) * t36;
t10 = t38 * t18 + t40 * t52;
t8 = qJD(4) * pkin(9) + t10;
t4 = t37 * t14 + t54 * t8;
t51 = qJD(2) * t40;
t25 = qJD(2) * qJ(3) + t46;
t50 = t25 * qJD(2);
t49 = t38 * qJD(2);
t48 = t25 ^ 2 / 0.2e1;
t47 = qJD(2) * qJD(4);
t3 = t54 * t14 - t37 * t8;
t45 = qJD(2) * t53;
t9 = t40 * t18 - t38 * t52;
t7 = -qJD(4) * pkin(4) - t9;
t31 = t36 ^ 2 * t55;
t30 = qJD(5) + t49;
t27 = t30 ^ 2 / 0.2e1;
t24 = t37 * qJD(4) + t54 * t51;
t22 = -t54 * qJD(4) + t37 * t51;
t21 = -qJD(2) * pkin(2) + t44;
t20 = t24 ^ 2 / 0.2e1;
t19 = t22 ^ 2 / 0.2e1;
t16 = t24 * t30;
t15 = t22 * t30;
t13 = t24 * t22;
t5 = t22 * pkin(5) + qJD(6) + t7;
t2 = -t22 * qJ(6) + t4;
t1 = t30 * pkin(5) - t24 * qJ(6) + t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, 0, 0, 0, 0, t34, t41 * t45, -t39 * t45, 0, t31 + (t39 ^ 2 / 0.2e1 + t41 ^ 2 / 0.2e1) * t43 * t35 ^ 2, t34, 0, 0, 0, 0, 0, 0, t21 * qJD(2), t50, t31 + t48 + t21 ^ 2 / 0.2e1, t40 ^ 2 * t34, -t40 * t42 * t38, t40 * t47, t38 ^ 2 * t34, -t38 * t47, qJD(4) ^ 2 / 0.2e1, t9 * qJD(4) + t25 * t49, -t10 * qJD(4) + t40 * t50 (-t10 * t38 - t40 * t9) * qJD(2), t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t48, t20, -t13, t16, t19, -t15, t27, t7 * t22 + t3 * t30, t7 * t24 - t4 * t30, -t4 * t22 - t3 * t24, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t20, -t13, t16, t19, -t15, t27, t1 * t30 + t5 * t22, -t2 * t30 + t5 * t24, -t1 * t24 - t2 * t22, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t6;
