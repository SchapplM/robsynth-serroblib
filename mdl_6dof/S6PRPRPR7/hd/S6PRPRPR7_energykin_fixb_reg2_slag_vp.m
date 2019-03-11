% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRPR7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:53:46
% EndTime: 2019-03-08 19:53:46
% DurationCPUTime: 0.15s
% Computational Cost: add. (211->52), mult. (454->104), div. (0->0), fcn. (247->8), ass. (0->49)
t37 = qJD(2) ^ 2;
t29 = t37 / 0.2e1;
t38 = qJD(1) ^ 2;
t55 = t38 / 0.2e1;
t54 = cos(qJ(6));
t36 = cos(qJ(2));
t30 = sin(pkin(6));
t53 = qJD(1) * t30;
t39 = -t36 * t53 + qJD(3);
t12 = (-pkin(2) - pkin(8)) * qJD(2) + t39;
t33 = sin(qJ(4));
t35 = cos(qJ(4));
t31 = cos(pkin(6));
t52 = qJD(1) * t31;
t9 = t33 * t12 + t35 * t52;
t51 = qJD(2) * t33;
t34 = sin(qJ(2));
t45 = t34 * t53;
t17 = qJD(2) * qJ(3) + t45;
t50 = t17 * qJD(2);
t49 = t35 * qJD(2);
t48 = t17 ^ 2 / 0.2e1;
t47 = qJD(2) * qJD(4);
t46 = t35 * t37 * t33;
t44 = qJD(2) * t53;
t43 = t33 * t47;
t42 = t35 * t47;
t19 = t33 * t52;
t8 = t35 * t12 - t19;
t41 = -qJ(5) * t35 + qJ(3);
t40 = pkin(4) * t51 + t45;
t6 = -qJD(4) * qJ(5) - t9;
t32 = sin(qJ(6));
t28 = qJD(4) ^ 2 / 0.2e1;
t25 = t35 ^ 2 * t29;
t24 = t33 ^ 2 * t29;
t23 = t31 ^ 2 * t55;
t22 = qJD(6) + t49;
t16 = t54 * qJD(4) + t32 * t51;
t14 = t32 * qJD(4) - t54 * t51;
t13 = -qJD(2) * pkin(2) + t39;
t10 = t41 * qJD(2) + t40;
t7 = (pkin(9) * t33 + t41) * qJD(2) + t40;
t5 = -qJD(4) * pkin(4) + qJD(5) - t8;
t4 = -pkin(5) * t51 - t6;
t3 = qJD(5) + t19 + (pkin(5) * qJD(2) - t12) * t35 + (-pkin(4) - pkin(9)) * qJD(4);
t2 = t32 * t3 + t54 * t7;
t1 = t54 * t3 - t32 * t7;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, 0, 0, 0, 0, t29, t36 * t44, -t34 * t44, 0, t23 + (t34 ^ 2 / 0.2e1 + t36 ^ 2 / 0.2e1) * t38 * t30 ^ 2, t29, 0, 0, 0, 0, 0, 0, t13 * qJD(2), t50, t23 + t48 + t13 ^ 2 / 0.2e1, t25, -t46, t42, t24, -t43, t28, t8 * qJD(4) + t33 * t50, -t9 * qJD(4) + t17 * t49 (-t33 * t9 - t35 * t8) * qJD(2), t9 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t48, t28, -t42, t43, t25, -t46, t24 (t33 * t6 + t35 * t5) * qJD(2), t5 * qJD(4) - t10 * t51, -t6 * qJD(4) - t10 * t49, t10 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1, t16 ^ 2 / 0.2e1, -t16 * t14, t16 * t22, t14 ^ 2 / 0.2e1, -t14 * t22, t22 ^ 2 / 0.2e1, t1 * t22 + t4 * t14, t4 * t16 - t2 * t22, -t1 * t16 - t2 * t14, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1;];
T_reg  = t11;
