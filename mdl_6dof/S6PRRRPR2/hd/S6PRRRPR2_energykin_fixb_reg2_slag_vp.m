% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRPR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR2_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:08:57
% EndTime: 2019-03-08 23:08:57
% DurationCPUTime: 0.17s
% Computational Cost: add. (672->62), mult. (1509->151), div. (0->0), fcn. (1124->12), ass. (0->52)
t53 = qJD(2) ^ 2;
t66 = t53 / 0.2e1;
t50 = sin(qJ(2));
t45 = sin(pkin(6));
t62 = qJD(1) * t45;
t36 = qJD(2) * pkin(8) + t50 * t62;
t51 = cos(qJ(3));
t46 = cos(pkin(6));
t61 = qJD(1) * t46;
t39 = t51 * t61;
t49 = sin(qJ(3));
t22 = qJD(3) * pkin(3) + t39 + (-pkin(9) * qJD(2) - t36) * t49;
t29 = t51 * t36 + t49 * t61;
t59 = qJD(2) * t51;
t24 = pkin(9) * t59 + t29;
t48 = sin(qJ(4));
t65 = cos(qJ(4));
t12 = t48 * t22 + t65 * t24;
t43 = qJD(3) + qJD(4);
t10 = t43 * qJ(5) + t12;
t52 = cos(qJ(2));
t56 = t52 * t62;
t30 = -t56 + (-pkin(3) * t51 - pkin(2)) * qJD(2);
t60 = qJD(2) * t49;
t32 = t48 * t60 - t65 * t59;
t34 = (t48 * t51 + t65 * t49) * qJD(2);
t18 = t32 * pkin(4) - t34 * qJ(5) + t30;
t44 = sin(pkin(12));
t63 = cos(pkin(12));
t6 = t63 * t10 + t44 * t18;
t64 = cos(qJ(6));
t58 = t32 ^ 2 / 0.2e1;
t57 = qJD(2) * qJD(3);
t55 = qJD(2) * t62;
t5 = -t44 * t10 + t63 * t18;
t11 = t65 * t22 - t48 * t24;
t9 = -t43 * pkin(4) + qJD(5) - t11;
t54 = qJD(1) ^ 2;
t47 = sin(qJ(6));
t37 = -qJD(2) * pkin(2) - t56;
t31 = qJD(6) + t32;
t28 = -t49 * t36 + t39;
t27 = t63 * t34 + t44 * t43;
t25 = t44 * t34 - t63 * t43;
t15 = -t47 * t25 + t64 * t27;
t13 = t64 * t25 + t47 * t27;
t7 = t25 * pkin(5) + t9;
t4 = -t25 * pkin(10) + t6;
t3 = t32 * pkin(5) - t27 * pkin(10) + t5;
t2 = t47 * t3 + t64 * t4;
t1 = t64 * t3 - t47 * t4;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t54 / 0.2e1, 0, 0, 0, 0, 0, t66, t52 * t55, -t50 * t55, 0 (t46 ^ 2 / 0.2e1 + (t50 ^ 2 / 0.2e1 + t52 ^ 2 / 0.2e1) * t45 ^ 2) * t54, t49 ^ 2 * t66, t49 * t53 * t51, t49 * t57, t51 ^ 2 * t66, t51 * t57, qJD(3) ^ 2 / 0.2e1, t28 * qJD(3) - t37 * t59, -t29 * qJD(3) + t37 * t60 (-t28 * t49 + t29 * t51) * qJD(2), t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t37 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, t34 * t43, t58, -t32 * t43, t43 ^ 2 / 0.2e1, t11 * t43 + t30 * t32, -t12 * t43 + t30 * t34, -t11 * t34 - t12 * t32, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t32, t25 ^ 2 / 0.2e1, -t25 * t32, t58, t9 * t25 + t5 * t32, t9 * t27 - t6 * t32, -t6 * t25 - t5 * t27, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t15 ^ 2 / 0.2e1, -t15 * t13, t15 * t31, t13 ^ 2 / 0.2e1, -t13 * t31, t31 ^ 2 / 0.2e1, t1 * t31 + t7 * t13, t7 * t15 - t2 * t31, -t1 * t15 - t2 * t13, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
