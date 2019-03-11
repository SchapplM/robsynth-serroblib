% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRPR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:13:59
% EndTime: 2019-03-08 23:14:00
% DurationCPUTime: 0.17s
% Computational Cost: add. (423->59), mult. (984->134), div. (0->0), fcn. (689->10), ass. (0->54)
t40 = sin(qJ(4));
t41 = sin(qJ(3));
t43 = cos(qJ(4));
t44 = cos(qJ(3));
t26 = (t40 * t44 + t41 * t43) * qJD(2);
t46 = qJD(2) ^ 2;
t65 = t46 / 0.2e1;
t64 = pkin(4) + pkin(10);
t63 = cos(qJ(6));
t56 = qJD(2) * t44;
t57 = qJD(2) * t41;
t24 = t40 * t57 - t43 * t56;
t62 = t26 * t24;
t36 = qJD(3) + qJD(4);
t61 = t26 * t36;
t60 = t36 * t24;
t42 = sin(qJ(2));
t37 = sin(pkin(6));
t59 = qJD(1) * t37;
t29 = qJD(2) * pkin(8) + t42 * t59;
t38 = cos(pkin(6));
t58 = qJD(1) * t38;
t32 = t44 * t58;
t14 = qJD(3) * pkin(3) + t32 + (-pkin(9) * qJD(2) - t29) * t41;
t20 = t44 * t29 + t41 * t58;
t15 = pkin(9) * t56 + t20;
t9 = t40 * t14 + t43 * t15;
t55 = t24 ^ 2 / 0.2e1;
t54 = t26 ^ 2 / 0.2e1;
t53 = qJD(2) * qJD(3);
t45 = cos(qJ(2));
t52 = t45 * t59;
t51 = qJD(2) * t59;
t8 = t43 * t14 - t40 * t15;
t6 = -t36 * qJ(5) - t9;
t50 = qJD(5) - t8;
t22 = -t52 + (-pkin(3) * t44 - pkin(2)) * qJD(2);
t48 = -t26 * qJ(5) + t22;
t47 = qJD(1) ^ 2;
t39 = sin(qJ(6));
t34 = t36 ^ 2 / 0.2e1;
t30 = -qJD(2) * pkin(2) - t52;
t23 = qJD(6) + t26;
t19 = -t41 * t29 + t32;
t18 = t39 * t24 + t63 * t36;
t16 = -t63 * t24 + t39 * t36;
t10 = t24 * pkin(4) + t48;
t7 = t64 * t24 + t48;
t5 = -t36 * pkin(4) + t50;
t4 = -t24 * pkin(5) - t6;
t3 = t26 * pkin(5) - t64 * t36 + t50;
t2 = t39 * t3 + t63 * t7;
t1 = t63 * t3 - t39 * t7;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t47 / 0.2e1, 0, 0, 0, 0, 0, t65, t45 * t51, -t42 * t51, 0 (t38 ^ 2 / 0.2e1 + (t42 ^ 2 / 0.2e1 + t45 ^ 2 / 0.2e1) * t37 ^ 2) * t47, t41 ^ 2 * t65, t41 * t46 * t44, t41 * t53, t44 ^ 2 * t65, t44 * t53, qJD(3) ^ 2 / 0.2e1, t19 * qJD(3) - t30 * t56, -t20 * qJD(3) + t30 * t57 (-t19 * t41 + t20 * t44) * qJD(2), t20 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1, t54, -t62, t61, t55, -t60, t34, t22 * t24 + t8 * t36, t22 * t26 - t9 * t36, -t9 * t24 - t8 * t26, t9 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t34, -t61, t60, t54, -t62, t55, t6 * t24 + t5 * t26, -t10 * t24 + t5 * t36, -t10 * t26 - t6 * t36, t10 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t23, t16 ^ 2 / 0.2e1, -t16 * t23, t23 ^ 2 / 0.2e1, t1 * t23 + t4 * t16, t4 * t18 - t2 * t23, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1;];
T_reg  = t11;
