% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRRPR4
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
% Datum: 2019-03-08 23:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRPR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR4_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:20:27
% EndTime: 2019-03-08 23:20:27
% DurationCPUTime: 0.17s
% Computational Cost: add. (711->61), mult. (1606->149), div. (0->0), fcn. (1185->12), ass. (0->51)
t56 = qJD(2) ^ 2;
t68 = t56 / 0.2e1;
t53 = sin(qJ(2));
t48 = sin(pkin(6));
t64 = qJD(1) * t48;
t37 = qJD(2) * pkin(8) + t53 * t64;
t52 = sin(qJ(3));
t54 = cos(qJ(3));
t49 = cos(pkin(6));
t63 = qJD(1) * t49;
t29 = t54 * t37 + t52 * t63;
t25 = qJD(3) * pkin(9) + t29;
t55 = cos(qJ(2));
t59 = t55 * t64;
t30 = -t59 + (-pkin(3) * t54 - pkin(9) * t52 - pkin(2)) * qJD(2);
t51 = sin(qJ(4));
t67 = cos(qJ(4));
t17 = t67 * t25 + t51 * t30;
t62 = qJD(2) * t52;
t34 = -t67 * qJD(3) + t51 * t62;
t15 = -t34 * qJ(5) + t17;
t47 = sin(pkin(12));
t65 = cos(pkin(12));
t16 = -t51 * t25 + t67 * t30;
t36 = t51 * qJD(3) + t67 * t62;
t61 = t54 * qJD(2);
t43 = -qJD(4) + t61;
t9 = -t43 * pkin(4) - t36 * qJ(5) + t16;
t6 = t65 * t15 + t47 * t9;
t66 = cos(qJ(6));
t60 = qJD(2) * qJD(3);
t5 = -t47 * t15 + t65 * t9;
t58 = qJD(2) * t64;
t28 = -t52 * t37 + t54 * t63;
t24 = -qJD(3) * pkin(3) - t28;
t18 = t34 * pkin(4) + qJD(5) + t24;
t57 = qJD(1) ^ 2;
t50 = sin(qJ(6));
t40 = -qJD(6) + t43;
t39 = t43 ^ 2 / 0.2e1;
t38 = -qJD(2) * pkin(2) - t59;
t22 = -t47 * t34 + t65 * t36;
t20 = t65 * t34 + t47 * t36;
t13 = -t50 * t20 + t66 * t22;
t11 = t66 * t20 + t50 * t22;
t10 = t20 * pkin(5) + t18;
t4 = -t20 * pkin(10) + t6;
t3 = -t43 * pkin(5) - t22 * pkin(10) + t5;
t2 = t50 * t3 + t66 * t4;
t1 = t66 * t3 - t50 * t4;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t57 / 0.2e1, 0, 0, 0, 0, 0, t68, t55 * t58, -t53 * t58, 0 (t49 ^ 2 / 0.2e1 + (t53 ^ 2 / 0.2e1 + t55 ^ 2 / 0.2e1) * t48 ^ 2) * t57, t52 ^ 2 * t68, t52 * t56 * t54, t52 * t60, t54 ^ 2 * t68, t54 * t60, qJD(3) ^ 2 / 0.2e1, t28 * qJD(3) - t38 * t61, -t29 * qJD(3) + t38 * t62 (-t28 * t52 + t29 * t54) * qJD(2), t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1, t36 ^ 2 / 0.2e1, -t36 * t34, -t36 * t43, t34 ^ 2 / 0.2e1, t34 * t43, t39, -t16 * t43 + t24 * t34, t17 * t43 + t24 * t36, -t16 * t36 - t17 * t34, t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, -t22 * t43, t20 ^ 2 / 0.2e1, t20 * t43, t39, t18 * t20 - t5 * t43, t18 * t22 + t6 * t43, -t6 * t20 - t5 * t22, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t13 * t11, -t13 * t40, t11 ^ 2 / 0.2e1, t11 * t40, t40 ^ 2 / 0.2e1, -t1 * t40 + t10 * t11, t10 * t13 + t2 * t40, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg  = t7;
