% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR4_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:10:15
% EndTime: 2019-03-09 05:10:15
% DurationCPUTime: 0.22s
% Computational Cost: add. (1071->67), mult. (2702->152), div. (0->0), fcn. (2082->10), ass. (0->50)
t57 = qJD(1) ^ 2;
t67 = t57 / 0.2e1;
t52 = sin(pkin(10));
t61 = qJD(1) * t52;
t63 = pkin(7) + qJ(2);
t42 = t63 * t61;
t53 = cos(pkin(10));
t60 = qJD(1) * t53;
t43 = t63 * t60;
t56 = sin(qJ(3));
t66 = cos(qJ(3));
t32 = -t66 * t42 - t56 * t43;
t41 = (t66 * t52 + t53 * t56) * qJD(1);
t23 = qJD(3) * pkin(3) - t41 * pkin(8) + t32;
t33 = -t56 * t42 + t66 * t43;
t39 = t56 * t61 - t66 * t60;
t27 = -t39 * pkin(8) + t33;
t55 = sin(qJ(4));
t65 = cos(qJ(4));
t14 = t55 * t23 + t65 * t27;
t50 = qJD(3) + qJD(4);
t10 = t50 * qJ(5) + t14;
t29 = t65 * t39 + t55 * t41;
t31 = -t55 * t39 + t65 * t41;
t44 = qJD(2) + (-pkin(2) * t53 - pkin(1)) * qJD(1);
t34 = t39 * pkin(3) + t44;
t15 = t29 * pkin(4) - t31 * qJ(5) + t34;
t51 = sin(pkin(11));
t62 = cos(pkin(11));
t6 = t62 * t10 + t51 * t15;
t64 = cos(qJ(6));
t59 = t29 ^ 2 / 0.2e1;
t5 = -t51 * t10 + t62 * t15;
t13 = t65 * t23 - t55 * t27;
t9 = -t50 * pkin(4) + qJD(5) - t13;
t54 = sin(qJ(6));
t49 = t53 ^ 2;
t48 = t52 ^ 2;
t47 = -qJD(1) * pkin(1) + qJD(2);
t28 = qJD(6) + t29;
t26 = t62 * t31 + t51 * t50;
t24 = t51 * t31 - t62 * t50;
t18 = -t54 * t24 + t64 * t26;
t16 = t64 * t24 + t54 * t26;
t7 = t24 * pkin(5) + t9;
t4 = -t24 * pkin(9) + t6;
t3 = t29 * pkin(5) - t26 * pkin(9) + t5;
t2 = t54 * t3 + t64 * t4;
t1 = t64 * t3 - t54 * t4;
t8 = [0, 0, 0, 0, 0, t67, 0, 0, 0, 0, t48 * t67, t52 * t57 * t53, 0, t49 * t67, 0, 0, -t47 * t60, t47 * t61 (t48 + t49) * t57 * qJ(2), t47 ^ 2 / 0.2e1 + (t49 / 0.2e1 + t48 / 0.2e1) * qJ(2) ^ 2 * t57, t41 ^ 2 / 0.2e1, -t41 * t39, t41 * qJD(3), t39 ^ 2 / 0.2e1, -t39 * qJD(3), qJD(3) ^ 2 / 0.2e1, t32 * qJD(3) + t44 * t39, -t33 * qJD(3) + t44 * t41, -t32 * t41 - t33 * t39, t33 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1 + t44 ^ 2 / 0.2e1, t31 ^ 2 / 0.2e1, -t31 * t29, t31 * t50, t59, -t29 * t50, t50 ^ 2 / 0.2e1, t13 * t50 + t34 * t29, -t14 * t50 + t34 * t31, -t13 * t31 - t14 * t29, t14 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1, t26 ^ 2 / 0.2e1, -t26 * t24, t26 * t29, t24 ^ 2 / 0.2e1, -t24 * t29, t59, t9 * t24 + t5 * t29, t9 * t26 - t6 * t29, -t6 * t24 - t5 * t26, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t28, t16 ^ 2 / 0.2e1, -t16 * t28, t28 ^ 2 / 0.2e1, t1 * t28 + t7 * t16, t7 * t18 - t2 * t28, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
