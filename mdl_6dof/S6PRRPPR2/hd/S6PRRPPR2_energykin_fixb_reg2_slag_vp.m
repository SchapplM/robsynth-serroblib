% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPPR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:06:30
% EndTime: 2019-03-08 21:06:30
% DurationCPUTime: 0.19s
% Computational Cost: add. (394->59), mult. (984->132), div. (0->0), fcn. (689->10), ass. (0->54)
t36 = sin(pkin(11));
t38 = cos(pkin(11));
t41 = sin(qJ(3));
t43 = cos(qJ(3));
t26 = (t36 * t43 + t38 * t41) * qJD(2);
t45 = qJD(2) ^ 2;
t65 = t45 / 0.2e1;
t64 = pkin(4) + pkin(9);
t63 = cos(qJ(6));
t58 = qJD(2) * t43;
t59 = qJD(2) * t41;
t24 = t36 * t59 - t38 * t58;
t62 = t26 * t24;
t42 = sin(qJ(2));
t37 = sin(pkin(6));
t61 = qJD(1) * t37;
t29 = qJD(2) * pkin(8) + t42 * t61;
t39 = cos(pkin(6));
t60 = qJD(1) * t39;
t33 = t43 * t60;
t53 = qJ(4) * qJD(2);
t14 = qJD(3) * pkin(3) + t33 + (-t29 - t53) * t41;
t20 = t43 * t29 + t41 * t60;
t15 = t43 * t53 + t20;
t9 = t36 * t14 + t38 * t15;
t57 = qJD(3) * t24;
t56 = t26 * qJD(3);
t55 = t24 ^ 2 / 0.2e1;
t54 = t26 ^ 2 / 0.2e1;
t52 = qJD(2) * qJD(3);
t44 = cos(qJ(2));
t51 = t44 * t61;
t50 = qJD(2) * t61;
t8 = t38 * t14 - t36 * t15;
t49 = qJD(5) - t8;
t6 = -qJD(3) * qJ(5) - t9;
t22 = -t51 + qJD(4) + (-pkin(3) * t43 - pkin(2)) * qJD(2);
t47 = -t26 * qJ(5) + t22;
t46 = qJD(1) ^ 2;
t40 = sin(qJ(6));
t35 = qJD(3) ^ 2 / 0.2e1;
t30 = -qJD(2) * pkin(2) - t51;
t23 = qJD(6) + t26;
t19 = -t41 * t29 + t33;
t18 = t63 * qJD(3) + t40 * t24;
t16 = t40 * qJD(3) - t63 * t24;
t10 = t24 * pkin(4) + t47;
t7 = t64 * t24 + t47;
t5 = -qJD(3) * pkin(4) + t49;
t4 = -t24 * pkin(5) - t6;
t3 = t26 * pkin(5) - t64 * qJD(3) + t49;
t2 = t40 * t3 + t63 * t7;
t1 = t63 * t3 - t40 * t7;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t46 / 0.2e1, 0, 0, 0, 0, 0, t65, t44 * t50, -t42 * t50, 0 (t39 ^ 2 / 0.2e1 + (t42 ^ 2 / 0.2e1 + t44 ^ 2 / 0.2e1) * t37 ^ 2) * t46, t41 ^ 2 * t65, t41 * t45 * t43, t41 * t52, t43 ^ 2 * t65, t43 * t52, t35, t19 * qJD(3) - t30 * t58, -t20 * qJD(3) + t30 * t59 (-t19 * t41 + t20 * t43) * qJD(2), t20 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1, t54, -t62, t56, t55, -t57, t35, t8 * qJD(3) + t22 * t24, -t9 * qJD(3) + t22 * t26, -t9 * t24 - t8 * t26, t9 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t35, -t56, t57, t54, -t62, t55, t6 * t24 + t5 * t26, t5 * qJD(3) - t10 * t24, -t6 * qJD(3) - t10 * t26, t10 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t23, t16 ^ 2 / 0.2e1, -t16 * t23, t23 ^ 2 / 0.2e1, t1 * t23 + t4 * t16, t4 * t18 - t2 * t23, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1;];
T_reg  = t11;
