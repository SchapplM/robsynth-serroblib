% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PPRRPR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:50:53
% EndTime: 2019-03-08 18:50:53
% DurationCPUTime: 0.16s
% Computational Cost: add. (316->50), mult. (803->113), div. (0->0), fcn. (627->12), ass. (0->50)
t27 = cos(pkin(6)) * qJD(1) + qJD(2);
t34 = sin(pkin(7));
t37 = cos(pkin(7));
t36 = cos(pkin(12));
t35 = sin(pkin(6));
t58 = qJD(1) * t35;
t53 = t36 * t58;
t62 = t27 * t34 + t37 * t53;
t40 = sin(qJ(3));
t43 = cos(qJ(3));
t33 = sin(pkin(12));
t54 = t33 * t58;
t15 = -t40 * t54 + t62 * t43;
t44 = qJD(3) ^ 2;
t61 = t44 / 0.2e1;
t60 = -pkin(4) - pkin(10);
t16 = t62 * t40 + t43 * t54;
t14 = qJD(3) * pkin(9) + t16;
t18 = t37 * t27 - t34 * t53;
t39 = sin(qJ(4));
t42 = cos(qJ(4));
t8 = t42 * t14 + t39 * t18;
t57 = qJD(3) * t42;
t56 = t39 * qJD(3);
t55 = qJD(3) * qJD(4);
t52 = t39 * t55;
t51 = t42 * t55;
t50 = -qJ(5) * t39 - pkin(3);
t7 = -t39 * t14 + t42 * t18;
t48 = qJD(5) - t7;
t6 = -qJD(4) * qJ(5) - t8;
t45 = qJD(1) ^ 2;
t41 = cos(qJ(6));
t38 = sin(qJ(6));
t32 = qJD(4) ^ 2 / 0.2e1;
t30 = t42 ^ 2 * t61;
t29 = t39 ^ 2 * t61;
t28 = qJD(6) + t56;
t26 = t39 * t44 * t42;
t23 = t41 * qJD(4) - t38 * t57;
t21 = t38 * qJD(4) + t41 * t57;
t13 = -qJD(3) * pkin(3) - t15;
t10 = (-pkin(4) * t42 + t50) * qJD(3) - t15;
t9 = (t60 * t42 + t50) * qJD(3) - t15;
t5 = -qJD(4) * pkin(4) + t48;
t4 = pkin(5) * t57 - t6;
t3 = pkin(5) * t56 + t60 * qJD(4) + t48;
t2 = t38 * t3 + t41 * t9;
t1 = t41 * t3 - t38 * t9;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t45 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27 ^ 2 / 0.2e1 + (t33 ^ 2 / 0.2e1 + t36 ^ 2 / 0.2e1) * t45 * t35 ^ 2, 0, 0, 0, 0, 0, t61, t15 * qJD(3), -t16 * qJD(3), 0, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t29, t26, t52, t30, t51, t32, t7 * qJD(4) - t13 * t57, -t8 * qJD(4) + t13 * t56 (-t39 * t7 + t42 * t8) * qJD(3), t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1, t32, -t52, -t51, t29, t26, t30 (t39 * t5 - t42 * t6) * qJD(3), t5 * qJD(4) + t10 * t57, -t6 * qJD(4) - t10 * t56, t10 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1, t23 ^ 2 / 0.2e1, -t23 * t21, t23 * t28, t21 ^ 2 / 0.2e1, -t21 * t28, t28 ^ 2 / 0.2e1, t1 * t28 + t4 * t21, -t2 * t28 + t4 * t23, -t1 * t23 - t2 * t21, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1;];
T_reg  = t11;
