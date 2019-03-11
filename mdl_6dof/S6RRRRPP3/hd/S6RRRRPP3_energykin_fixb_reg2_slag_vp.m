% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPP3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:56:28
% EndTime: 2019-03-09 20:56:28
% DurationCPUTime: 0.17s
% Computational Cost: add. (595->59), mult. (1310->126), div. (0->0), fcn. (879->6), ass. (0->49)
t45 = qJD(1) ^ 2;
t60 = t45 / 0.2e1;
t59 = -pkin(8) - pkin(7);
t42 = sin(qJ(3));
t44 = cos(qJ(2));
t52 = qJD(1) * t44;
t43 = sin(qJ(2));
t53 = qJD(1) * t43;
t58 = cos(qJ(3));
t28 = t42 * t53 - t58 * t52;
t30 = (t42 * t44 + t58 * t43) * qJD(1);
t35 = (-pkin(2) * t44 - pkin(1)) * qJD(1);
t12 = t28 * pkin(3) - t30 * pkin(9) + t35;
t33 = qJD(2) * pkin(2) + t59 * t53;
t34 = t59 * t52;
t18 = t42 * t33 - t58 * t34;
t38 = qJD(2) + qJD(3);
t16 = t38 * pkin(9) + t18;
t41 = sin(qJ(4));
t57 = cos(qJ(4));
t8 = t41 * t12 + t57 * t16;
t21 = t41 * t30 - t57 * t38;
t23 = t57 * t30 + t41 * t38;
t56 = t21 * t23;
t27 = qJD(4) + t28;
t11 = t23 * t27;
t10 = t27 * t21;
t55 = t44 * t45;
t54 = pkin(4) + qJ(6);
t19 = t21 ^ 2 / 0.2e1;
t20 = t23 ^ 2 / 0.2e1;
t51 = qJD(1) * qJD(2);
t50 = t43 * t51;
t49 = t44 * t51;
t5 = -t27 * qJ(5) - t8;
t7 = t57 * t12 - t41 * t16;
t17 = t58 * t33 + t42 * t34;
t48 = qJD(5) - t7;
t15 = -t38 * pkin(3) - t17;
t47 = -t23 * qJ(5) + t15;
t40 = t44 ^ 2;
t39 = t43 ^ 2;
t25 = t27 ^ 2 / 0.2e1;
t6 = t21 * pkin(4) + t47;
t4 = -t27 * pkin(4) + t48;
t3 = t54 * t21 + t47;
t2 = -t21 * pkin(5) + qJD(6) - t5;
t1 = t23 * pkin(5) - t54 * t27 + t48;
t9 = [0, 0, 0, 0, 0, t60, 0, 0, 0, 0, t39 * t60, t43 * t55, t50, t40 * t60, t49, qJD(2) ^ 2 / 0.2e1, pkin(1) * t55 - pkin(7) * t50, -t45 * pkin(1) * t43 - pkin(7) * t49 (t39 + t40) * t45 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t40 / 0.2e1 + t39 / 0.2e1) * pkin(7) ^ 2) * t45, t30 ^ 2 / 0.2e1, -t30 * t28, t30 * t38, t28 ^ 2 / 0.2e1, -t28 * t38, t38 ^ 2 / 0.2e1, t17 * t38 + t35 * t28, -t18 * t38 + t35 * t30, -t17 * t30 - t18 * t28, t18 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t35 ^ 2 / 0.2e1, t20, -t56, t11, t19, -t10, t25, t15 * t21 + t7 * t27, t15 * t23 - t8 * t27, -t8 * t21 - t7 * t23, t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t25, -t11, t10, t20, -t56, t19, t5 * t21 + t4 * t23, -t6 * t21 + t4 * t27, -t6 * t23 - t5 * t27, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1, t25, t10, t11, t19, t56, t20, t1 * t23 - t2 * t21, t2 * t27 - t3 * t23, -t1 * t27 + t3 * t21, t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1;];
T_reg  = t9;
