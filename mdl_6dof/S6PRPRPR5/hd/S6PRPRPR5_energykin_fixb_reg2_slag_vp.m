% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRPR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR5_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:45:02
% EndTime: 2019-03-08 19:45:02
% DurationCPUTime: 0.16s
% Computational Cost: add. (371->57), mult. (943->126), div. (0->0), fcn. (677->10), ass. (0->52)
t36 = sin(pkin(11));
t38 = cos(pkin(11));
t41 = sin(qJ(4));
t43 = cos(qJ(4));
t26 = (t36 * t43 + t38 * t41) * qJD(2);
t45 = qJD(2) ^ 2;
t63 = t45 / 0.2e1;
t62 = pkin(4) + pkin(9);
t61 = cos(qJ(6));
t56 = qJD(2) * t38;
t57 = qJD(2) * t36;
t24 = t41 * t57 - t43 * t56;
t60 = t26 * t24;
t42 = sin(qJ(2));
t37 = sin(pkin(6));
t59 = qJD(1) * t37;
t30 = qJD(2) * qJ(3) + t42 * t59;
t39 = cos(pkin(6));
t58 = qJD(1) * t39;
t32 = t38 * t58;
t14 = t32 + (-pkin(8) * qJD(2) - t30) * t36;
t20 = t38 * t30 + t36 * t58;
t15 = pkin(8) * t56 + t20;
t9 = t41 * t14 + t43 * t15;
t55 = qJD(4) * t24;
t54 = t26 * qJD(4);
t53 = t24 ^ 2 / 0.2e1;
t52 = t26 ^ 2 / 0.2e1;
t51 = qJD(2) * t59;
t8 = t43 * t14 - t41 * t15;
t50 = qJD(5) - t8;
t7 = -qJD(4) * qJ(5) - t9;
t44 = cos(qJ(2));
t49 = -t44 * t59 + qJD(3);
t22 = (-pkin(3) * t38 - pkin(2)) * qJD(2) + t49;
t47 = -t26 * qJ(5) + t22;
t46 = qJD(1) ^ 2;
t40 = sin(qJ(6));
t35 = qJD(4) ^ 2 / 0.2e1;
t29 = -qJD(2) * pkin(2) + t49;
t23 = qJD(6) + t26;
t19 = -t36 * t30 + t32;
t18 = t61 * qJD(4) + t40 * t24;
t16 = t40 * qJD(4) - t61 * t24;
t10 = t24 * pkin(4) + t47;
t6 = -qJD(4) * pkin(4) + t50;
t5 = t62 * t24 + t47;
t4 = -t24 * pkin(5) - t7;
t3 = t26 * pkin(5) - t62 * qJD(4) + t50;
t2 = t40 * t3 + t61 * t5;
t1 = t61 * t3 - t40 * t5;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t46 / 0.2e1, 0, 0, 0, 0, 0, t63, t44 * t51, -t42 * t51, 0 (t39 ^ 2 / 0.2e1 + (t42 ^ 2 / 0.2e1 + t44 ^ 2 / 0.2e1) * t37 ^ 2) * t46, t36 ^ 2 * t63, t36 * t45 * t38, 0, t38 ^ 2 * t63, 0, 0, -t29 * t56, t29 * t57 (-t19 * t36 + t20 * t38) * qJD(2), t20 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1, t52, -t60, t54, t53, -t55, t35, t8 * qJD(4) + t22 * t24, -t9 * qJD(4) + t22 * t26, -t9 * t24 - t8 * t26, t9 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t35, -t54, t55, t52, -t60, t53, t7 * t24 + t6 * t26, t6 * qJD(4) - t10 * t24, -t7 * qJD(4) - t10 * t26, t10 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t23, t16 ^ 2 / 0.2e1, -t16 * t23, t23 ^ 2 / 0.2e1, t1 * t23 + t4 * t16, t4 * t18 - t2 * t23, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1;];
T_reg  = t11;
