% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPPR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR5_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:21:04
% EndTime: 2019-03-08 21:21:04
% DurationCPUTime: 0.15s
% Computational Cost: add. (385->59), mult. (879->131), div. (0->0), fcn. (560->10), ass. (0->52)
t48 = qJD(2) ^ 2;
t63 = t48 / 0.2e1;
t45 = sin(qJ(2));
t40 = sin(pkin(6));
t60 = qJD(1) * t40;
t28 = qJD(2) * pkin(8) + t45 * t60;
t44 = sin(qJ(3));
t46 = cos(qJ(3));
t42 = cos(pkin(6));
t59 = qJD(1) * t42;
t19 = -t44 * t28 + t46 * t59;
t50 = qJD(4) - t19;
t57 = t44 * qJD(2);
t61 = -pkin(3) - qJ(5);
t10 = pkin(4) * t57 + t61 * qJD(3) + t50;
t51 = -qJ(4) * t44 - pkin(2);
t47 = cos(qJ(2));
t55 = t47 * t60;
t16 = -t55 + (t61 * t46 + t51) * qJD(2);
t39 = sin(pkin(11));
t41 = cos(pkin(11));
t6 = t39 * t10 + t41 * t16;
t62 = cos(qJ(6));
t20 = t46 * t28 + t44 * t59;
t58 = qJD(2) * t46;
t56 = qJD(2) * qJD(3);
t18 = -qJD(3) * qJ(4) - t20;
t5 = t41 * t10 - t39 * t16;
t54 = qJD(2) * t60;
t53 = t44 * t56;
t52 = t46 * t56;
t14 = pkin(4) * t58 + qJD(5) - t18;
t49 = qJD(1) ^ 2;
t43 = sin(qJ(6));
t37 = qJD(3) ^ 2 / 0.2e1;
t34 = t46 ^ 2 * t63;
t33 = t44 ^ 2 * t63;
t32 = qJD(6) + t57;
t31 = t44 * t48 * t46;
t29 = -qJD(2) * pkin(2) - t55;
t27 = t41 * qJD(3) - t39 * t58;
t25 = t39 * qJD(3) + t41 * t58;
t21 = -t55 + (-pkin(3) * t46 + t51) * qJD(2);
t17 = -qJD(3) * pkin(3) + t50;
t13 = -t43 * t25 + t62 * t27;
t11 = t62 * t25 + t43 * t27;
t7 = t25 * pkin(5) + t14;
t4 = -t25 * pkin(9) + t6;
t3 = pkin(5) * t57 - t27 * pkin(9) + t5;
t2 = t43 * t3 + t62 * t4;
t1 = t62 * t3 - t43 * t4;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t49 / 0.2e1, 0, 0, 0, 0, 0, t63, t47 * t54, -t45 * t54, 0 (t42 ^ 2 / 0.2e1 + (t45 ^ 2 / 0.2e1 + t47 ^ 2 / 0.2e1) * t40 ^ 2) * t49, t33, t31, t53, t34, t52, t37, t19 * qJD(3) - t29 * t58, -t20 * qJD(3) + t29 * t57 (-t19 * t44 + t20 * t46) * qJD(2), t20 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1, t37, -t53, -t52, t33, t31, t34 (t17 * t44 - t18 * t46) * qJD(2), t17 * qJD(3) + t21 * t58, -t18 * qJD(3) - t21 * t57, t21 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t57, t25 ^ 2 / 0.2e1, -t25 * t57, t33, t14 * t25 + t5 * t57, t14 * t27 - t6 * t57, -t6 * t25 - t5 * t27, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t13 * t11, t13 * t32, t11 ^ 2 / 0.2e1, -t11 * t32, t32 ^ 2 / 0.2e1, t1 * t32 + t7 * t11, t7 * t13 - t2 * t32, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
