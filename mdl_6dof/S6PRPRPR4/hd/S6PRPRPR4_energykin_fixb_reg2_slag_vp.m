% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRPR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR4_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:41:16
% EndTime: 2019-03-08 19:41:16
% DurationCPUTime: 0.16s
% Computational Cost: add. (587->60), mult. (1459->143), div. (0->0), fcn. (1112->12), ass. (0->50)
t52 = qJD(2) ^ 2;
t64 = t52 / 0.2e1;
t50 = sin(qJ(2));
t45 = sin(pkin(6));
t60 = qJD(1) * t45;
t37 = qJD(2) * qJ(3) + t50 * t60;
t46 = cos(pkin(11));
t47 = cos(pkin(6));
t59 = qJD(1) * t47;
t39 = t46 * t59;
t44 = sin(pkin(11));
t22 = t39 + (-pkin(8) * qJD(2) - t37) * t44;
t29 = t46 * t37 + t44 * t59;
t57 = qJD(2) * t46;
t23 = pkin(8) * t57 + t29;
t49 = sin(qJ(4));
t63 = cos(qJ(4));
t12 = t49 * t22 + t63 * t23;
t10 = qJD(4) * qJ(5) + t12;
t51 = cos(qJ(2));
t54 = -t51 * t60 + qJD(3);
t30 = (-pkin(3) * t46 - pkin(2)) * qJD(2) + t54;
t58 = qJD(2) * t44;
t32 = t49 * t58 - t63 * t57;
t34 = (t63 * t44 + t46 * t49) * qJD(2);
t18 = t32 * pkin(4) - t34 * qJ(5) + t30;
t43 = sin(pkin(12));
t61 = cos(pkin(12));
t6 = t61 * t10 + t43 * t18;
t62 = cos(qJ(6));
t56 = t32 ^ 2 / 0.2e1;
t55 = qJD(2) * t60;
t5 = -t43 * t10 + t61 * t18;
t11 = t63 * t22 - t49 * t23;
t9 = -qJD(4) * pkin(4) + qJD(5) - t11;
t53 = qJD(1) ^ 2;
t48 = sin(qJ(6));
t36 = -qJD(2) * pkin(2) + t54;
t31 = qJD(6) + t32;
t28 = -t44 * t37 + t39;
t27 = t43 * qJD(4) + t61 * t34;
t25 = -t61 * qJD(4) + t43 * t34;
t15 = -t48 * t25 + t62 * t27;
t13 = t62 * t25 + t48 * t27;
t7 = t25 * pkin(5) + t9;
t4 = -t25 * pkin(9) + t6;
t3 = t32 * pkin(5) - t27 * pkin(9) + t5;
t2 = t48 * t3 + t62 * t4;
t1 = t62 * t3 - t48 * t4;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t53 / 0.2e1, 0, 0, 0, 0, 0, t64, t51 * t55, -t50 * t55, 0 (t47 ^ 2 / 0.2e1 + (t50 ^ 2 / 0.2e1 + t51 ^ 2 / 0.2e1) * t45 ^ 2) * t53, t44 ^ 2 * t64, t44 * t52 * t46, 0, t46 ^ 2 * t64, 0, 0, -t36 * t57, t36 * t58 (-t28 * t44 + t29 * t46) * qJD(2), t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t36 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, t34 * qJD(4), t56, -t32 * qJD(4), qJD(4) ^ 2 / 0.2e1, t11 * qJD(4) + t30 * t32, -t12 * qJD(4) + t30 * t34, -t11 * t34 - t12 * t32, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t32, t25 ^ 2 / 0.2e1, -t25 * t32, t56, t9 * t25 + t5 * t32, t9 * t27 - t6 * t32, -t6 * t25 - t5 * t27, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t15 ^ 2 / 0.2e1, -t15 * t13, t15 * t31, t13 ^ 2 / 0.2e1, -t13 * t31, t31 ^ 2 / 0.2e1, t1 * t31 + t7 * t13, t7 * t15 - t2 * t31, -t1 * t15 - t2 * t13, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
