% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPRP2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP2_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:31:56
% EndTime: 2019-03-09 08:31:56
% DurationCPUTime: 0.18s
% Computational Cost: add. (483->61), mult. (1201->121), div. (0->0), fcn. (782->6), ass. (0->54)
t43 = sin(pkin(9));
t44 = cos(pkin(9));
t46 = sin(qJ(2));
t47 = cos(qJ(2));
t31 = (t43 * t47 + t44 * t46) * qJD(1);
t48 = qJD(1) ^ 2;
t67 = t48 / 0.2e1;
t66 = pkin(3) + pkin(8);
t61 = qJD(1) * t46;
t62 = pkin(7) + qJ(3);
t35 = qJD(2) * pkin(2) - t62 * t61;
t60 = qJD(1) * t47;
t36 = t62 * t60;
t18 = t44 * t35 - t43 * t36;
t52 = qJD(4) - t18;
t10 = t31 * pkin(4) - t66 * qJD(2) + t52;
t45 = sin(qJ(5));
t65 = cos(qJ(5));
t29 = t43 * t61 - t44 * t60;
t37 = qJD(3) + (-pkin(2) * t47 - pkin(1)) * qJD(1);
t50 = -t31 * qJ(4) + t37;
t7 = t66 * t29 + t50;
t4 = t45 * t10 + t65 * t7;
t64 = t31 * t29;
t63 = t47 * t48;
t19 = t43 * t35 + t44 * t36;
t59 = qJD(2) * t29;
t58 = t31 * qJD(2);
t57 = t29 ^ 2 / 0.2e1;
t56 = t31 ^ 2 / 0.2e1;
t55 = qJD(1) * qJD(2);
t17 = -qJD(2) * qJ(4) - t19;
t3 = t65 * t10 - t45 * t7;
t54 = t46 * t55;
t53 = t47 * t55;
t11 = -t29 * pkin(4) - t17;
t42 = t47 ^ 2;
t41 = t46 ^ 2;
t39 = qJD(2) ^ 2 / 0.2e1;
t28 = qJD(5) + t31;
t25 = t28 ^ 2 / 0.2e1;
t24 = t65 * qJD(2) + t45 * t29;
t22 = t45 * qJD(2) - t65 * t29;
t21 = t24 ^ 2 / 0.2e1;
t20 = t22 ^ 2 / 0.2e1;
t16 = -qJD(2) * pkin(3) + t52;
t15 = t24 * t28;
t14 = t22 * t28;
t13 = t29 * pkin(3) + t50;
t12 = t24 * t22;
t5 = t22 * pkin(5) + qJD(6) + t11;
t2 = -t22 * qJ(6) + t4;
t1 = t28 * pkin(5) - t24 * qJ(6) + t3;
t6 = [0, 0, 0, 0, 0, t67, 0, 0, 0, 0, t41 * t67, t46 * t63, t54, t42 * t67, t53, t39, pkin(1) * t63 - pkin(7) * t54, -t48 * pkin(1) * t46 - pkin(7) * t53 (t41 + t42) * t48 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t42 / 0.2e1 + t41 / 0.2e1) * pkin(7) ^ 2) * t48, t56, -t64, t58, t57, -t59, t39, t18 * qJD(2) + t37 * t29, -t19 * qJD(2) + t37 * t31, -t18 * t31 - t19 * t29, t19 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1 + t37 ^ 2 / 0.2e1, t39, -t58, t59, t56, -t64, t57, t16 * t31 + t17 * t29, t16 * qJD(2) - t13 * t29, -t17 * qJD(2) - t13 * t31, t13 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1, t21, -t12, t15, t20, -t14, t25, t11 * t22 + t3 * t28, t11 * t24 - t4 * t28, -t4 * t22 - t3 * t24, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1, t21, -t12, t15, t20, -t14, t25, t1 * t28 + t5 * t22, -t2 * t28 + t5 * t24, -t1 * t24 - t2 * t22, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t6;
