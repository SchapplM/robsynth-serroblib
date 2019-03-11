% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPRR7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR7_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:34:10
% EndTime: 2019-03-08 22:34:10
% DurationCPUTime: 0.15s
% Computational Cost: add. (395->59), mult. (878->133), div. (0->0), fcn. (560->10), ass. (0->53)
t49 = qJD(2) ^ 2;
t64 = t49 / 0.2e1;
t63 = -pkin(3) - pkin(9);
t45 = sin(qJ(2));
t40 = sin(pkin(6));
t61 = qJD(1) * t40;
t28 = qJD(2) * pkin(8) + t45 * t61;
t44 = sin(qJ(3));
t47 = cos(qJ(3));
t41 = cos(pkin(6));
t60 = qJD(1) * t41;
t19 = -t44 * t28 + t47 * t60;
t51 = qJD(4) - t19;
t58 = t44 * qJD(2);
t10 = pkin(4) * t58 + t63 * qJD(3) + t51;
t52 = -qJ(4) * t44 - pkin(2);
t48 = cos(qJ(2));
t56 = t48 * t61;
t16 = -t56 + (t63 * t47 + t52) * qJD(2);
t43 = sin(qJ(5));
t46 = cos(qJ(5));
t6 = t43 * t10 + t46 * t16;
t62 = cos(qJ(6));
t20 = t47 * t28 + t44 * t60;
t59 = qJD(2) * t47;
t57 = qJD(2) * qJD(3);
t18 = -qJD(3) * qJ(4) - t20;
t5 = t46 * t10 - t43 * t16;
t55 = qJD(2) * t61;
t54 = t44 * t57;
t53 = t47 * t57;
t14 = pkin(4) * t59 - t18;
t33 = qJD(5) + t58;
t50 = qJD(1) ^ 2;
t42 = sin(qJ(6));
t38 = qJD(3) ^ 2 / 0.2e1;
t35 = t47 ^ 2 * t64;
t34 = t44 ^ 2 * t64;
t32 = t44 * t49 * t47;
t30 = qJD(6) + t33;
t29 = -qJD(2) * pkin(2) - t56;
t27 = t46 * qJD(3) - t43 * t59;
t25 = t43 * qJD(3) + t46 * t59;
t21 = -t56 + (-pkin(3) * t47 + t52) * qJD(2);
t17 = -qJD(3) * pkin(3) + t51;
t13 = -t42 * t25 + t62 * t27;
t11 = t62 * t25 + t42 * t27;
t7 = t25 * pkin(5) + t14;
t4 = -t25 * pkin(10) + t6;
t3 = t33 * pkin(5) - t27 * pkin(10) + t5;
t2 = t42 * t3 + t62 * t4;
t1 = t62 * t3 - t42 * t4;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t50 / 0.2e1, 0, 0, 0, 0, 0, t64, t48 * t55, -t45 * t55, 0 (t41 ^ 2 / 0.2e1 + (t45 ^ 2 / 0.2e1 + t48 ^ 2 / 0.2e1) * t40 ^ 2) * t50, t34, t32, t54, t35, t53, t38, t19 * qJD(3) - t29 * t59, -t20 * qJD(3) + t29 * t58 (-t19 * t44 + t20 * t47) * qJD(2), t20 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1, t38, -t54, -t53, t34, t32, t35 (t17 * t44 - t18 * t47) * qJD(2), t17 * qJD(3) + t21 * t59, -t18 * qJD(3) - t21 * t58, t21 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t33, t25 ^ 2 / 0.2e1, -t25 * t33, t33 ^ 2 / 0.2e1, t14 * t25 + t5 * t33, t14 * t27 - t6 * t33, -t6 * t25 - t5 * t27, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t13 * t11, t13 * t30, t11 ^ 2 / 0.2e1, -t11 * t30, t30 ^ 2 / 0.2e1, t1 * t30 + t7 * t11, t7 * t13 - t2 * t30, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
