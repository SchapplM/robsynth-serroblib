% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRRP6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP6_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:34:24
% EndTime: 2019-03-09 00:34:24
% DurationCPUTime: 0.23s
% Computational Cost: add. (724->58), mult. (1764->140), div. (0->0), fcn. (1400->12), ass. (0->56)
t55 = cos(qJ(2));
t47 = sin(pkin(6));
t67 = qJD(1) * t47;
t37 = qJD(2) * pkin(2) + t55 * t67;
t46 = sin(pkin(7));
t48 = cos(pkin(7));
t49 = cos(pkin(6));
t66 = qJD(1) * t49;
t74 = t37 * t48 + t46 * t66;
t53 = sin(qJ(2));
t65 = qJD(2) * t46;
t35 = pkin(9) * t65 + t53 * t67;
t52 = sin(qJ(3));
t54 = cos(qJ(3));
t18 = -t52 * t35 + t74 * t54;
t43 = t48 * qJD(2) + qJD(3);
t16 = -t43 * pkin(3) - t18;
t51 = sin(qJ(4));
t62 = t52 * t65;
t73 = cos(qJ(4));
t29 = -t73 * t43 + t51 * t62;
t31 = t51 * t43 + t73 * t62;
t12 = t29 * pkin(4) - t31 * pkin(11) + t16;
t50 = sin(qJ(5));
t72 = cos(qJ(5));
t19 = t54 * t35 + t74 * t52;
t17 = t43 * pkin(10) + t19;
t42 = t48 * t66;
t21 = t42 + (-t37 + (-pkin(3) * t54 - pkin(10) * t52) * qJD(2)) * t46;
t11 = t73 * t17 + t51 * t21;
t61 = t54 * t65;
t40 = -qJD(4) + t61;
t8 = -t40 * pkin(11) + t11;
t4 = t50 * t12 + t72 * t8;
t23 = t50 * t31 + t72 * t40;
t25 = t72 * t31 - t50 * t40;
t71 = t25 * t23;
t28 = qJD(5) + t29;
t70 = t28 * t23;
t56 = qJD(2) ^ 2;
t68 = t46 ^ 2 * t56;
t64 = t23 ^ 2 / 0.2e1;
t60 = t68 / 0.2e1;
t59 = qJD(2) * t67;
t10 = -t51 * t17 + t73 * t21;
t3 = t72 * t12 - t50 * t8;
t7 = t40 * pkin(4) - t10;
t57 = qJD(1) ^ 2;
t27 = -t46 * t37 + t42;
t26 = t28 ^ 2 / 0.2e1;
t22 = t25 ^ 2 / 0.2e1;
t13 = t25 * t28;
t5 = t23 * pkin(5) - t25 * qJ(6) + t7;
t2 = t28 * qJ(6) + t4;
t1 = -t28 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t57 / 0.2e1, 0, 0, 0, 0, 0, t56 / 0.2e1, t55 * t59, -t53 * t59, 0 (t49 ^ 2 / 0.2e1 + (t53 ^ 2 / 0.2e1 + t55 ^ 2 / 0.2e1) * t47 ^ 2) * t57, t52 ^ 2 * t60, t52 * t54 * t68, t43 * t62, t54 ^ 2 * t60, t43 * t61, t43 ^ 2 / 0.2e1, t18 * t43 - t27 * t61, -t19 * t43 + t27 * t62 (-t18 * t52 + t19 * t54) * t65, t19 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1, t31 ^ 2 / 0.2e1, -t31 * t29, -t31 * t40, t29 ^ 2 / 0.2e1, t29 * t40, t40 ^ 2 / 0.2e1, -t10 * t40 + t16 * t29, t11 * t40 + t16 * t31, -t10 * t31 - t11 * t29, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1, t22, -t71, t13, t64, -t70, t26, t7 * t23 + t3 * t28, t7 * t25 - t4 * t28, -t4 * t23 - t3 * t25, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t22, t13, t71, t26, t70, t64, -t1 * t28 + t5 * t23, t1 * t25 - t2 * t23, t2 * t28 - t5 * t25, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
