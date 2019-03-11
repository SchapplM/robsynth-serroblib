% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPRR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR5_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:10:58
% EndTime: 2019-03-09 09:10:58
% DurationCPUTime: 0.20s
% Computational Cost: add. (538->64), mult. (1316->137), div. (0->0), fcn. (869->8), ass. (0->49)
t55 = cos(qJ(2));
t50 = sin(pkin(6));
t66 = qJD(1) * t50;
t61 = t55 * t66;
t54 = sin(qJ(2));
t62 = t54 * t66;
t21 = -pkin(1) * t66 - pkin(2) * t61 - qJ(3) * t62;
t18 = pkin(3) * t61 + qJD(4) - t21;
t10 = (pkin(4) * t55 - pkin(9) * t54) * t66 + t18;
t65 = cos(pkin(6)) * qJD(1);
t63 = pkin(1) * t65;
t27 = pkin(8) * t61 + t54 * t63;
t46 = qJD(2) + t65;
t20 = t46 * qJ(3) + t27;
t59 = qJ(4) * t66;
t17 = -t55 * t59 + t20;
t13 = -t46 * pkin(9) + t17;
t53 = sin(qJ(5));
t69 = cos(qJ(5));
t6 = t53 * t10 + t69 * t13;
t68 = cos(qJ(6));
t56 = qJD(1) ^ 2;
t67 = t50 ^ 2 * t56;
t26 = -pkin(8) * t62 + t55 * t63;
t32 = t46 ^ 2 / 0.2e1;
t64 = t55 * t67;
t60 = t67 / 0.2e1;
t58 = t54 * t64;
t28 = t46 * t62;
t57 = t46 * t61;
t19 = -t46 * pkin(2) + qJD(3) - t26;
t23 = t69 * t46 + t53 * t62;
t12 = -t46 * pkin(3) - t54 * t59 + t19;
t5 = t69 * t10 - t53 * t13;
t9 = t46 * pkin(4) - t12;
t52 = sin(qJ(6));
t34 = t55 ^ 2 * t60;
t33 = t54 ^ 2 * t60;
t31 = qJD(5) + t61;
t25 = -t53 * t46 + t69 * t62;
t22 = qJD(6) + t23;
t16 = t68 * t25 + t52 * t31;
t14 = t52 * t25 - t68 * t31;
t7 = t23 * pkin(5) - t25 * pkin(10) + t9;
t4 = t31 * pkin(10) + t6;
t3 = -t31 * pkin(5) - t5;
t2 = t68 * t4 + t52 * t7;
t1 = -t52 * t4 + t68 * t7;
t8 = [0, 0, 0, 0, 0, t56 / 0.2e1, 0, 0, 0, 0, t33, t58, t28, t34, t57, t32, pkin(1) * t64 + t26 * t46, -pkin(1) * t54 * t67 - t27 * t46 (-t26 * t54 + t27 * t55) * t66, t27 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t60, t33, t28, -t58, t32, -t57, t34, -t19 * t46 - t21 * t61 (t19 * t54 + t20 * t55) * t66, t20 * t46 - t21 * t62, t20 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t33, -t58, -t28, t34, t57, t32, -t12 * t46 + t18 * t61, t17 * t46 + t18 * t62 (-t12 * t54 - t17 * t55) * t66, t17 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t25 ^ 2 / 0.2e1, -t25 * t23, t25 * t31, t23 ^ 2 / 0.2e1, -t23 * t31, t31 ^ 2 / 0.2e1, t9 * t23 + t5 * t31, t9 * t25 - t6 * t31, -t6 * t23 - t5 * t25, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t16 ^ 2 / 0.2e1, -t16 * t14, t16 * t22, t14 ^ 2 / 0.2e1, -t14 * t22, t22 ^ 2 / 0.2e1, t1 * t22 + t3 * t14, t3 * t16 - t2 * t22, -t1 * t16 - t2 * t14, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
