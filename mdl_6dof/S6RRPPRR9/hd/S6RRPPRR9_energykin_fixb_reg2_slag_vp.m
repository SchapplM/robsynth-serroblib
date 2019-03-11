% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPPRR9
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
% Datum: 2019-03-09 09:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPRR9_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR9_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:32:06
% EndTime: 2019-03-09 09:32:06
% DurationCPUTime: 0.18s
% Computational Cost: add. (529->68), mult. (1325->137), div. (0->0), fcn. (878->8), ass. (0->51)
t50 = sin(qJ(2));
t51 = cos(qJ(2));
t46 = sin(pkin(6));
t62 = qJD(1) * t46;
t57 = t51 * t62;
t61 = cos(pkin(6)) * qJD(1);
t59 = pkin(1) * t61;
t27 = pkin(8) * t57 + t50 * t59;
t44 = qJD(2) + t61;
t20 = -t44 * qJ(3) - t27;
t17 = pkin(3) * t57 + qJD(4) - t20;
t10 = pkin(4) * t57 - t44 * pkin(9) + t17;
t54 = -pkin(1) + (-pkin(2) - qJ(4)) * t51;
t13 = ((pkin(9) - qJ(3)) * t50 + t54) * t62;
t49 = sin(qJ(5));
t66 = cos(qJ(5));
t6 = t49 * t10 + t66 * t13;
t65 = cos(qJ(6));
t52 = qJD(1) ^ 2;
t64 = t46 ^ 2 * t52;
t58 = t50 * t62;
t26 = -pkin(8) * t58 + t51 * t59;
t63 = qJ(3) * t50;
t60 = t51 * t64;
t56 = t64 / 0.2e1;
t23 = t49 * t44 - t66 * t58;
t55 = t44 * t57;
t19 = -t44 * pkin(2) + qJD(3) - t26;
t53 = t44 * qJ(4) - t19;
t5 = t66 * t10 - t49 * t13;
t9 = (-pkin(3) - pkin(4)) * t58 + t53;
t48 = sin(qJ(6));
t35 = t51 ^ 2 * t56;
t34 = t50 ^ 2 * t56;
t33 = t44 ^ 2 / 0.2e1;
t32 = qJD(5) + t57;
t30 = t50 * t60;
t28 = t44 * t58;
t25 = t66 * t44 + t49 * t58;
t22 = qJD(6) + t23;
t21 = (-pkin(2) * t51 - pkin(1) - t63) * t62;
t18 = (t54 - t63) * t62;
t16 = t65 * t25 + t48 * t32;
t14 = t48 * t25 - t65 * t32;
t11 = pkin(3) * t58 - t53;
t7 = t23 * pkin(5) - t25 * pkin(10) + t9;
t4 = t32 * pkin(10) + t6;
t3 = -t32 * pkin(5) - t5;
t2 = t65 * t4 + t48 * t7;
t1 = -t48 * t4 + t65 * t7;
t8 = [0, 0, 0, 0, 0, t52 / 0.2e1, 0, 0, 0, 0, t34, t30, t28, t35, t55, t33, pkin(1) * t60 + t26 * t44, -pkin(1) * t50 * t64 - t27 * t44 (-t26 * t50 + t27 * t51) * t62, t27 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t56, t33, -t28, -t55, t34, t30, t35 (t19 * t50 - t20 * t51) * t62, t19 * t44 + t21 * t57, -t20 * t44 - t21 * t58, t21 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t33, -t55, t28, t35, -t30, t34 (t11 * t50 + t17 * t51) * t62, t17 * t44 - t18 * t58, -t11 * t44 - t18 * t57, t18 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t25 ^ 2 / 0.2e1, -t25 * t23, t25 * t32, t23 ^ 2 / 0.2e1, -t23 * t32, t32 ^ 2 / 0.2e1, t9 * t23 + t5 * t32, t9 * t25 - t6 * t32, -t6 * t23 - t5 * t25, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t16 ^ 2 / 0.2e1, -t16 * t14, t16 * t22, t14 ^ 2 / 0.2e1, -t14 * t22, t22 ^ 2 / 0.2e1, t1 * t22 + t3 * t14, t3 * t16 - t2 * t22, -t1 * t16 - t2 * t14, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
