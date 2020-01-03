% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRR14_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR14_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:19:41
% EndTime: 2019-12-31 19:19:41
% DurationCPUTime: 0.28s
% Computational Cost: add. (991->66), mult. (3293->164), div. (0->0), fcn. (2655->12), ass. (0->56)
t49 = sin(pkin(11));
t52 = cos(pkin(11));
t51 = sin(pkin(5));
t67 = qJD(1) * t51;
t61 = qJ(2) * t67;
t54 = cos(pkin(5));
t66 = qJD(1) * t54;
t64 = pkin(1) * t66;
t40 = t49 * t64 + t52 * t61;
t50 = sin(pkin(6));
t53 = cos(pkin(6));
t69 = t51 * t52;
t27 = (t50 * t54 + t53 * t69) * qJD(1) * pkin(8) + t40;
t45 = t52 * t64;
t71 = t49 * t51;
t29 = t45 + (pkin(2) * t54 + (-pkin(8) * t53 - qJ(2)) * t71) * qJD(1);
t35 = qJD(2) + (-pkin(8) * t49 * t50 - pkin(2) * t52 - pkin(1)) * t67;
t57 = sin(qJ(3));
t58 = cos(qJ(3));
t13 = -t57 * t27 + (t29 * t53 + t35 * t50) * t58;
t59 = qJD(1) ^ 2;
t75 = t59 / 0.2e1;
t68 = t53 * t57;
t70 = t50 * t57;
t14 = t58 * t27 + t29 * t68 + t35 * t70;
t63 = t52 * t67;
t37 = t50 * t63 - t53 * t66 - qJD(3);
t12 = -t37 * pkin(9) + t14;
t56 = sin(qJ(4));
t74 = cos(qJ(4));
t18 = -t50 * t29 + t53 * t35;
t30 = t57 * t49 * t67 + (-t50 * t66 - t53 * t63) * t58;
t32 = (t54 * t70 + (t49 * t58 + t52 * t68) * t51) * qJD(1);
t9 = t30 * pkin(3) - t32 * pkin(9) + t18;
t6 = t74 * t12 + t56 * t9;
t73 = cos(qJ(5));
t72 = t51 ^ 2 * t59;
t65 = t51 * t54 * t59;
t62 = t72 / 0.2e1;
t20 = t56 * t32 + t74 * t37;
t5 = -t56 * t12 + t74 * t9;
t11 = t37 * pkin(3) - t13;
t55 = sin(qJ(5));
t46 = -pkin(1) * t67 + qJD(2);
t39 = -t49 * t61 + t45;
t28 = qJD(4) + t30;
t22 = t74 * t32 - t56 * t37;
t19 = qJD(5) + t20;
t17 = t73 * t22 + t55 * t28;
t15 = t55 * t22 - t73 * t28;
t7 = t20 * pkin(4) - t22 * pkin(10) + t11;
t4 = t28 * pkin(10) + t6;
t3 = -t28 * pkin(4) - t5;
t2 = t73 * t4 + t55 * t7;
t1 = -t55 * t4 + t73 * t7;
t8 = [0, 0, 0, 0, 0, t75, 0, 0, 0, 0, t49 ^ 2 * t62, t49 * t52 * t72, t49 * t65, t52 ^ 2 * t62, t52 * t65, t54 ^ 2 * t75, (t39 * t54 - t46 * t69) * qJD(1), (-t40 * t54 + t46 * t71) * qJD(1), (-t39 * t49 + t40 * t52) * t67, t40 ^ 2 / 0.2e1 + t39 ^ 2 / 0.2e1 + t46 ^ 2 / 0.2e1, t32 ^ 2 / 0.2e1, -t32 * t30, -t32 * t37, t30 ^ 2 / 0.2e1, t30 * t37, t37 ^ 2 / 0.2e1, -t13 * t37 + t18 * t30, t14 * t37 + t18 * t32, -t13 * t32 - t14 * t30, t14 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * t28, t20 ^ 2 / 0.2e1, -t20 * t28, t28 ^ 2 / 0.2e1, t11 * t20 + t5 * t28, t11 * t22 - t6 * t28, -t6 * t20 - t5 * t22, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1, t17 ^ 2 / 0.2e1, -t17 * t15, t17 * t19, t15 ^ 2 / 0.2e1, -t15 * t19, t19 ^ 2 / 0.2e1, t1 * t19 + t3 * t15, t3 * t17 - t2 * t19, -t1 * t17 - t2 * t15, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t8;
