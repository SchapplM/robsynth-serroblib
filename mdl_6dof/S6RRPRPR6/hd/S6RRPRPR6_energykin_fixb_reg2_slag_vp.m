% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPR6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR6_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:42:54
% EndTime: 2019-03-09 10:42:54
% DurationCPUTime: 0.24s
% Computational Cost: add. (913->70), mult. (2569->148), div. (0->0), fcn. (1988->10), ass. (0->59)
t75 = pkin(4) + pkin(10);
t74 = cos(qJ(6));
t49 = sin(pkin(11));
t51 = cos(pkin(11));
t55 = sin(qJ(2));
t57 = cos(qJ(2));
t50 = sin(pkin(6));
t69 = qJD(1) * t50;
t37 = (t49 * t57 + t51 * t55) * t69;
t68 = cos(pkin(6)) * qJD(1);
t47 = qJD(2) + t68;
t54 = sin(qJ(4));
t56 = cos(qJ(4));
t25 = t54 * t37 - t56 * t47;
t27 = t56 * t37 + t54 * t47;
t73 = t27 * t25;
t62 = t57 * t69;
t63 = t55 * t69;
t35 = t49 * t63 - t51 * t62;
t34 = qJD(4) + t35;
t72 = t27 * t34;
t71 = t34 * t25;
t58 = qJD(1) ^ 2;
t70 = t50 ^ 2 * t58;
t64 = pkin(1) * t68;
t46 = t57 * t64;
t29 = t47 * pkin(2) + t46 + (-pkin(8) - qJ(3)) * t63;
t39 = pkin(8) * t62 + t55 * t64;
t32 = qJ(3) * t62 + t39;
t16 = t49 * t29 + t51 * t32;
t14 = t47 * pkin(9) + t16;
t40 = qJD(3) + (-pkin(2) * t57 - pkin(1)) * t69;
t21 = t35 * pkin(3) - t37 * pkin(9) + t40;
t10 = t56 * t14 + t54 * t21;
t67 = t25 ^ 2 / 0.2e1;
t66 = t27 ^ 2 / 0.2e1;
t65 = t57 * t70;
t61 = t70 / 0.2e1;
t9 = -t54 * t14 + t56 * t21;
t15 = t51 * t29 - t49 * t32;
t7 = -t34 * qJ(5) - t10;
t60 = qJD(5) - t9;
t13 = -t47 * pkin(3) - t15;
t59 = -t27 * qJ(5) + t13;
t53 = sin(qJ(6));
t43 = t47 ^ 2 / 0.2e1;
t38 = -pkin(8) * t63 + t46;
t33 = t34 ^ 2 / 0.2e1;
t24 = qJD(6) + t27;
t19 = t53 * t25 + t74 * t34;
t17 = -t74 * t25 + t53 * t34;
t8 = t25 * pkin(4) + t59;
t6 = -t34 * pkin(4) + t60;
t5 = t75 * t25 + t59;
t4 = -t25 * pkin(5) - t7;
t3 = t27 * pkin(5) - t75 * t34 + t60;
t2 = t53 * t3 + t74 * t5;
t1 = t74 * t3 - t53 * t5;
t11 = [0, 0, 0, 0, 0, t58 / 0.2e1, 0, 0, 0, 0, t55 ^ 2 * t61, t55 * t65, t47 * t63, t57 ^ 2 * t61, t47 * t62, t43, pkin(1) * t65 + t38 * t47, -pkin(1) * t55 * t70 - t39 * t47 (-t38 * t55 + t39 * t57) * t69, t39 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t61, t37 ^ 2 / 0.2e1, -t37 * t35, t37 * t47, t35 ^ 2 / 0.2e1, -t35 * t47, t43, t15 * t47 + t40 * t35, -t16 * t47 + t40 * t37, -t15 * t37 - t16 * t35, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t40 ^ 2 / 0.2e1, t66, -t73, t72, t67, -t71, t33, t13 * t25 + t9 * t34, -t10 * t34 + t13 * t27, -t10 * t25 - t9 * t27, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1, t33, -t72, t71, t66, -t73, t67, t7 * t25 + t6 * t27, -t8 * t25 + t6 * t34, -t8 * t27 - t7 * t34, t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t19 ^ 2 / 0.2e1, -t19 * t17, t19 * t24, t17 ^ 2 / 0.2e1, -t17 * t24, t24 ^ 2 / 0.2e1, t1 * t24 + t4 * t17, t4 * t19 - t2 * t24, -t1 * t19 - t2 * t17, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1;];
T_reg  = t11;
