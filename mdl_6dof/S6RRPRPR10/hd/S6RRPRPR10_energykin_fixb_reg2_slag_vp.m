% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRPR10
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
% Datum: 2019-03-09 11:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPR10_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR10_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:09:23
% EndTime: 2019-03-09 11:09:23
% DurationCPUTime: 0.24s
% Computational Cost: add. (1075->69), mult. (2639->147), div. (0->0), fcn. (2044->10), ass. (0->58)
t76 = pkin(4) + pkin(10);
t75 = cos(qJ(6));
t68 = cos(pkin(6)) * qJD(1);
t48 = qJD(2) + t68;
t50 = sin(pkin(11));
t55 = sin(qJ(2));
t51 = sin(pkin(6));
t69 = qJD(1) * t51;
t63 = t55 * t69;
t70 = cos(pkin(11));
t35 = -t70 * t48 + t50 * t63;
t37 = t50 * t48 + t70 * t63;
t54 = sin(qJ(4));
t56 = cos(qJ(4));
t24 = t56 * t35 + t54 * t37;
t26 = -t54 * t35 + t56 * t37;
t74 = t26 * t24;
t57 = cos(qJ(2));
t62 = t57 * t69;
t42 = -qJD(4) + t62;
t73 = t26 * t42;
t72 = t42 * t24;
t58 = qJD(1) ^ 2;
t71 = t51 ^ 2 * t58;
t64 = pkin(1) * t68;
t39 = pkin(8) * t62 + t55 * t64;
t32 = t48 * qJ(3) + t39;
t34 = (-pkin(2) * t57 - qJ(3) * t55 - pkin(1)) * t69;
t19 = -t50 * t32 + t70 * t34;
t12 = -pkin(3) * t62 - t37 * pkin(9) + t19;
t20 = t70 * t32 + t50 * t34;
t15 = -t35 * pkin(9) + t20;
t9 = t54 * t12 + t56 * t15;
t67 = t24 ^ 2 / 0.2e1;
t66 = t26 ^ 2 / 0.2e1;
t65 = t57 * t71;
t61 = t71 / 0.2e1;
t8 = t56 * t12 - t54 * t15;
t7 = t42 * qJ(5) - t9;
t60 = qJD(5) - t8;
t38 = -pkin(8) * t63 + t57 * t64;
t29 = -t48 * pkin(2) + qJD(3) - t38;
t23 = t35 * pkin(3) + t29;
t59 = -t26 * qJ(5) + t23;
t53 = sin(qJ(6));
t44 = t57 ^ 2 * t61;
t40 = t42 ^ 2 / 0.2e1;
t22 = qJD(6) + t26;
t18 = t53 * t24 - t75 * t42;
t16 = -t75 * t24 - t53 * t42;
t10 = t24 * pkin(4) + t59;
t6 = t42 * pkin(4) + t60;
t5 = t76 * t24 + t59;
t4 = -t24 * pkin(5) - t7;
t3 = t26 * pkin(5) + t76 * t42 + t60;
t2 = t53 * t3 + t75 * t5;
t1 = t75 * t3 - t53 * t5;
t11 = [0, 0, 0, 0, 0, t58 / 0.2e1, 0, 0, 0, 0, t55 ^ 2 * t61, t55 * t65, t48 * t63, t44, t48 * t62, t48 ^ 2 / 0.2e1, pkin(1) * t65 + t38 * t48, -pkin(1) * t55 * t71 - t39 * t48 (-t38 * t55 + t39 * t57) * t69, t39 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t61, t37 ^ 2 / 0.2e1, -t37 * t35, -t37 * t62, t35 ^ 2 / 0.2e1, t35 * t62, t44, -t19 * t62 + t29 * t35, t20 * t62 + t29 * t37, -t19 * t37 - t20 * t35, t20 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1, t66, -t74, -t73, t67, t72, t40, t23 * t24 - t8 * t42, t23 * t26 + t9 * t42, -t9 * t24 - t8 * t26, t9 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1, t40, t73, -t72, t66, -t74, t67, t7 * t24 + t6 * t26, -t10 * t24 - t6 * t42, -t10 * t26 + t7 * t42, t10 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t22, t16 ^ 2 / 0.2e1, -t16 * t22, t22 ^ 2 / 0.2e1, t1 * t22 + t4 * t16, t4 * t18 - t2 * t22, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1;];
T_reg  = t11;
