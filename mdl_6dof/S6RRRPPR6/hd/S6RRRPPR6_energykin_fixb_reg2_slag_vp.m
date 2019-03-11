% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPPR6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR6_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:53:53
% EndTime: 2019-03-09 15:53:53
% DurationCPUTime: 0.24s
% Computational Cost: add. (1096->69), mult. (2637->148), div. (0->0), fcn. (2044->10), ass. (0->57)
t75 = pkin(4) + pkin(10);
t74 = cos(qJ(3));
t73 = cos(qJ(6));
t67 = cos(pkin(6)) * qJD(1);
t47 = qJD(2) + t67;
t54 = sin(qJ(3));
t55 = sin(qJ(2));
t50 = sin(pkin(6));
t68 = qJD(1) * t50;
t62 = t55 * t68;
t35 = -t74 * t47 + t54 * t62;
t37 = t54 * t47 + t74 * t62;
t49 = sin(pkin(11));
t51 = cos(pkin(11));
t24 = t51 * t35 + t49 * t37;
t26 = -t49 * t35 + t51 * t37;
t72 = t26 * t24;
t56 = cos(qJ(2));
t61 = t56 * t68;
t42 = -qJD(3) + t61;
t71 = t26 * t42;
t70 = t42 * t24;
t57 = qJD(1) ^ 2;
t69 = t50 ^ 2 * t57;
t63 = pkin(1) * t67;
t39 = pkin(8) * t61 + t55 * t63;
t32 = t47 * pkin(9) + t39;
t34 = (-pkin(2) * t56 - pkin(9) * t55 - pkin(1)) * t68;
t19 = -t54 * t32 + t74 * t34;
t12 = -t42 * pkin(3) - t37 * qJ(4) + t19;
t20 = t74 * t32 + t54 * t34;
t15 = -t35 * qJ(4) + t20;
t9 = t49 * t12 + t51 * t15;
t66 = t24 ^ 2 / 0.2e1;
t65 = t26 ^ 2 / 0.2e1;
t64 = t56 * t69;
t60 = t69 / 0.2e1;
t8 = t51 * t12 - t49 * t15;
t7 = t42 * qJ(5) - t9;
t59 = qJD(5) - t8;
t38 = -pkin(8) * t62 + t56 * t63;
t31 = -t47 * pkin(2) - t38;
t23 = t35 * pkin(3) + qJD(4) + t31;
t58 = -t26 * qJ(5) + t23;
t53 = sin(qJ(6));
t40 = t42 ^ 2 / 0.2e1;
t22 = qJD(6) + t26;
t18 = t53 * t24 - t73 * t42;
t16 = -t73 * t24 - t53 * t42;
t10 = t24 * pkin(4) + t58;
t6 = t42 * pkin(4) + t59;
t5 = t75 * t24 + t58;
t4 = -t24 * pkin(5) - t7;
t3 = t26 * pkin(5) + t75 * t42 + t59;
t2 = t53 * t3 + t73 * t5;
t1 = t73 * t3 - t53 * t5;
t11 = [0, 0, 0, 0, 0, t57 / 0.2e1, 0, 0, 0, 0, t55 ^ 2 * t60, t55 * t64, t47 * t62, t56 ^ 2 * t60, t47 * t61, t47 ^ 2 / 0.2e1, pkin(1) * t64 + t38 * t47, -pkin(1) * t55 * t69 - t39 * t47 (-t38 * t55 + t39 * t56) * t68, t39 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t60, t37 ^ 2 / 0.2e1, -t37 * t35, -t37 * t42, t35 ^ 2 / 0.2e1, t35 * t42, t40, -t19 * t42 + t31 * t35, t20 * t42 + t31 * t37, -t19 * t37 - t20 * t35, t20 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1, t65, -t72, -t71, t66, t70, t40, t23 * t24 - t8 * t42, t23 * t26 + t9 * t42, -t9 * t24 - t8 * t26, t9 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1, t40, t71, -t70, t65, -t72, t66, t7 * t24 + t6 * t26, -t10 * t24 - t6 * t42, -t10 * t26 + t7 * t42, t10 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t22, t16 ^ 2 / 0.2e1, -t16 * t22, t22 ^ 2 / 0.2e1, t1 * t22 + t4 * t16, t4 * t18 - t2 * t22, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1;];
T_reg  = t11;
