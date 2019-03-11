% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRP11_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP11_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:49:01
% EndTime: 2019-03-09 17:49:02
% DurationCPUTime: 0.20s
% Computational Cost: add. (769->65), mult. (1814->130), div. (0->0), fcn. (1341->8), ass. (0->57)
t74 = pkin(3) + pkin(10);
t67 = cos(pkin(6)) * qJD(1);
t48 = qJD(2) + t67;
t53 = sin(qJ(3));
t55 = cos(qJ(3));
t54 = sin(qJ(2));
t50 = sin(pkin(6));
t68 = qJD(1) * t50;
t62 = t54 * t68;
t34 = -t55 * t48 + t53 * t62;
t56 = cos(qJ(2));
t63 = pkin(1) * t67;
t38 = -pkin(8) * t62 + t56 * t63;
t28 = -t48 * pkin(2) - t38;
t36 = t53 * t48 + t55 * t62;
t58 = -t36 * qJ(4) + t28;
t10 = t74 * t34 + t58;
t52 = sin(qJ(5));
t73 = cos(qJ(5));
t61 = t56 * t68;
t43 = -qJD(3) + t61;
t39 = pkin(8) * t61 + t54 * t63;
t29 = t48 * pkin(9) + t39;
t32 = (-pkin(2) * t56 - pkin(9) * t54 - pkin(1)) * t68;
t18 = -t53 * t29 + t55 * t32;
t59 = qJD(4) - t18;
t8 = t36 * pkin(4) + t74 * t43 + t59;
t4 = t73 * t10 + t52 * t8;
t72 = t36 * t34;
t71 = t36 * t43;
t70 = t43 * t34;
t57 = qJD(1) ^ 2;
t69 = t50 ^ 2 * t57;
t19 = t55 * t29 + t53 * t32;
t66 = t34 ^ 2 / 0.2e1;
t65 = t36 ^ 2 / 0.2e1;
t64 = t56 * t69;
t15 = t43 * qJ(4) - t19;
t60 = t69 / 0.2e1;
t3 = -t52 * t10 + t73 * t8;
t11 = -t34 * pkin(4) - t15;
t40 = t43 ^ 2 / 0.2e1;
t33 = qJD(5) + t36;
t30 = t33 ^ 2 / 0.2e1;
t24 = t52 * t34 - t73 * t43;
t22 = -t73 * t34 - t52 * t43;
t21 = t24 ^ 2 / 0.2e1;
t20 = t22 ^ 2 / 0.2e1;
t17 = t24 * t33;
t16 = t22 * t33;
t14 = t43 * pkin(3) + t59;
t13 = t34 * pkin(3) + t58;
t12 = t24 * t22;
t5 = t22 * pkin(5) + qJD(6) + t11;
t2 = -t22 * qJ(6) + t4;
t1 = t33 * pkin(5) - t24 * qJ(6) + t3;
t6 = [0, 0, 0, 0, 0, t57 / 0.2e1, 0, 0, 0, 0, t54 ^ 2 * t60, t54 * t64, t48 * t62, t56 ^ 2 * t60, t48 * t61, t48 ^ 2 / 0.2e1, pkin(1) * t64 + t38 * t48, -pkin(1) * t54 * t69 - t39 * t48 (-t38 * t54 + t39 * t56) * t68, t39 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t60, t65, -t72, -t71, t66, t70, t40, -t18 * t43 + t28 * t34, t19 * t43 + t28 * t36, -t18 * t36 - t19 * t34, t19 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1, t40, t71, -t70, t65, -t72, t66, t14 * t36 + t15 * t34, -t13 * t34 - t14 * t43, -t13 * t36 + t15 * t43, t13 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t21, -t12, t17, t20, -t16, t30, t11 * t22 + t3 * t33, t11 * t24 - t4 * t33, -t4 * t22 - t3 * t24, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1, t21, -t12, t17, t20, -t16, t30, t1 * t33 + t5 * t22, -t2 * t33 + t5 * t24, -t1 * t24 - t2 * t22, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t6;
