% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPR8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR8_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:44:47
% EndTime: 2019-03-09 22:44:47
% DurationCPUTime: 0.22s
% Computational Cost: add. (814->66), mult. (1732->144), div. (0->0), fcn. (1181->8), ass. (0->55)
t54 = qJD(1) ^ 2;
t72 = t54 / 0.2e1;
t71 = pkin(4) + pkin(5);
t70 = cos(qJ(3));
t69 = cos(qJ(4));
t68 = cos(qJ(6));
t51 = sin(qJ(3));
t52 = sin(qJ(2));
t64 = qJD(1) * t52;
t32 = -t70 * qJD(2) + t51 * t64;
t34 = t51 * qJD(2) + t70 * t64;
t50 = sin(qJ(4));
t22 = t69 * t32 + t50 * t34;
t24 = -t50 * t32 + t69 * t34;
t67 = t24 * t22;
t53 = cos(qJ(2));
t63 = t53 * qJD(1);
t43 = -qJD(3) + t63;
t41 = -qJD(4) + t43;
t66 = t41 * t22;
t65 = t53 * t54;
t31 = (-pkin(2) * t53 - pkin(8) * t52 - pkin(1)) * qJD(1);
t40 = pkin(7) * t63 + qJD(2) * pkin(8);
t25 = t70 * t31 - t51 * t40;
t15 = -t43 * pkin(3) - t34 * pkin(9) + t25;
t26 = t51 * t31 + t70 * t40;
t18 = -t32 * pkin(9) + t26;
t9 = t50 * t15 + t69 * t18;
t62 = t22 ^ 2 / 0.2e1;
t61 = qJD(1) * qJD(2);
t7 = -t41 * qJ(5) + t9;
t60 = t52 * t61;
t59 = t53 * t61;
t39 = -qJD(2) * pkin(2) + pkin(7) * t64;
t8 = t69 * t15 - t50 * t18;
t58 = -t32 * pkin(3) - t39;
t57 = qJD(5) - t8;
t56 = t24 * qJ(5) + t58;
t49 = sin(qJ(6));
t48 = t53 ^ 2;
t47 = t52 ^ 2;
t38 = qJD(6) + t41;
t36 = t41 ^ 2 / 0.2e1;
t21 = t24 ^ 2 / 0.2e1;
t19 = t24 * t41;
t13 = t49 * t22 + t68 * t24;
t11 = -t68 * t22 + t49 * t24;
t10 = t22 * pkin(4) - t56;
t6 = t41 * pkin(4) + t57;
t5 = -t71 * t22 + t56;
t4 = t22 * pkin(10) + t7;
t3 = -t24 * pkin(10) + t71 * t41 + t57;
t2 = t49 * t3 + t68 * t4;
t1 = t68 * t3 - t49 * t4;
t12 = [0, 0, 0, 0, 0, t72, 0, 0, 0, 0, t47 * t72, t52 * t65, t60, t48 * t72, t59, qJD(2) ^ 2 / 0.2e1, pkin(1) * t65 - pkin(7) * t60, -t54 * pkin(1) * t52 - pkin(7) * t59 (t47 + t48) * t54 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t48 / 0.2e1 + t47 / 0.2e1) * pkin(7) ^ 2) * t54, t34 ^ 2 / 0.2e1, -t34 * t32, -t34 * t43, t32 ^ 2 / 0.2e1, t32 * t43, t43 ^ 2 / 0.2e1, -t25 * t43 + t39 * t32, t26 * t43 + t39 * t34, -t25 * t34 - t26 * t32, t26 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1 + t39 ^ 2 / 0.2e1, t21, -t67, -t19, t62, t66, t36, -t22 * t58 - t8 * t41, -t24 * t58 + t9 * t41, -t9 * t22 - t8 * t24, t9 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t58 ^ 2 / 0.2e1, t21, -t19, t67, t36, -t66, t62, t10 * t22 + t6 * t41, -t7 * t22 + t6 * t24, -t10 * t24 - t7 * t41, t7 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t13 * t11, t13 * t38, t11 ^ 2 / 0.2e1, -t11 * t38, t38 ^ 2 / 0.2e1, t1 * t38 + t5 * t11, t5 * t13 - t2 * t38, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t12;
