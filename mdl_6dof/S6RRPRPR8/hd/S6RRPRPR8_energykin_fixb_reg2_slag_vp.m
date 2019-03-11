% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPR8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR8_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:53:18
% EndTime: 2019-03-09 10:53:19
% DurationCPUTime: 0.22s
% Computational Cost: add. (760->66), mult. (1733->142), div. (0->0), fcn. (1181->8), ass. (0->55)
t53 = qJD(1) ^ 2;
t72 = t53 / 0.2e1;
t71 = pkin(4) + pkin(5);
t70 = cos(qJ(4));
t69 = cos(qJ(6));
t48 = sin(pkin(10));
t51 = sin(qJ(2));
t64 = qJD(1) * t51;
t65 = cos(pkin(10));
t32 = -t65 * qJD(2) + t48 * t64;
t34 = t48 * qJD(2) + t65 * t64;
t50 = sin(qJ(4));
t22 = t70 * t32 + t50 * t34;
t24 = -t50 * t32 + t70 * t34;
t68 = t24 * t22;
t52 = cos(qJ(2));
t63 = t52 * qJD(1);
t41 = -qJD(4) + t63;
t67 = t41 * t22;
t66 = t52 * t53;
t31 = (-pkin(2) * t52 - qJ(3) * t51 - pkin(1)) * qJD(1);
t37 = pkin(7) * t63 + qJD(2) * qJ(3);
t25 = t65 * t31 - t48 * t37;
t15 = -pkin(3) * t63 - t34 * pkin(8) + t25;
t26 = t48 * t31 + t65 * t37;
t18 = -t32 * pkin(8) + t26;
t9 = t50 * t15 + t70 * t18;
t62 = t22 ^ 2 / 0.2e1;
t61 = qJD(1) * qJD(2);
t7 = -t41 * qJ(5) + t9;
t60 = t51 * t61;
t59 = t52 * t61;
t8 = t70 * t15 - t50 * t18;
t58 = qJD(2) * pkin(2) - pkin(7) * t64 - qJD(3);
t57 = qJD(5) - t8;
t56 = -t32 * pkin(3) + t58;
t55 = t24 * qJ(5) + t56;
t49 = sin(qJ(6));
t47 = t52 ^ 2;
t46 = t51 ^ 2;
t43 = t47 * t72;
t40 = qJD(6) + t41;
t38 = t41 ^ 2 / 0.2e1;
t21 = t24 ^ 2 / 0.2e1;
t19 = t24 * t41;
t13 = t49 * t22 + t69 * t24;
t11 = -t69 * t22 + t49 * t24;
t10 = t22 * pkin(4) - t55;
t6 = t41 * pkin(4) + t57;
t5 = -t71 * t22 + t55;
t4 = t22 * pkin(9) + t7;
t3 = -t24 * pkin(9) + t71 * t41 + t57;
t2 = t49 * t3 + t69 * t4;
t1 = t69 * t3 - t49 * t4;
t12 = [0, 0, 0, 0, 0, t72, 0, 0, 0, 0, t46 * t72, t51 * t66, t60, t43, t59, qJD(2) ^ 2 / 0.2e1, pkin(1) * t66 - pkin(7) * t60, -t53 * pkin(1) * t51 - pkin(7) * t59 (t46 + t47) * t53 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t47 / 0.2e1 + t46 / 0.2e1) * pkin(7) ^ 2) * t53, t34 ^ 2 / 0.2e1, -t34 * t32, -t34 * t63, t32 ^ 2 / 0.2e1, t32 * t63, t43, -t25 * t63 - t32 * t58, t26 * t63 - t34 * t58, -t25 * t34 - t26 * t32, t26 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1 + t58 ^ 2 / 0.2e1, t21, -t68, -t19, t62, t67, t38, -t22 * t56 - t8 * t41, -t24 * t56 + t9 * t41, -t9 * t22 - t8 * t24, t9 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t56 ^ 2 / 0.2e1, t21, -t19, t68, t38, -t67, t62, t10 * t22 + t6 * t41, -t7 * t22 + t6 * t24, -t10 * t24 - t7 * t41, t7 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t13 * t11, t13 * t40, t11 ^ 2 / 0.2e1, -t11 * t40, t40 ^ 2 / 0.2e1, t1 * t40 + t5 * t11, t5 * t13 - t2 * t40, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t12;
