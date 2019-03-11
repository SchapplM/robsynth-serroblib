% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRP10_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP10_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:27:24
% EndTime: 2019-03-10 02:27:25
% DurationCPUTime: 0.25s
% Computational Cost: add. (1296->66), mult. (2974->150), div. (0->0), fcn. (2321->10), ass. (0->55)
t55 = sin(qJ(5));
t60 = cos(qJ(2));
t58 = sin(qJ(2));
t53 = sin(pkin(6));
t69 = qJD(1) * t53;
t64 = t58 * t69;
t68 = cos(pkin(6)) * qJD(1);
t65 = pkin(1) * t68;
t42 = -pkin(8) * t64 + t60 * t65;
t51 = qJD(2) + t68;
t33 = -t51 * pkin(2) - t42;
t57 = sin(qJ(3));
t59 = cos(qJ(3));
t39 = -t59 * t51 + t57 * t64;
t41 = t57 * t51 + t59 * t64;
t20 = t39 * pkin(3) - t41 * pkin(10) + t33;
t63 = t60 * t69;
t43 = pkin(8) * t63 + t58 * t65;
t34 = t51 * pkin(9) + t43;
t37 = (-pkin(2) * t60 - pkin(9) * t58 - pkin(1)) * t69;
t25 = t59 * t34 + t57 * t37;
t46 = -qJD(3) + t63;
t23 = -t46 * pkin(10) + t25;
t56 = sin(qJ(4));
t74 = cos(qJ(4));
t10 = t74 * t20 - t56 * t23;
t29 = t74 * t41 - t56 * t46;
t38 = qJD(4) + t39;
t7 = t38 * pkin(4) - t29 * pkin(11) + t10;
t73 = cos(qJ(5));
t11 = t56 * t20 + t74 * t23;
t27 = t56 * t41 + t74 * t46;
t9 = -t27 * pkin(11) + t11;
t4 = t55 * t7 + t73 * t9;
t15 = t73 * t27 + t55 * t29;
t17 = -t55 * t27 + t73 * t29;
t72 = t17 * t15;
t36 = qJD(5) + t38;
t71 = t36 * t15;
t61 = qJD(1) ^ 2;
t70 = t53 ^ 2 * t61;
t67 = t15 ^ 2 / 0.2e1;
t66 = t60 * t70;
t62 = t70 / 0.2e1;
t24 = -t57 * t34 + t59 * t37;
t22 = t46 * pkin(3) - t24;
t3 = -t55 * t9 + t73 * t7;
t13 = t27 * pkin(4) + t22;
t35 = t36 ^ 2 / 0.2e1;
t14 = t17 ^ 2 / 0.2e1;
t12 = t17 * t36;
t5 = t15 * pkin(5) - t17 * qJ(6) + t13;
t2 = t36 * qJ(6) + t4;
t1 = -t36 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t61 / 0.2e1, 0, 0, 0, 0, t58 ^ 2 * t62, t58 * t66, t51 * t64, t60 ^ 2 * t62, t51 * t63, t51 ^ 2 / 0.2e1, pkin(1) * t66 + t42 * t51, -pkin(1) * t58 * t70 - t43 * t51 (-t42 * t58 + t43 * t60) * t69, t43 ^ 2 / 0.2e1 + t42 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t62, t41 ^ 2 / 0.2e1, -t41 * t39, -t41 * t46, t39 ^ 2 / 0.2e1, t39 * t46, t46 ^ 2 / 0.2e1, -t24 * t46 + t33 * t39, t25 * t46 + t33 * t41, -t24 * t41 - t25 * t39, t25 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1, t29 ^ 2 / 0.2e1, -t29 * t27, t29 * t38, t27 ^ 2 / 0.2e1, -t27 * t38, t38 ^ 2 / 0.2e1, t10 * t38 + t22 * t27, -t11 * t38 + t22 * t29, -t10 * t29 - t11 * t27, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t14, -t72, t12, t67, -t71, t35, t13 * t15 + t3 * t36, t13 * t17 - t4 * t36, -t4 * t15 - t3 * t17, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1, t14, t12, t72, t35, t71, t67, -t1 * t36 + t5 * t15, t1 * t17 - t2 * t15, -t5 * t17 + t2 * t36, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
