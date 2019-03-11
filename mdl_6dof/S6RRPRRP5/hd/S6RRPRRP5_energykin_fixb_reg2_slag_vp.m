% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRP5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP5_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:05:24
% EndTime: 2019-03-09 12:05:24
% DurationCPUTime: 0.23s
% Computational Cost: add. (1061->69), mult. (2964->148), div. (0->0), fcn. (2323->10), ass. (0->56)
t63 = cos(qJ(2));
t70 = cos(pkin(6)) * qJD(1);
t68 = pkin(1) * t70;
t53 = t63 * t68;
t54 = qJD(2) + t70;
t62 = sin(qJ(2));
t57 = sin(pkin(6));
t71 = qJD(1) * t57;
t67 = t62 * t71;
t36 = t54 * pkin(2) + t53 + (-pkin(8) - qJ(3)) * t67;
t66 = t63 * t71;
t46 = pkin(8) * t66 + t62 * t68;
t39 = qJ(3) * t66 + t46;
t56 = sin(pkin(11));
t58 = cos(pkin(11));
t21 = t58 * t36 - t56 * t39;
t19 = -t54 * pkin(3) - t21;
t44 = (t56 * t63 + t58 * t62) * t71;
t61 = sin(qJ(4));
t74 = cos(qJ(4));
t33 = t61 * t44 - t74 * t54;
t35 = t74 * t44 + t61 * t54;
t11 = t33 * pkin(4) - t35 * pkin(10) + t19;
t60 = sin(qJ(5));
t73 = cos(qJ(5));
t22 = t56 * t36 + t58 * t39;
t20 = t54 * pkin(9) + t22;
t42 = t56 * t67 - t58 * t66;
t47 = qJD(3) + (-pkin(2) * t63 - pkin(1)) * t71;
t29 = t42 * pkin(3) - t44 * pkin(9) + t47;
t13 = t74 * t20 + t61 * t29;
t41 = qJD(4) + t42;
t8 = t41 * pkin(10) + t13;
t4 = t60 * t11 + t73 * t8;
t64 = qJD(1) ^ 2;
t72 = t57 ^ 2 * t64;
t69 = t63 * t72;
t65 = t72 / 0.2e1;
t3 = t73 * t11 - t60 * t8;
t12 = -t61 * t20 + t74 * t29;
t7 = -t41 * pkin(4) - t12;
t50 = t54 ^ 2 / 0.2e1;
t45 = -pkin(8) * t67 + t53;
t32 = qJD(5) + t33;
t30 = t32 ^ 2 / 0.2e1;
t27 = t73 * t35 + t60 * t41;
t25 = t60 * t35 - t73 * t41;
t24 = t27 ^ 2 / 0.2e1;
t23 = t25 ^ 2 / 0.2e1;
t16 = t27 * t32;
t15 = t25 * t32;
t14 = t27 * t25;
t5 = t25 * pkin(5) + qJD(6) + t7;
t2 = -t25 * qJ(6) + t4;
t1 = t32 * pkin(5) - t27 * qJ(6) + t3;
t6 = [0, 0, 0, 0, 0, t64 / 0.2e1, 0, 0, 0, 0, t62 ^ 2 * t65, t62 * t69, t54 * t67, t63 ^ 2 * t65, t54 * t66, t50, pkin(1) * t69 + t45 * t54, -pkin(1) * t62 * t72 - t46 * t54 (-t45 * t62 + t46 * t63) * t71, t46 ^ 2 / 0.2e1 + t45 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t65, t44 ^ 2 / 0.2e1, -t44 * t42, t44 * t54, t42 ^ 2 / 0.2e1, -t42 * t54, t50, t21 * t54 + t47 * t42, -t22 * t54 + t47 * t44, -t21 * t44 - t22 * t42, t22 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1 + t47 ^ 2 / 0.2e1, t35 ^ 2 / 0.2e1, -t35 * t33, t35 * t41, t33 ^ 2 / 0.2e1, -t33 * t41, t41 ^ 2 / 0.2e1, t12 * t41 + t19 * t33, -t13 * t41 + t19 * t35, -t12 * t35 - t13 * t33, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t24, -t14, t16, t23, -t15, t30, t7 * t25 + t3 * t32, t7 * t27 - t4 * t32, -t4 * t25 - t3 * t27, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t24, -t14, t16, t23, -t15, t30, t1 * t32 + t5 * t25, -t2 * t32 + t5 * t27, -t1 * t27 - t2 * t25, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t6;
