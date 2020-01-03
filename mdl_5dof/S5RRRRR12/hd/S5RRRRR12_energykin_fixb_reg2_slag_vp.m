% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRR12_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR12_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:55:10
% EndTime: 2019-12-31 22:55:10
% DurationCPUTime: 0.29s
% Computational Cost: add. (1190->64), mult. (3294->153), div. (0->0), fcn. (2653->12), ass. (0->53)
t67 = cos(pkin(5)) * qJD(1);
t47 = qJD(2) + t67;
t49 = sin(pkin(6));
t51 = cos(pkin(6));
t58 = cos(qJ(2));
t50 = sin(pkin(5));
t68 = qJD(1) * t50;
t63 = t58 * t68;
t75 = t47 * t49 + t51 * t63;
t56 = sin(qJ(2));
t65 = pkin(1) * t67;
t40 = pkin(8) * t63 + t56 * t65;
t27 = t75 * pkin(9) + t40;
t46 = t58 * t65;
t64 = t56 * t68;
t29 = t47 * pkin(2) + t46 + (-pkin(9) * t51 - pkin(8)) * t64;
t36 = (-pkin(9) * t49 * t56 - pkin(2) * t58 - pkin(1)) * t68;
t55 = sin(qJ(3));
t57 = cos(qJ(3));
t13 = -t55 * t27 + (t29 * t51 + t36 * t49) * t57;
t69 = t51 * t55;
t70 = t49 * t55;
t14 = t57 * t27 + t29 * t69 + t36 * t70;
t37 = -t51 * t47 + t49 * t63 - qJD(3);
t12 = -t37 * pkin(10) + t14;
t54 = sin(qJ(4));
t74 = cos(qJ(4));
t18 = -t49 * t29 + t51 * t36;
t30 = t55 * t64 - t75 * t57;
t32 = t47 * t70 + (t56 * t57 + t58 * t69) * t68;
t9 = t30 * pkin(3) - t32 * pkin(10) + t18;
t6 = t74 * t12 + t54 * t9;
t73 = cos(qJ(5));
t59 = qJD(1) ^ 2;
t71 = t50 ^ 2 * t59;
t66 = t58 * t71;
t62 = t71 / 0.2e1;
t20 = t54 * t32 + t74 * t37;
t5 = -t54 * t12 + t74 * t9;
t11 = t37 * pkin(3) - t13;
t53 = sin(qJ(5));
t39 = -pkin(8) * t64 + t46;
t28 = qJD(4) + t30;
t22 = t74 * t32 - t54 * t37;
t19 = qJD(5) + t20;
t17 = t73 * t22 + t53 * t28;
t15 = t53 * t22 - t73 * t28;
t7 = t20 * pkin(4) - t22 * pkin(11) + t11;
t4 = t28 * pkin(11) + t6;
t3 = -t28 * pkin(4) - t5;
t2 = t73 * t4 + t53 * t7;
t1 = -t53 * t4 + t73 * t7;
t8 = [0, 0, 0, 0, 0, t59 / 0.2e1, 0, 0, 0, 0, t56 ^ 2 * t62, t56 * t66, t47 * t64, t58 ^ 2 * t62, t47 * t63, t47 ^ 2 / 0.2e1, pkin(1) * t66 + t39 * t47, -pkin(1) * t56 * t71 - t40 * t47, (-t39 * t56 + t40 * t58) * t68, t40 ^ 2 / 0.2e1 + t39 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t62, t32 ^ 2 / 0.2e1, -t32 * t30, -t32 * t37, t30 ^ 2 / 0.2e1, t30 * t37, t37 ^ 2 / 0.2e1, -t13 * t37 + t18 * t30, t14 * t37 + t18 * t32, -t13 * t32 - t14 * t30, t14 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * t28, t20 ^ 2 / 0.2e1, -t20 * t28, t28 ^ 2 / 0.2e1, t11 * t20 + t5 * t28, t11 * t22 - t6 * t28, -t6 * t20 - t5 * t22, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1, t17 ^ 2 / 0.2e1, -t17 * t15, t17 * t19, t15 ^ 2 / 0.2e1, -t15 * t19, t19 ^ 2 / 0.2e1, t1 * t19 + t3 * t15, t3 * t17 - t2 * t19, -t1 * t17 - t2 * t15, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t8;
