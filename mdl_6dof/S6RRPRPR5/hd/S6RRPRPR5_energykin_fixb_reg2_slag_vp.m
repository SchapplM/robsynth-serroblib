% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR5_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:35:10
% EndTime: 2019-03-09 10:35:10
% DurationCPUTime: 0.25s
% Computational Cost: add. (1408->73), mult. (3892->165), div. (0->0), fcn. (3104->12), ass. (0->57)
t64 = cos(qJ(2));
t72 = cos(pkin(6)) * qJD(1);
t69 = pkin(1) * t72;
t53 = t64 * t69;
t54 = qJD(2) + t72;
t63 = sin(qJ(2));
t58 = sin(pkin(6));
t73 = qJD(1) * t58;
t68 = t63 * t73;
t36 = t54 * pkin(2) + t53 + (-pkin(8) - qJ(3)) * t68;
t67 = t64 * t73;
t46 = pkin(8) * t67 + t63 * t69;
t40 = qJ(3) * t67 + t46;
t57 = sin(pkin(11));
t59 = cos(pkin(11));
t25 = t36 * t57 + t40 * t59;
t22 = pkin(9) * t54 + t25;
t42 = t57 * t68 - t59 * t67;
t44 = (t57 * t64 + t59 * t63) * t73;
t47 = qJD(3) + (-pkin(2) * t64 - pkin(1)) * t73;
t30 = t42 * pkin(3) - t44 * pkin(9) + t47;
t62 = sin(qJ(4));
t77 = cos(qJ(4));
t15 = t22 * t77 + t30 * t62;
t41 = qJD(4) + t42;
t10 = qJ(5) * t41 + t15;
t24 = t36 * t59 - t40 * t57;
t21 = -pkin(3) * t54 - t24;
t33 = t44 * t62 - t54 * t77;
t35 = t44 * t77 + t54 * t62;
t13 = pkin(4) * t33 - qJ(5) * t35 + t21;
t56 = sin(pkin(12));
t74 = cos(pkin(12));
t6 = t10 * t74 + t13 * t56;
t76 = cos(qJ(6));
t65 = qJD(1) ^ 2;
t75 = t58 ^ 2 * t65;
t71 = t33 ^ 2 / 0.2e1;
t70 = t64 * t75;
t66 = t75 / 0.2e1;
t5 = -t10 * t56 + t13 * t74;
t14 = -t22 * t62 + t30 * t77;
t9 = -pkin(4) * t41 + qJD(5) - t14;
t61 = sin(qJ(6));
t50 = t54 ^ 2 / 0.2e1;
t45 = -pkin(8) * t68 + t53;
t32 = qJD(6) + t33;
t28 = t35 * t74 + t41 * t56;
t26 = t35 * t56 - t41 * t74;
t18 = -t26 * t61 + t28 * t76;
t16 = t26 * t76 + t28 * t61;
t7 = pkin(5) * t26 + t9;
t4 = -pkin(10) * t26 + t6;
t3 = pkin(5) * t33 - pkin(10) * t28 + t5;
t2 = t3 * t61 + t4 * t76;
t1 = t3 * t76 - t4 * t61;
t8 = [0, 0, 0, 0, 0, t65 / 0.2e1, 0, 0, 0, 0, t63 ^ 2 * t66, t63 * t70, t54 * t68, t64 ^ 2 * t66, t54 * t67, t50, pkin(1) * t70 + t45 * t54, -pkin(1) * t63 * t75 - t46 * t54 (-t45 * t63 + t46 * t64) * t73, t46 ^ 2 / 0.2e1 + t45 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t66, t44 ^ 2 / 0.2e1, -t44 * t42, t44 * t54, t42 ^ 2 / 0.2e1, -t42 * t54, t50, t24 * t54 + t42 * t47, -t25 * t54 + t44 * t47, -t24 * t44 - t25 * t42, t25 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1 + t47 ^ 2 / 0.2e1, t35 ^ 2 / 0.2e1, -t35 * t33, t35 * t41, t71, -t33 * t41, t41 ^ 2 / 0.2e1, t14 * t41 + t21 * t33, -t15 * t41 + t21 * t35, -t14 * t35 - t15 * t33, t15 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t28 ^ 2 / 0.2e1, -t28 * t26, t28 * t33, t26 ^ 2 / 0.2e1, -t26 * t33, t71, t26 * t9 + t33 * t5, t28 * t9 - t33 * t6, -t26 * t6 - t28 * t5, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t32, t16 ^ 2 / 0.2e1, -t16 * t32, t32 ^ 2 / 0.2e1, t1 * t32 + t16 * t7, t18 * t7 - t2 * t32, -t1 * t18 - t16 * t2, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
