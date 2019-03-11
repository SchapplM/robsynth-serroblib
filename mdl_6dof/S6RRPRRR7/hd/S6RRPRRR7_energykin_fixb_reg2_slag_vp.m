% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRR7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR7_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:59:44
% EndTime: 2019-03-09 13:59:44
% DurationCPUTime: 0.22s
% Computational Cost: add. (677->66), mult. (1388->149), div. (0->0), fcn. (865->8), ass. (0->51)
t55 = sin(qJ(4));
t56 = sin(qJ(2));
t57 = cos(qJ(4));
t58 = cos(qJ(2));
t30 = (t55 * t56 + t57 * t58) * qJD(1);
t59 = qJD(1) ^ 2;
t71 = t59 / 0.2e1;
t70 = cos(qJ(5));
t69 = cos(qJ(6));
t68 = t58 * t59;
t66 = qJD(1) * t58;
t67 = qJD(1) * t56;
t34 = -qJD(1) * pkin(1) - pkin(2) * t66 - qJ(3) * t67;
t24 = pkin(3) * t66 - t34;
t32 = (-t55 * t58 + t56 * t57) * qJD(1);
t13 = t30 * pkin(4) - t32 * pkin(9) + t24;
t65 = pkin(7) * t67 + qJD(3);
t25 = -pkin(8) * t67 + (-pkin(2) - pkin(3)) * qJD(2) + t65;
t36 = pkin(7) * t66 + qJD(2) * qJ(3);
t33 = -pkin(8) * t66 + t36;
t18 = t55 * t25 + t57 * t33;
t46 = qJD(2) - qJD(4);
t16 = -t46 * pkin(9) + t18;
t54 = sin(qJ(5));
t6 = t54 * t13 + t70 * t16;
t64 = qJD(1) * qJD(2);
t63 = t56 * t68;
t40 = t56 * t64;
t62 = t58 * t64;
t5 = t70 * t13 - t54 * t16;
t17 = t57 * t25 - t55 * t33;
t15 = t46 * pkin(4) - t17;
t29 = qJD(5) + t30;
t53 = sin(qJ(6));
t51 = t58 ^ 2;
t50 = t56 ^ 2;
t48 = qJD(2) ^ 2 / 0.2e1;
t39 = t51 * t71;
t38 = t50 * t71;
t35 = -qJD(2) * pkin(2) + t65;
t28 = qJD(6) + t29;
t22 = t70 * t32 - t54 * t46;
t20 = t54 * t32 + t70 * t46;
t10 = -t53 * t20 + t69 * t22;
t8 = t69 * t20 + t53 * t22;
t7 = t20 * pkin(5) + t15;
t4 = -t20 * pkin(10) + t6;
t3 = t29 * pkin(5) - t22 * pkin(10) + t5;
t2 = t53 * t3 + t69 * t4;
t1 = t69 * t3 - t53 * t4;
t9 = [0, 0, 0, 0, 0, t71, 0, 0, 0, 0, t38, t63, t40, t39, t62, t48, pkin(1) * t68 - pkin(7) * t40, -t59 * pkin(1) * t56 - pkin(7) * t62 (t50 + t51) * t59 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t51 / 0.2e1 + t50 / 0.2e1) * pkin(7) ^ 2) * t59, t38, t40, -t63, t48, -t62, t39, -t35 * qJD(2) - t34 * t66 (t35 * t56 + t36 * t58) * qJD(1), t36 * qJD(2) - t34 * t67, t36 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1 + t35 ^ 2 / 0.2e1, t32 ^ 2 / 0.2e1, -t32 * t30, -t32 * t46, t30 ^ 2 / 0.2e1, t30 * t46, t46 ^ 2 / 0.2e1, -t17 * t46 + t24 * t30, t18 * t46 + t24 * t32, -t17 * t32 - t18 * t30, t18 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * t29, t20 ^ 2 / 0.2e1, -t20 * t29, t29 ^ 2 / 0.2e1, t15 * t20 + t5 * t29, t15 * t22 - t6 * t29, -t6 * t20 - t5 * t22, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, t10 * t28, t8 ^ 2 / 0.2e1, -t8 * t28, t28 ^ 2 / 0.2e1, t1 * t28 + t7 * t8, t7 * t10 - t2 * t28, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t9;
