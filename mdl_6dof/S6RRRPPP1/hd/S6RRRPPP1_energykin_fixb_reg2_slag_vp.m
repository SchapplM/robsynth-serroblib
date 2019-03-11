% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPPP1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPP1_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:18:27
% EndTime: 2019-03-09 15:18:28
% DurationCPUTime: 0.24s
% Computational Cost: add. (1006->69), mult. (2286->141), div. (0->0), fcn. (1666->8), ass. (0->55)
t54 = qJD(1) ^ 2;
t72 = t54 / 0.2e1;
t50 = sin(qJ(3));
t52 = cos(qJ(3));
t51 = sin(qJ(2));
t65 = qJD(1) * t51;
t36 = -t52 * qJD(2) + t50 * t65;
t38 = t50 * qJD(2) + t52 * t65;
t53 = cos(qJ(2));
t64 = t53 * qJD(1);
t43 = qJD(3) - t64;
t48 = sin(pkin(10));
t66 = cos(pkin(10));
t67 = cos(pkin(6));
t59 = t67 * t66;
t49 = sin(pkin(6));
t60 = t49 * t66;
t19 = t36 * t59 + t48 * t38 - t43 * t60;
t58 = t67 * t36 - t43 * t49;
t21 = t66 * t38 - t58 * t48;
t71 = t19 * t21;
t29 = -t49 * t36 - t67 * t43;
t10 = t21 * t29;
t11 = t29 * t19;
t70 = t53 * t54;
t69 = pkin(4) + qJ(6);
t34 = (-pkin(2) * t53 - pkin(9) * t51 - pkin(1)) * qJD(1);
t42 = pkin(8) * t64 + qJD(2) * pkin(9);
t27 = t50 * t34 + t52 * t42;
t68 = qJ(4) * t38;
t17 = t19 ^ 2 / 0.2e1;
t18 = t21 ^ 2 / 0.2e1;
t63 = qJD(1) * qJD(2);
t15 = -t58 * qJ(4) + t27;
t26 = t52 * t34 - t50 * t42;
t16 = t43 * pkin(3) - t67 * t68 + t26;
t41 = -qJD(2) * pkin(2) + pkin(8) * t65;
t25 = t36 * pkin(3) - t49 * t68 + t41;
t8 = t66 * t15 + (t16 * t67 + t25 * t49) * t48;
t62 = t51 * t63;
t61 = t53 * t63;
t9 = -t49 * t16 + t67 * t25 + qJD(4);
t6 = t29 * qJ(5) - t8;
t57 = -t21 * qJ(5) + t9;
t7 = -t48 * t15 + t16 * t59 + t25 * t60;
t56 = qJD(5) - t7;
t47 = t53 ^ 2;
t46 = t51 ^ 2;
t28 = t29 ^ 2 / 0.2e1;
t5 = t29 * pkin(4) + t56;
t4 = t19 * pkin(4) + t57;
t3 = -t19 * pkin(5) + qJD(6) - t6;
t2 = t69 * t19 + t57;
t1 = t21 * pkin(5) + t69 * t29 + t56;
t12 = [0, 0, 0, 0, 0, t72, 0, 0, 0, 0, t46 * t72, t51 * t70, t62, t47 * t72, t61, qJD(2) ^ 2 / 0.2e1, pkin(1) * t70 - pkin(8) * t62, -t54 * pkin(1) * t51 - pkin(8) * t61 (t46 + t47) * t54 * pkin(8) (pkin(1) ^ 2 / 0.2e1 + (t47 / 0.2e1 + t46 / 0.2e1) * pkin(8) ^ 2) * t54, t38 ^ 2 / 0.2e1, -t38 * t36, t38 * t43, t36 ^ 2 / 0.2e1, -t36 * t43, t43 ^ 2 / 0.2e1, t26 * t43 + t41 * t36, -t27 * t43 + t41 * t38, -t26 * t38 - t27 * t36, t27 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1 + t41 ^ 2 / 0.2e1, t18, -t71, -t10, t17, t11, t28, t9 * t19 - t7 * t29, t9 * t21 + t8 * t29, -t8 * t19 - t7 * t21, t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t28, t10, -t11, t18, -t71, t17, t6 * t19 + t5 * t21, -t4 * t19 - t5 * t29, -t4 * t21 + t6 * t29, t4 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1, t28, -t11, -t10, t17, t71, t18, t1 * t21 - t3 * t19, -t2 * t21 - t3 * t29, t1 * t29 + t2 * t19, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t12;
