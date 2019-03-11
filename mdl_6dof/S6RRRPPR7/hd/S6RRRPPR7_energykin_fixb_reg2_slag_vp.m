% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPPR7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR7_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:01:04
% EndTime: 2019-03-09 16:01:04
% DurationCPUTime: 0.18s
% Computational Cost: add. (737->66), mult. (1550->142), div. (0->0), fcn. (1020->8), ass. (0->52)
t55 = qJD(1) ^ 2;
t70 = t55 / 0.2e1;
t69 = cos(qJ(3));
t68 = cos(qJ(6));
t52 = sin(qJ(3));
t53 = sin(qJ(2));
t63 = qJD(1) * t53;
t33 = -t69 * qJD(2) + t52 * t63;
t35 = t52 * qJD(2) + t69 * t63;
t67 = t35 * t33;
t54 = cos(qJ(2));
t62 = t54 * qJD(1);
t43 = -qJD(3) + t62;
t66 = t43 * t33;
t65 = t54 * t55;
t29 = (-pkin(2) * t54 - pkin(8) * t53 - pkin(1)) * qJD(1);
t39 = pkin(7) * t62 + qJD(2) * pkin(8);
t24 = t69 * t29 - t52 * t39;
t57 = qJD(4) - t24;
t13 = -t35 * qJ(5) + (pkin(3) + pkin(4)) * t43 + t57;
t25 = t52 * t29 + t69 * t39;
t18 = -t43 * qJ(4) + t25;
t15 = t33 * qJ(5) + t18;
t49 = sin(pkin(10));
t64 = cos(pkin(10));
t6 = t49 * t13 + t64 * t15;
t38 = -qJD(2) * pkin(2) + pkin(7) * t63;
t61 = t33 ^ 2 / 0.2e1;
t40 = t43 ^ 2 / 0.2e1;
t60 = qJD(1) * qJD(2);
t59 = t53 * t60;
t58 = t54 * t60;
t5 = t64 * t13 - t49 * t15;
t19 = t33 * pkin(3) - t35 * qJ(4) + t38;
t16 = -t33 * pkin(4) + qJD(5) - t19;
t51 = sin(qJ(6));
t48 = t54 ^ 2;
t47 = t53 ^ 2;
t42 = qJD(6) + t43;
t30 = t35 ^ 2 / 0.2e1;
t26 = t35 * t43;
t23 = t49 * t33 + t64 * t35;
t21 = -t64 * t33 + t49 * t35;
t17 = t43 * pkin(3) + t57;
t10 = -t51 * t21 + t68 * t23;
t8 = t68 * t21 + t51 * t23;
t7 = t21 * pkin(5) + t16;
t4 = -t21 * pkin(9) + t6;
t3 = t43 * pkin(5) - t23 * pkin(9) + t5;
t2 = t51 * t3 + t68 * t4;
t1 = t68 * t3 - t51 * t4;
t9 = [0, 0, 0, 0, 0, t70, 0, 0, 0, 0, t47 * t70, t53 * t65, t59, t48 * t70, t58, qJD(2) ^ 2 / 0.2e1, pkin(1) * t65 - pkin(7) * t59, -t55 * pkin(1) * t53 - pkin(7) * t58 (t47 + t48) * t55 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t48 / 0.2e1 + t47 / 0.2e1) * pkin(7) ^ 2) * t55, t30, -t67, -t26, t61, t66, t40, -t24 * t43 + t38 * t33, t25 * t43 + t38 * t35, -t24 * t35 - t25 * t33, t25 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1, t30, -t26, t67, t40, -t66, t61, t17 * t43 + t19 * t33, t17 * t35 - t18 * t33, -t18 * t43 - t19 * t35, t18 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t23 ^ 2 / 0.2e1, -t23 * t21, t23 * t43, t21 ^ 2 / 0.2e1, -t21 * t43, t40, t16 * t21 + t5 * t43, t16 * t23 - t6 * t43, -t6 * t21 - t5 * t23, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, t10 * t42, t8 ^ 2 / 0.2e1, -t8 * t42, t42 ^ 2 / 0.2e1, t1 * t42 + t7 * t8, t7 * t10 - t2 * t42, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t9;
