% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPP4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP4_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:02:33
% EndTime: 2019-03-09 21:02:34
% DurationCPUTime: 0.21s
% Computational Cost: add. (987->65), mult. (2165->142), div. (0->0), fcn. (1533->8), ass. (0->51)
t54 = qJD(1) ^ 2;
t68 = t54 / 0.2e1;
t49 = sin(pkin(10));
t62 = cos(pkin(10));
t52 = sin(qJ(2));
t53 = cos(qJ(2));
t33 = (-pkin(2) * t53 - pkin(8) * t52 - pkin(1)) * qJD(1);
t60 = t53 * qJD(1);
t40 = pkin(7) * t60 + qJD(2) * pkin(8);
t51 = sin(qJ(3));
t67 = cos(qJ(3));
t27 = t67 * t33 - t51 * t40;
t61 = qJD(1) * t52;
t36 = t51 * qJD(2) + t67 * t61;
t43 = -qJD(3) + t60;
t20 = -t43 * pkin(3) - t36 * pkin(9) + t27;
t28 = t51 * t33 + t67 * t40;
t34 = -t67 * qJD(2) + t51 * t61;
t22 = -t34 * pkin(9) + t28;
t50 = sin(qJ(4));
t66 = cos(qJ(4));
t10 = t66 * t20 - t50 * t22;
t26 = -t50 * t34 + t66 * t36;
t41 = -qJD(4) + t43;
t7 = -t41 * pkin(4) - t26 * qJ(5) + t10;
t11 = t50 * t20 + t66 * t22;
t24 = t66 * t34 + t50 * t36;
t9 = -t24 * qJ(5) + t11;
t4 = t49 * t7 + t62 * t9;
t14 = t62 * t24 + t49 * t26;
t16 = -t49 * t24 + t62 * t26;
t65 = t16 * t14;
t64 = t41 * t14;
t63 = t53 * t54;
t59 = t14 ^ 2 / 0.2e1;
t58 = qJD(1) * qJD(2);
t57 = t52 * t58;
t56 = t53 * t58;
t39 = -qJD(2) * pkin(2) + pkin(7) * t61;
t29 = t34 * pkin(3) + t39;
t3 = -t49 * t9 + t62 * t7;
t17 = t24 * pkin(4) + qJD(5) + t29;
t48 = t53 ^ 2;
t47 = t52 ^ 2;
t38 = t41 ^ 2 / 0.2e1;
t13 = t16 ^ 2 / 0.2e1;
t12 = t16 * t41;
t5 = t14 * pkin(5) - t16 * qJ(6) + t17;
t2 = -t41 * qJ(6) + t4;
t1 = t41 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t68, 0, 0, 0, 0, t47 * t68, t52 * t63, t57, t48 * t68, t56, qJD(2) ^ 2 / 0.2e1, pkin(1) * t63 - pkin(7) * t57, -t54 * pkin(1) * t52 - pkin(7) * t56 (t47 + t48) * t54 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t48 / 0.2e1 + t47 / 0.2e1) * pkin(7) ^ 2) * t54, t36 ^ 2 / 0.2e1, -t36 * t34, -t36 * t43, t34 ^ 2 / 0.2e1, t34 * t43, t43 ^ 2 / 0.2e1, -t27 * t43 + t39 * t34, t28 * t43 + t39 * t36, -t27 * t36 - t28 * t34, t28 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1 + t39 ^ 2 / 0.2e1, t26 ^ 2 / 0.2e1, -t26 * t24, -t26 * t41, t24 ^ 2 / 0.2e1, t24 * t41, t38, -t10 * t41 + t29 * t24, t11 * t41 + t29 * t26, -t10 * t26 - t11 * t24, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1, t13, -t65, -t12, t59, t64, t38, t17 * t14 - t3 * t41, t17 * t16 + t4 * t41, -t4 * t14 - t3 * t16, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t13, -t12, t65, t38, -t64, t59, t1 * t41 + t5 * t14, t1 * t16 - t2 * t14, -t5 * t16 - t2 * t41, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
