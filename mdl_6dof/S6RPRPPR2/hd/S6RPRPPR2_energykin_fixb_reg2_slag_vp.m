% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPPR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR2_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:42:38
% EndTime: 2019-03-09 02:42:38
% DurationCPUTime: 0.17s
% Computational Cost: add. (394->58), mult. (934->128), div. (0->0), fcn. (573->8), ass. (0->49)
t36 = sin(pkin(10));
t38 = cos(pkin(10));
t41 = sin(qJ(3));
t42 = cos(qJ(3));
t26 = (t36 * t42 + t38 * t41) * qJD(1);
t43 = qJD(1) ^ 2;
t35 = t43 / 0.2e1;
t60 = pkin(4) + pkin(8);
t59 = pkin(1) * t43;
t58 = cos(qJ(6));
t55 = qJD(1) * t42;
t56 = qJD(1) * t41;
t24 = t36 * t56 - t38 * t55;
t57 = t26 * t24;
t37 = sin(pkin(9));
t29 = (pkin(1) * t37 + pkin(7)) * qJD(1);
t33 = t42 * qJD(2);
t50 = qJ(4) * qJD(1);
t14 = qJD(3) * pkin(3) + t33 + (-t29 - t50) * t41;
t21 = t41 * qJD(2) + t42 * t29;
t18 = t42 * t50 + t21;
t9 = t36 * t14 + t38 * t18;
t54 = qJD(3) * t24;
t53 = t26 * qJD(3);
t52 = t24 ^ 2 / 0.2e1;
t51 = t26 ^ 2 / 0.2e1;
t49 = qJD(1) * qJD(3);
t39 = cos(pkin(9));
t48 = -pkin(1) * t39 - pkin(2);
t8 = t38 * t14 - t36 * t18;
t47 = qJD(5) - t8;
t7 = -qJD(3) * qJ(5) - t9;
t23 = qJD(4) + (-pkin(3) * t42 + t48) * qJD(1);
t45 = -t26 * qJ(5) + t23;
t40 = sin(qJ(6));
t34 = qJD(3) ^ 2 / 0.2e1;
t30 = t48 * qJD(1);
t22 = qJD(6) + t26;
t20 = -t41 * t29 + t33;
t17 = t58 * qJD(3) + t40 * t24;
t15 = t40 * qJD(3) - t58 * t24;
t10 = t24 * pkin(4) + t45;
t6 = -qJD(3) * pkin(4) + t47;
t5 = t60 * t24 + t45;
t4 = -t24 * pkin(5) - t7;
t3 = t26 * pkin(5) - t60 * qJD(3) + t47;
t2 = t40 * t3 + t58 * t5;
t1 = t58 * t3 - t40 * t5;
t11 = [0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t39 * t59, -t37 * t59, 0, qJD(2) ^ 2 / 0.2e1 + (t37 ^ 2 / 0.2e1 + t39 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t43, t41 ^ 2 * t35, t41 * t43 * t42, t41 * t49, t42 ^ 2 * t35, t42 * t49, t34, t20 * qJD(3) - t30 * t55, -t21 * qJD(3) + t30 * t56 (-t20 * t41 + t21 * t42) * qJD(1), t21 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1, t51, -t57, t53, t52, -t54, t34, t8 * qJD(3) + t23 * t24, -t9 * qJD(3) + t23 * t26, -t9 * t24 - t8 * t26, t9 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1, t34, -t53, t54, t51, -t57, t52, t7 * t24 + t6 * t26, t6 * qJD(3) - t10 * t24, -t7 * qJD(3) - t10 * t26, t10 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t17 ^ 2 / 0.2e1, -t17 * t15, t17 * t22, t15 ^ 2 / 0.2e1, -t15 * t22, t22 ^ 2 / 0.2e1, t1 * t22 + t4 * t15, t4 * t17 - t2 * t22, -t1 * t17 - t2 * t15, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1;];
T_reg  = t11;
