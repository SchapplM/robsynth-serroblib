% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRPR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR2_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:42:36
% EndTime: 2019-03-09 01:42:36
% DurationCPUTime: 0.17s
% Computational Cost: add. (371->57), mult. (895->122), div. (0->0), fcn. (565->8), ass. (0->47)
t36 = sin(pkin(10));
t38 = cos(pkin(10));
t41 = sin(qJ(4));
t42 = cos(qJ(4));
t26 = (t36 * t42 + t38 * t41) * qJD(1);
t43 = qJD(1) ^ 2;
t35 = t43 / 0.2e1;
t58 = pkin(4) + pkin(8);
t57 = pkin(1) * t43;
t56 = cos(qJ(6));
t53 = qJD(1) * t38;
t54 = qJD(1) * t36;
t24 = t41 * t54 - t42 * t53;
t55 = t26 * t24;
t37 = sin(pkin(9));
t30 = (pkin(1) * t37 + qJ(3)) * qJD(1);
t33 = t38 * qJD(2);
t14 = t33 + (-pkin(7) * qJD(1) - t30) * t36;
t20 = t36 * qJD(2) + t38 * t30;
t15 = pkin(7) * t53 + t20;
t9 = t41 * t14 + t42 * t15;
t52 = qJD(4) * t24;
t51 = t26 * qJD(4);
t50 = t24 ^ 2 / 0.2e1;
t49 = t26 ^ 2 / 0.2e1;
t39 = cos(pkin(9));
t48 = -pkin(1) * t39 - pkin(2);
t8 = t42 * t14 - t41 * t15;
t47 = qJD(5) - t8;
t7 = -qJD(4) * qJ(5) - t9;
t23 = qJD(3) + (-pkin(3) * t38 + t48) * qJD(1);
t45 = -t26 * qJ(5) + t23;
t40 = sin(qJ(6));
t34 = qJD(4) ^ 2 / 0.2e1;
t29 = t48 * qJD(1) + qJD(3);
t22 = qJD(6) + t26;
t19 = -t36 * t30 + t33;
t18 = t56 * qJD(4) + t40 * t24;
t16 = t40 * qJD(4) - t56 * t24;
t10 = t24 * pkin(4) + t45;
t6 = -qJD(4) * pkin(4) + t47;
t5 = t58 * t24 + t45;
t4 = -t24 * pkin(5) - t7;
t3 = t26 * pkin(5) - t58 * qJD(4) + t47;
t2 = t40 * t3 + t56 * t5;
t1 = t56 * t3 - t40 * t5;
t11 = [0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t39 * t57, -t37 * t57, 0, qJD(2) ^ 2 / 0.2e1 + (t37 ^ 2 / 0.2e1 + t39 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t43, t36 ^ 2 * t35, t36 * t43 * t38, 0, t38 ^ 2 * t35, 0, 0, -t29 * t53, t29 * t54 (-t19 * t36 + t20 * t38) * qJD(1), t20 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1, t49, -t55, t51, t50, -t52, t34, t8 * qJD(4) + t23 * t24, -t9 * qJD(4) + t23 * t26, -t9 * t24 - t8 * t26, t9 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1, t34, -t51, t52, t49, -t55, t50, t7 * t24 + t6 * t26, t6 * qJD(4) - t10 * t24, -t7 * qJD(4) - t10 * t26, t10 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t22, t16 ^ 2 / 0.2e1, -t16 * t22, t22 ^ 2 / 0.2e1, t1 * t22 + t4 * t16, t4 * t18 - t2 * t22, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1;];
T_reg  = t11;
