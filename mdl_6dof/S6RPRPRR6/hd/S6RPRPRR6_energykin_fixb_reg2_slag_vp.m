% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRR6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR6_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:52:51
% EndTime: 2019-03-09 03:52:52
% DurationCPUTime: 0.21s
% Computational Cost: add. (1072->67), mult. (2673->152), div. (0->0), fcn. (2023->10), ass. (0->50)
t58 = qJD(1) ^ 2;
t67 = t58 / 0.2e1;
t56 = sin(qJ(3));
t57 = cos(qJ(3));
t53 = cos(pkin(10));
t61 = qJD(1) * t53;
t52 = sin(pkin(10));
t62 = qJD(1) * t52;
t40 = t56 * t62 - t57 * t61;
t42 = (t52 * t57 + t53 * t56) * qJD(1);
t45 = qJD(2) + (-pkin(2) * t53 - pkin(1)) * qJD(1);
t25 = t40 * pkin(3) - t42 * qJ(4) + t45;
t64 = pkin(7) + qJ(2);
t43 = t64 * t62;
t44 = t64 * t61;
t30 = -t56 * t43 + t57 * t44;
t28 = qJD(3) * qJ(4) + t30;
t51 = sin(pkin(11));
t63 = cos(pkin(11));
t17 = t51 * t25 + t63 * t28;
t32 = -t63 * qJD(3) + t51 * t42;
t15 = -t32 * pkin(8) + t17;
t55 = sin(qJ(5));
t66 = cos(qJ(5));
t16 = t63 * t25 - t51 * t28;
t34 = t51 * qJD(3) + t63 * t42;
t9 = t40 * pkin(4) - t34 * pkin(8) + t16;
t6 = t66 * t15 + t55 * t9;
t65 = cos(qJ(6));
t60 = t40 ^ 2 / 0.2e1;
t5 = -t55 * t15 + t66 * t9;
t29 = -t57 * t43 - t56 * t44;
t36 = qJD(5) + t40;
t27 = -qJD(3) * pkin(3) + qJD(4) - t29;
t18 = t32 * pkin(4) + t27;
t54 = sin(qJ(6));
t50 = t53 ^ 2;
t49 = t52 ^ 2;
t47 = -qJD(1) * pkin(1) + qJD(2);
t35 = qJD(6) + t36;
t22 = -t55 * t32 + t66 * t34;
t20 = t66 * t32 + t55 * t34;
t13 = -t54 * t20 + t65 * t22;
t11 = t65 * t20 + t54 * t22;
t10 = t20 * pkin(5) + t18;
t4 = -t20 * pkin(9) + t6;
t3 = t36 * pkin(5) - t22 * pkin(9) + t5;
t2 = t54 * t3 + t65 * t4;
t1 = t65 * t3 - t54 * t4;
t7 = [0, 0, 0, 0, 0, t67, 0, 0, 0, 0, t49 * t67, t52 * t58 * t53, 0, t50 * t67, 0, 0, -t47 * t61, t47 * t62 (t49 + t50) * t58 * qJ(2), t47 ^ 2 / 0.2e1 + (t50 / 0.2e1 + t49 / 0.2e1) * qJ(2) ^ 2 * t58, t42 ^ 2 / 0.2e1, -t42 * t40, t42 * qJD(3), t60, -t40 * qJD(3), qJD(3) ^ 2 / 0.2e1, t29 * qJD(3) + t45 * t40, -t30 * qJD(3) + t45 * t42, -t29 * t42 - t30 * t40, t30 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1 + t45 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, t34 * t40, t32 ^ 2 / 0.2e1, -t32 * t40, t60, t16 * t40 + t27 * t32, -t17 * t40 + t27 * t34, -t16 * t34 - t17 * t32, t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * t36, t20 ^ 2 / 0.2e1, -t20 * t36, t36 ^ 2 / 0.2e1, t18 * t20 + t5 * t36, t18 * t22 - t6 * t36, -t6 * t20 - t5 * t22, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t13 * t11, t13 * t35, t11 ^ 2 / 0.2e1, -t11 * t35, t35 ^ 2 / 0.2e1, t1 * t35 + t10 * t11, t10 * t13 - t2 * t35, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg  = t7;
