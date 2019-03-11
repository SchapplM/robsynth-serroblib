% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPP5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP5_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:06:16
% EndTime: 2019-03-09 10:06:16
% DurationCPUTime: 0.15s
% Computational Cost: add. (330->57), mult. (718->111), div. (0->0), fcn. (363->4), ass. (0->48)
t44 = qJD(1) ^ 2;
t59 = t44 / 0.2e1;
t58 = -pkin(2) - pkin(8);
t57 = -pkin(4) - pkin(5);
t40 = sin(qJ(4));
t42 = cos(qJ(4));
t43 = cos(qJ(2));
t54 = qJD(1) * t43;
t20 = t40 * qJD(2) + t42 * t54;
t22 = t42 * qJD(2) - t40 * t54;
t9 = t22 * t20;
t41 = sin(qJ(2));
t53 = t41 * qJD(1);
t29 = qJD(4) + t53;
t14 = t22 * t29;
t56 = t29 * t20;
t55 = t43 * t44;
t48 = -qJ(3) * t41 - pkin(1);
t13 = (t58 * t43 + t48) * qJD(1);
t52 = pkin(7) * t53 + qJD(3);
t15 = pkin(3) * t53 + t58 * qJD(2) + t52;
t7 = t42 * t13 + t40 * t15;
t24 = -pkin(7) * t54 - qJD(2) * qJ(3);
t18 = t20 ^ 2 / 0.2e1;
t25 = t29 ^ 2 / 0.2e1;
t51 = qJD(1) * qJD(2);
t5 = t29 * qJ(5) + t7;
t16 = pkin(3) * t54 - t24;
t50 = t41 * t51;
t49 = t43 * t51;
t6 = -t40 * t13 + t42 * t15;
t47 = qJD(5) - t6;
t46 = t22 * qJ(5) - t16;
t39 = t43 ^ 2;
t38 = t41 ^ 2;
t36 = qJD(2) ^ 2 / 0.2e1;
t32 = t39 * t59;
t31 = t38 * t59;
t28 = t41 * t55;
t23 = -qJD(2) * pkin(2) + t52;
t19 = t22 ^ 2 / 0.2e1;
t17 = (-pkin(2) * t43 + t48) * qJD(1);
t8 = t20 * pkin(4) - t46;
t4 = -t29 * pkin(4) + t47;
t3 = t57 * t20 + qJD(6) + t46;
t2 = t20 * qJ(6) + t5;
t1 = -t22 * qJ(6) + t57 * t29 + t47;
t10 = [0, 0, 0, 0, 0, t59, 0, 0, 0, 0, t31, t28, t50, t32, t49, t36, pkin(1) * t55 - pkin(7) * t50, -t44 * pkin(1) * t41 - pkin(7) * t49 (t38 + t39) * t44 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t39 / 0.2e1 + t38 / 0.2e1) * pkin(7) ^ 2) * t44, t36, -t50, -t49, t31, t28, t32 (t23 * t41 - t24 * t43) * qJD(1), t23 * qJD(2) + t17 * t54, -t24 * qJD(2) - t17 * t53, t17 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1, t19, -t9, t14, t18, -t56, t25, t16 * t20 + t6 * t29, t16 * t22 - t7 * t29, -t7 * t20 - t6 * t22, t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1, t19, t14, t9, t25, t56, t18, t8 * t20 - t4 * t29, -t5 * t20 + t4 * t22, -t8 * t22 + t5 * t29, t5 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1, t19, t9, -t14, t18, -t56, t25, -t1 * t29 - t3 * t20, t2 * t29 + t3 * t22, -t1 * t22 + t2 * t20, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t10;
