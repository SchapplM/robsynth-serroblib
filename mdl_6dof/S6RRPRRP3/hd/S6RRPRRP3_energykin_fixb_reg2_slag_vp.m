% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRP3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:50:33
% EndTime: 2019-03-09 11:50:33
% DurationCPUTime: 0.19s
% Computational Cost: add. (837->64), mult. (2016->140), div. (0->0), fcn. (1447->8), ass. (0->53)
t58 = qJD(1) ^ 2;
t68 = t58 / 0.2e1;
t51 = sin(pkin(10));
t52 = cos(pkin(10));
t57 = cos(qJ(2));
t63 = qJD(1) * t57;
t55 = sin(qJ(2));
t64 = qJD(1) * t55;
t38 = t51 * t64 - t52 * t63;
t40 = (t51 * t57 + t52 * t55) * qJD(1);
t45 = qJD(3) + (-pkin(2) * t57 - pkin(1)) * qJD(1);
t24 = t38 * pkin(3) - t40 * pkin(8) + t45;
t65 = pkin(7) + qJ(3);
t43 = qJD(2) * pkin(2) - t65 * t64;
t44 = t65 * t63;
t29 = t51 * t43 + t52 * t44;
t27 = qJD(2) * pkin(8) + t29;
t54 = sin(qJ(4));
t56 = cos(qJ(4));
t13 = t54 * t24 + t56 * t27;
t31 = -t56 * qJD(2) + t54 * t40;
t10 = -t31 * pkin(9) + t13;
t53 = sin(qJ(5));
t67 = cos(qJ(5));
t12 = t56 * t24 - t54 * t27;
t33 = t54 * qJD(2) + t56 * t40;
t37 = qJD(4) + t38;
t7 = t37 * pkin(4) - t33 * pkin(9) + t12;
t4 = t67 * t10 + t53 * t7;
t66 = t57 * t58;
t62 = qJD(1) * qJD(2);
t3 = -t53 * t10 + t67 * t7;
t61 = t55 * t62;
t60 = t57 * t62;
t28 = t52 * t43 - t51 * t44;
t26 = -qJD(2) * pkin(3) - t28;
t16 = t31 * pkin(4) + t26;
t50 = t57 ^ 2;
t49 = t55 ^ 2;
t48 = qJD(2) ^ 2 / 0.2e1;
t35 = qJD(5) + t37;
t34 = t35 ^ 2 / 0.2e1;
t21 = -t53 * t31 + t67 * t33;
t19 = t67 * t31 + t53 * t33;
t18 = t21 ^ 2 / 0.2e1;
t17 = t19 ^ 2 / 0.2e1;
t15 = t21 * t35;
t14 = t19 * t35;
t11 = t21 * t19;
t8 = t19 * pkin(5) + qJD(6) + t16;
t2 = -t19 * qJ(6) + t4;
t1 = t35 * pkin(5) - t21 * qJ(6) + t3;
t5 = [0, 0, 0, 0, 0, t68, 0, 0, 0, 0, t49 * t68, t55 * t66, t61, t50 * t68, t60, t48, pkin(1) * t66 - pkin(7) * t61, -t58 * pkin(1) * t55 - pkin(7) * t60 (t49 + t50) * t58 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t50 / 0.2e1 + t49 / 0.2e1) * pkin(7) ^ 2) * t58, t40 ^ 2 / 0.2e1, -t40 * t38, t40 * qJD(2), t38 ^ 2 / 0.2e1, -t38 * qJD(2), t48, t28 * qJD(2) + t45 * t38, -t29 * qJD(2) + t45 * t40, -t28 * t40 - t29 * t38, t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t45 ^ 2 / 0.2e1, t33 ^ 2 / 0.2e1, -t33 * t31, t33 * t37, t31 ^ 2 / 0.2e1, -t31 * t37, t37 ^ 2 / 0.2e1, t12 * t37 + t26 * t31, -t13 * t37 + t26 * t33, -t12 * t33 - t13 * t31, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1, t18, -t11, t15, t17, -t14, t34, t16 * t19 + t3 * t35, t16 * t21 - t4 * t35, -t4 * t19 - t3 * t21, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1, t18, -t11, t15, t17, -t14, t34, t1 * t35 + t8 * t19, -t2 * t35 + t8 * t21, -t1 * t21 - t2 * t19, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1;];
T_reg  = t5;
