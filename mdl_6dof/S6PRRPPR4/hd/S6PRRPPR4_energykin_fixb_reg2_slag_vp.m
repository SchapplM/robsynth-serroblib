% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPPR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR4_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:16:25
% EndTime: 2019-03-08 21:16:25
% DurationCPUTime: 0.16s
% Computational Cost: add. (394->58), mult. (939->130), div. (0->0), fcn. (636->10), ass. (0->50)
t46 = qJD(2) ^ 2;
t63 = t46 / 0.2e1;
t62 = cos(qJ(6));
t38 = sin(pkin(11));
t42 = sin(qJ(3));
t57 = qJD(2) * t42;
t60 = cos(pkin(11));
t26 = -t60 * qJD(3) + t38 * t57;
t28 = t38 * qJD(3) + t60 * t57;
t61 = t28 * t26;
t43 = sin(qJ(2));
t39 = sin(pkin(6));
t59 = qJD(1) * t39;
t29 = qJD(2) * pkin(8) + t43 * t59;
t44 = cos(qJ(3));
t40 = cos(pkin(6));
t58 = qJD(1) * t40;
t20 = t44 * t29 + t42 * t58;
t18 = qJD(3) * qJ(4) + t20;
t45 = cos(qJ(2));
t51 = t45 * t59;
t21 = -t51 + (-pkin(3) * t44 - qJ(4) * t42 - pkin(2)) * qJD(2);
t10 = t60 * t18 + t38 * t21;
t19 = -t42 * t29 + t44 * t58;
t56 = t44 * qJD(2);
t55 = t26 ^ 2 / 0.2e1;
t54 = qJD(2) * qJD(3);
t53 = t26 * t56;
t52 = t28 * t56;
t50 = qJD(2) * t59;
t9 = -t38 * t18 + t60 * t21;
t7 = -qJ(5) * t56 + t10;
t49 = qJD(3) * pkin(3) - qJD(4) + t19;
t6 = pkin(4) * t56 + qJD(5) - t9;
t48 = t28 * qJ(5) + t49;
t47 = qJD(1) ^ 2;
t41 = sin(qJ(6));
t34 = t44 ^ 2 * t63;
t33 = qJD(6) + t56;
t30 = -qJD(2) * pkin(2) - t51;
t23 = t28 ^ 2 / 0.2e1;
t13 = t41 * t26 + t62 * t28;
t11 = -t62 * t26 + t41 * t28;
t8 = t26 * pkin(4) - t48;
t5 = (-pkin(4) - pkin(5)) * t26 + t48;
t4 = t26 * pkin(9) + t7;
t3 = pkin(5) * t56 - t28 * pkin(9) + t6;
t2 = t41 * t3 + t62 * t4;
t1 = t62 * t3 - t41 * t4;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t47 / 0.2e1, 0, 0, 0, 0, 0, t63, t45 * t50, -t43 * t50, 0 (t40 ^ 2 / 0.2e1 + (t43 ^ 2 / 0.2e1 + t45 ^ 2 / 0.2e1) * t39 ^ 2) * t47, t42 ^ 2 * t63, t42 * t46 * t44, t42 * t54, t34, t44 * t54, qJD(3) ^ 2 / 0.2e1, t19 * qJD(3) - t30 * t56, -t20 * qJD(3) + t30 * t57 (-t19 * t42 + t20 * t44) * qJD(2), t20 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1, t23, -t61, -t52, t55, t53, t34, -t26 * t49 - t9 * t56, t10 * t56 - t28 * t49, -t10 * t26 - t9 * t28, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t49 ^ 2 / 0.2e1, t23, -t52, t61, t34, -t53, t55, t8 * t26 + t6 * t56, -t7 * t26 + t6 * t28, -t8 * t28 - t7 * t56, t7 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t13 * t11, t13 * t33, t11 ^ 2 / 0.2e1, -t11 * t33, t33 ^ 2 / 0.2e1, t1 * t33 + t5 * t11, t5 * t13 - t2 * t33, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t12;
