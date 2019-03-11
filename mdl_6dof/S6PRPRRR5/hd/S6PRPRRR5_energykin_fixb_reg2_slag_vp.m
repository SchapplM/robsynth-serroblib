% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRRR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR5_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:43:34
% EndTime: 2019-03-08 20:43:34
% DurationCPUTime: 0.18s
% Computational Cost: add. (357->51), mult. (737->123), div. (0->0), fcn. (481->10), ass. (0->45)
t38 = sin(qJ(5));
t39 = sin(qJ(4));
t41 = cos(qJ(5));
t42 = cos(qJ(4));
t22 = (t38 * t42 + t39 * t41) * qJD(2);
t44 = qJD(2) ^ 2;
t34 = t44 / 0.2e1;
t45 = qJD(1) ^ 2;
t56 = t45 / 0.2e1;
t43 = cos(qJ(2));
t35 = sin(pkin(6));
t54 = qJD(1) * t35;
t47 = -t43 * t54 + qJD(3);
t21 = (-pkin(2) - pkin(8)) * qJD(2) + t47;
t36 = cos(pkin(6));
t53 = qJD(1) * t36;
t12 = t42 * t21 - t39 * t53;
t10 = -t42 * qJD(2) * pkin(9) + qJD(4) * pkin(4) + t12;
t13 = t39 * t21 + t42 * t53;
t52 = qJD(2) * t39;
t11 = -pkin(9) * t52 + t13;
t6 = t38 * t10 + t41 * t11;
t55 = cos(qJ(6));
t40 = sin(qJ(2));
t26 = qJD(2) * qJ(3) + t40 * t54;
t51 = t26 * qJD(2);
t50 = t26 ^ 2 / 0.2e1;
t49 = qJD(2) * qJD(4);
t48 = qJD(2) * t54;
t5 = t41 * t10 - t38 * t11;
t19 = pkin(4) * t52 + t26;
t37 = sin(qJ(6));
t33 = qJD(4) + qJD(5);
t30 = t36 ^ 2 * t56;
t25 = -qJD(2) * pkin(2) + t47;
t24 = (-t38 * t39 + t41 * t42) * qJD(2);
t20 = qJD(6) + t22;
t16 = t55 * t24 + t37 * t33;
t14 = t37 * t24 - t55 * t33;
t7 = t22 * pkin(5) - t24 * pkin(10) + t19;
t4 = t33 * pkin(10) + t6;
t3 = -t33 * pkin(5) - t5;
t2 = t37 * t7 + t55 * t4;
t1 = -t37 * t4 + t55 * t7;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t56, 0, 0, 0, 0, 0, t34, t43 * t48, -t40 * t48, 0, t30 + (t40 ^ 2 / 0.2e1 + t43 ^ 2 / 0.2e1) * t45 * t35 ^ 2, t34, 0, 0, 0, 0, 0, 0, t25 * qJD(2), t51, t30 + t50 + t25 ^ 2 / 0.2e1, t42 ^ 2 * t34, -t42 * t44 * t39, t42 * t49, t39 ^ 2 * t34, -t39 * t49, qJD(4) ^ 2 / 0.2e1, t12 * qJD(4) + t39 * t51, -t13 * qJD(4) + t42 * t51 (-t12 * t42 - t13 * t39) * qJD(2), t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t50, t24 ^ 2 / 0.2e1, -t24 * t22, t24 * t33, t22 ^ 2 / 0.2e1, -t22 * t33, t33 ^ 2 / 0.2e1, t19 * t22 + t5 * t33, t19 * t24 - t6 * t33, -t6 * t22 - t5 * t24, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t16 ^ 2 / 0.2e1, -t16 * t14, t16 * t20, t14 ^ 2 / 0.2e1, -t14 * t20, t20 ^ 2 / 0.2e1, t1 * t20 + t3 * t14, t3 * t16 - t2 * t20, -t1 * t16 - t2 * t14, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
