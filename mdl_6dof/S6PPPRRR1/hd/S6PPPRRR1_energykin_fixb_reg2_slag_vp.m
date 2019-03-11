% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PPPRRR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_energykin_fixb_reg2_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:40:33
% EndTime: 2019-03-08 18:40:33
% DurationCPUTime: 0.18s
% Computational Cost: add. (528->44), mult. (1358->116), div. (0->0), fcn. (1241->16), ass. (0->45)
t30 = sin(pkin(14));
t31 = sin(pkin(13));
t35 = cos(pkin(14));
t34 = sin(pkin(6));
t51 = qJD(1) * t34;
t36 = cos(pkin(13));
t38 = cos(pkin(7));
t52 = t36 * t38;
t25 = cos(pkin(6)) * qJD(1) + qJD(2);
t33 = sin(pkin(7));
t53 = t25 * t33;
t18 = t35 * t53 + (-t30 * t31 + t35 * t52) * t51;
t21 = -t33 * t36 * t51 + t38 * t25 + qJD(3);
t32 = sin(pkin(8));
t37 = cos(pkin(8));
t58 = t18 * t37 + t21 * t32;
t19 = t30 * t53 + (t30 * t52 + t31 * t35) * t51;
t41 = sin(qJ(4));
t43 = cos(qJ(4));
t11 = -t41 * t19 + t58 * t43;
t44 = qJD(4) ^ 2;
t57 = t44 / 0.2e1;
t12 = t43 * t19 + t58 * t41;
t10 = qJD(4) * pkin(10) + t12;
t14 = -t32 * t18 + t37 * t21;
t40 = sin(qJ(5));
t42 = cos(qJ(5));
t6 = t42 * t10 + t40 * t14;
t56 = cos(qJ(6));
t50 = qJD(4) * t40;
t49 = t42 * qJD(4);
t48 = qJD(4) * qJD(5);
t5 = -t40 * t10 + t42 * t14;
t45 = qJD(1) ^ 2;
t39 = sin(qJ(6));
t26 = -qJD(6) + t49;
t24 = t39 * qJD(5) + t56 * t50;
t22 = -t56 * qJD(5) + t39 * t50;
t9 = -qJD(4) * pkin(4) - t11;
t7 = (-pkin(5) * t42 - pkin(11) * t40 - pkin(4)) * qJD(4) - t11;
t4 = qJD(5) * pkin(11) + t6;
t3 = -qJD(5) * pkin(5) - t5;
t2 = t39 * t7 + t56 * t4;
t1 = -t39 * t4 + t56 * t7;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t45 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25 ^ 2 / 0.2e1 + (t31 ^ 2 / 0.2e1 + t36 ^ 2 / 0.2e1) * t45 * t34 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t57, t11 * qJD(4), -t12 * qJD(4), 0, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t40 ^ 2 * t57, t40 * t44 * t42, t40 * t48, t42 ^ 2 * t57, t42 * t48, qJD(5) ^ 2 / 0.2e1, t5 * qJD(5) - t9 * t49, -t6 * qJD(5) + t9 * t50 (-t40 * t5 + t42 * t6) * qJD(4), t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t24 ^ 2 / 0.2e1, -t22 * t24, -t26 * t24, t22 ^ 2 / 0.2e1, t22 * t26, t26 ^ 2 / 0.2e1, -t1 * t26 + t3 * t22, t2 * t26 + t3 * t24, -t1 * t24 - t2 * t22, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
