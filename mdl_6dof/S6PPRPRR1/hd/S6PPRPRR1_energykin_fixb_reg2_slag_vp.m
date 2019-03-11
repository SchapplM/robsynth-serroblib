% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PPRPRR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_energykin_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:43:37
% EndTime: 2019-03-08 18:43:37
% DurationCPUTime: 0.14s
% Computational Cost: add. (386->44), mult. (994->115), div. (0->0), fcn. (845->14), ass. (0->43)
t46 = qJD(3) ^ 2;
t33 = t46 / 0.2e1;
t35 = sin(pkin(12));
t40 = cos(pkin(7));
t43 = sin(qJ(3));
t45 = cos(qJ(3));
t39 = cos(pkin(12));
t37 = sin(pkin(6));
t52 = qJD(1) * t37;
t48 = t39 * t52;
t28 = cos(pkin(6)) * qJD(1) + qJD(2);
t36 = sin(pkin(7));
t53 = t28 * t36;
t17 = -t43 * t35 * t52 + (t40 * t48 + t53) * t45;
t16 = qJD(3) * pkin(3) + t17;
t18 = t43 * t53 + (t39 * t40 * t43 + t35 * t45) * t52;
t34 = sin(pkin(13));
t38 = cos(pkin(13));
t12 = t34 * t16 + t38 * t18;
t10 = qJD(3) * pkin(9) + t12;
t21 = t40 * t28 - t36 * t48;
t20 = qJD(4) + t21;
t42 = sin(qJ(5));
t44 = cos(qJ(5));
t6 = t44 * t10 + t42 * t20;
t54 = cos(qJ(6));
t51 = qJD(3) * t42;
t50 = t44 * qJD(3);
t49 = qJD(3) * qJD(5);
t11 = t38 * t16 - t34 * t18;
t5 = -t42 * t10 + t44 * t20;
t47 = qJD(1) ^ 2;
t41 = sin(qJ(6));
t29 = -qJD(6) + t50;
t26 = t41 * qJD(5) + t54 * t51;
t24 = -t54 * qJD(5) + t41 * t51;
t9 = -qJD(3) * pkin(4) - t11;
t7 = (-pkin(5) * t44 - pkin(10) * t42 - pkin(4)) * qJD(3) - t11;
t4 = qJD(5) * pkin(10) + t6;
t3 = -qJD(5) * pkin(5) - t5;
t2 = t54 * t4 + t41 * t7;
t1 = -t41 * t4 + t54 * t7;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t47 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 ^ 2 / 0.2e1 + (t35 ^ 2 / 0.2e1 + t39 ^ 2 / 0.2e1) * t47 * t37 ^ 2, 0, 0, 0, 0, 0, t33, t17 * qJD(3), -t18 * qJD(3), 0, t18 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t33, t11 * qJD(3), -t12 * qJD(3), 0, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1, t42 ^ 2 * t33, t42 * t46 * t44, t42 * t49, t44 ^ 2 * t33, t44 * t49, qJD(5) ^ 2 / 0.2e1, t5 * qJD(5) - t9 * t50, -t6 * qJD(5) + t9 * t51 (-t42 * t5 + t44 * t6) * qJD(3), t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t26 ^ 2 / 0.2e1, -t26 * t24, -t26 * t29, t24 ^ 2 / 0.2e1, t24 * t29, t29 ^ 2 / 0.2e1, -t1 * t29 + t3 * t24, t2 * t29 + t3 * t26, -t1 * t26 - t2 * t24, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
