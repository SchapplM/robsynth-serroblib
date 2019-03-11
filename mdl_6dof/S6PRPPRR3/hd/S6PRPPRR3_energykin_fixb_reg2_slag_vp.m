% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPPRR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:23:01
% EndTime: 2019-03-08 19:23:02
% DurationCPUTime: 0.12s
% Computational Cost: add. (264->44), mult. (516->106), div. (0->0), fcn. (301->10), ass. (0->40)
t39 = qJD(2) ^ 2;
t28 = t39 / 0.2e1;
t40 = qJD(1) ^ 2;
t48 = t40 / 0.2e1;
t38 = cos(qJ(2));
t30 = sin(pkin(6));
t46 = qJD(1) * t30;
t41 = -t38 * t46 + qJD(3);
t14 = (-pkin(2) - pkin(3)) * qJD(2) + t41;
t35 = sin(qJ(2));
t21 = qJD(2) * qJ(3) + t35 * t46;
t29 = sin(pkin(11));
t31 = cos(pkin(11));
t12 = t29 * t14 + t31 * t21;
t10 = -qJD(2) * pkin(8) + t12;
t32 = cos(pkin(6));
t23 = -t32 * qJD(1) + qJD(4);
t34 = sin(qJ(5));
t37 = cos(qJ(5));
t6 = t37 * t10 + t34 * t23;
t47 = sin(qJ(6));
t45 = qJD(2) * t34;
t44 = t37 * qJD(2);
t43 = qJD(2) * qJD(5);
t42 = qJD(2) * t46;
t11 = t31 * t14 - t29 * t21;
t9 = qJD(2) * pkin(4) - t11;
t5 = -t34 * t10 + t37 * t23;
t36 = cos(qJ(6));
t25 = t32 ^ 2 * t48;
t24 = qJD(6) + t44;
t19 = -t47 * qJD(5) + t36 * t45;
t18 = t36 * qJD(5) + t47 * t45;
t17 = -qJD(2) * pkin(2) + t41;
t7 = (pkin(5) * t37 + pkin(9) * t34) * qJD(2) + t9;
t4 = qJD(5) * pkin(9) + t6;
t3 = -qJD(5) * pkin(5) - t5;
t2 = t36 * t4 + t47 * t7;
t1 = t36 * t7 - t47 * t4;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, 0, 0, 0, 0, t28, t38 * t42, -t35 * t42, 0, t25 + (t35 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1) * t40 * t30 ^ 2, 0, 0, 0, t28, 0, 0, -t17 * qJD(2), 0, t21 * qJD(2), t21 ^ 2 / 0.2e1 + t25 + t17 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t28, -t11 * qJD(2), t12 * qJD(2), 0, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1, t34 ^ 2 * t28, t34 * t39 * t37, -t34 * t43, t37 ^ 2 * t28, -t37 * t43, qJD(5) ^ 2 / 0.2e1, t5 * qJD(5) + t9 * t44, -t6 * qJD(5) - t9 * t45 (t34 * t5 - t37 * t6) * qJD(2), t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t19 ^ 2 / 0.2e1, -t19 * t18, -t19 * t24, t18 ^ 2 / 0.2e1, t18 * t24, t24 ^ 2 / 0.2e1, t1 * t24 - t3 * t18, -t3 * t19 - t2 * t24, t1 * t19 + t2 * t18, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
