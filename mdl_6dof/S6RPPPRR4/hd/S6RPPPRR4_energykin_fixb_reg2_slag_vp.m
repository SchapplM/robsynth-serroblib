% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPPRR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR4_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:35:51
% EndTime: 2019-03-09 01:35:51
% DurationCPUTime: 0.13s
% Computational Cost: add. (235->45), mult. (402->93), div. (0->0), fcn. (163->6), ass. (0->35)
t34 = qJD(1) ^ 2;
t26 = t34 / 0.2e1;
t31 = sin(qJ(5));
t33 = cos(qJ(5));
t18 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t28 = sin(pkin(9));
t29 = cos(pkin(9));
t36 = qJ(2) * qJD(1);
t13 = t29 * t18 - t28 * t36;
t10 = qJD(1) * pkin(3) + qJD(4) - t13;
t9 = qJD(1) * pkin(7) + t10;
t6 = t33 * qJD(3) + t31 * t9;
t42 = sin(qJ(6));
t41 = t28 * t18;
t40 = qJD(1) * t33;
t14 = t29 * t36 + t41;
t27 = qJD(1) * qJ(4);
t11 = t14 - t27;
t39 = t11 * qJD(1);
t38 = t31 * qJD(1);
t37 = t11 ^ 2 / 0.2e1;
t35 = qJD(1) * qJD(5);
t5 = -t31 * qJD(3) + t33 * t9;
t32 = cos(qJ(6));
t25 = qJD(3) ^ 2 / 0.2e1;
t22 = -pkin(1) * qJD(1) + qJD(2);
t19 = -qJD(6) + t38;
t16 = -t42 * qJD(5) + t32 * t40;
t15 = t32 * qJD(5) + t42 * t40;
t7 = t41 - t27 + (-pkin(5) * t31 + pkin(8) * t33 + qJ(2) * t29) * qJD(1);
t4 = qJD(5) * pkin(8) + t6;
t3 = -qJD(5) * pkin(5) - t5;
t2 = t32 * t4 + t42 * t7;
t1 = t32 * t7 - t42 * t4;
t8 = [0, 0, 0, 0, 0, t26, 0, 0, 0, 0, 0, 0, 0, t26, 0, 0, -t22 * qJD(1), 0, t34 * qJ(2), qJ(2) ^ 2 * t26 + t22 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t26, -t13 * qJD(1), t14 * qJD(1), 0, t14 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1 + t25, t26, 0, 0, 0, 0, 0, 0, -t10 * qJD(1), -t39, t25 + t37 + t10 ^ 2 / 0.2e1, t33 ^ 2 * t26, -t33 * t34 * t31, -t33 * t35, t31 ^ 2 * t26, t31 * t35, qJD(5) ^ 2 / 0.2e1, t5 * qJD(5) - t11 * t38, -t6 * qJD(5) - t33 * t39 (t31 * t6 + t33 * t5) * qJD(1), t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t37, t16 ^ 2 / 0.2e1, -t16 * t15, t16 * t19, t15 ^ 2 / 0.2e1, -t15 * t19, t19 ^ 2 / 0.2e1, -t1 * t19 - t3 * t15, -t3 * t16 + t2 * t19, t1 * t16 + t2 * t15, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
