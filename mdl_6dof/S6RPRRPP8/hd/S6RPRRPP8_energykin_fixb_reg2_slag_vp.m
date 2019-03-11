% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPP8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP8_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:55:51
% EndTime: 2019-03-09 04:55:51
% DurationCPUTime: 0.15s
% Computational Cost: add. (298->50), mult. (588->98), div. (0->0), fcn. (296->4), ass. (0->39)
t36 = qJD(1) ^ 2;
t30 = t36 / 0.2e1;
t45 = cos(qJ(4));
t33 = sin(qJ(4));
t35 = cos(qJ(3));
t41 = qJD(1) * t35;
t19 = -t45 * qJD(3) + t33 * t41;
t21 = t33 * qJD(3) + t45 * t41;
t44 = t19 * t21;
t34 = sin(qJ(3));
t26 = t34 * qJD(1) + qJD(4);
t9 = t21 * t26;
t10 = t26 * t19;
t43 = pkin(4) + qJ(6);
t14 = (pkin(3) * t34 - pkin(8) * t35 + qJ(2)) * qJD(1);
t25 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t15 = qJD(3) * pkin(8) + t34 * t25;
t8 = t33 * t14 + t45 * t15;
t42 = t36 * qJ(2);
t40 = qJD(3) * t25;
t17 = t19 ^ 2 / 0.2e1;
t18 = t21 ^ 2 / 0.2e1;
t39 = qJD(1) * qJD(3);
t5 = -t26 * qJ(5) - t8;
t7 = t45 * t14 - t33 * t15;
t16 = -qJD(3) * pkin(3) - t35 * t25;
t38 = qJD(5) - t7;
t37 = -t21 * qJ(5) + t16;
t32 = t35 ^ 2;
t31 = t34 ^ 2;
t28 = qJ(2) ^ 2 * t30;
t27 = -pkin(1) * qJD(1) + qJD(2);
t23 = t26 ^ 2 / 0.2e1;
t6 = t19 * pkin(4) + t37;
t4 = -t26 * pkin(4) + t38;
t3 = t43 * t19 + t37;
t2 = -t19 * pkin(5) + qJD(6) - t5;
t1 = t21 * pkin(5) - t43 * t26 + t38;
t11 = [0, 0, 0, 0, 0, t30, 0, 0, 0, 0, t30, 0, 0, 0, 0, 0, 0, t27 * qJD(1), t42, t28 + t27 ^ 2 / 0.2e1, t32 * t30, -t35 * t36 * t34, t35 * t39, t31 * t30, -t34 * t39, qJD(3) ^ 2 / 0.2e1, t34 * t42 + t35 * t40, -t34 * t40 + t35 * t42 (-t31 - t32) * t25 * qJD(1), t28 + (t31 / 0.2e1 + t32 / 0.2e1) * t25 ^ 2, t18, -t44, t9, t17, -t10, t23, t16 * t19 + t7 * t26, t16 * t21 - t8 * t26, -t8 * t19 - t7 * t21, t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1, t23, -t9, t10, t18, -t44, t17, t5 * t19 + t4 * t21, -t6 * t19 + t4 * t26, -t6 * t21 - t5 * t26, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1, t23, t10, t9, t17, t44, t18, t1 * t21 - t2 * t19, t2 * t26 - t3 * t21, -t1 * t26 + t3 * t19, t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1;];
T_reg  = t11;
