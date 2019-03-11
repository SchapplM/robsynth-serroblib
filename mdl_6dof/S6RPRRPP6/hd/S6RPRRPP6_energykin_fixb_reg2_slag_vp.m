% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPP6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP6_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:48:26
% EndTime: 2019-03-09 04:48:27
% DurationCPUTime: 0.14s
% Computational Cost: add. (470->53), mult. (940->113), div. (0->0), fcn. (558->6), ass. (0->42)
t42 = qJD(1) ^ 2;
t35 = t42 / 0.2e1;
t38 = sin(pkin(9));
t47 = cos(pkin(9));
t40 = sin(qJ(3));
t41 = cos(qJ(3));
t22 = (pkin(3) * t40 - pkin(8) * t41 + qJ(2)) * qJD(1);
t30 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t23 = qJD(3) * pkin(8) + t40 * t30;
t39 = sin(qJ(4));
t51 = cos(qJ(4));
t11 = t51 * t22 - t39 * t23;
t46 = qJD(1) * t41;
t27 = t39 * qJD(3) + t51 * t46;
t31 = t40 * qJD(1) + qJD(4);
t7 = t31 * pkin(4) - t27 * qJ(5) + t11;
t12 = t39 * t22 + t51 * t23;
t25 = -t51 * qJD(3) + t39 * t46;
t9 = -t25 * qJ(5) + t12;
t4 = t38 * t7 + t47 * t9;
t14 = t47 * t25 + t38 * t27;
t16 = -t38 * t25 + t47 * t27;
t50 = t16 * t14;
t49 = t31 * t14;
t48 = t42 * qJ(2);
t45 = qJD(3) * t30;
t44 = t14 ^ 2 / 0.2e1;
t43 = qJD(1) * qJD(3);
t24 = -qJD(3) * pkin(3) - t41 * t30;
t3 = -t38 * t9 + t47 * t7;
t17 = t25 * pkin(4) + qJD(5) + t24;
t37 = t41 ^ 2;
t36 = t40 ^ 2;
t33 = qJ(2) ^ 2 * t35;
t32 = -pkin(1) * qJD(1) + qJD(2);
t28 = t31 ^ 2 / 0.2e1;
t13 = t16 ^ 2 / 0.2e1;
t10 = t16 * t31;
t5 = t14 * pkin(5) - t16 * qJ(6) + t17;
t2 = t31 * qJ(6) + t4;
t1 = -t31 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t35, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, t32 * qJD(1), t48, t33 + t32 ^ 2 / 0.2e1, t37 * t35, -t41 * t42 * t40, t41 * t43, t36 * t35, -t40 * t43, qJD(3) ^ 2 / 0.2e1, t40 * t48 + t41 * t45, -t40 * t45 + t41 * t48 (-t36 - t37) * t30 * qJD(1), t33 + (t36 / 0.2e1 + t37 / 0.2e1) * t30 ^ 2, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t31, t25 ^ 2 / 0.2e1, -t25 * t31, t28, t11 * t31 + t24 * t25, -t12 * t31 + t24 * t27, -t11 * t27 - t12 * t25, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1, t13, -t50, t10, t44, -t49, t28, t17 * t14 + t3 * t31, t17 * t16 - t4 * t31, -t4 * t14 - t3 * t16, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t13, t10, t50, t28, t49, t44, -t1 * t31 + t5 * t14, t1 * t16 - t2 * t14, -t5 * t16 + t2 * t31, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
