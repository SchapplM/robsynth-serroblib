% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PPRRRP1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP1_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:54:23
% EndTime: 2019-03-08 18:54:23
% DurationCPUTime: 0.18s
% Computational Cost: add. (374->48), mult. (946->112), div. (0->0), fcn. (766->12), ass. (0->47)
t35 = cos(pkin(6)) * qJD(1) + qJD(2);
t41 = sin(pkin(7));
t44 = cos(pkin(7));
t43 = cos(pkin(12));
t42 = sin(pkin(6));
t60 = qJD(1) * t42;
t55 = t43 * t60;
t64 = t35 * t41 + t44 * t55;
t47 = sin(qJ(3));
t49 = cos(qJ(3));
t40 = sin(pkin(12));
t56 = t40 * t60;
t18 = -t47 * t56 + t49 * t64;
t50 = qJD(3) ^ 2;
t63 = t50 / 0.2e1;
t46 = sin(qJ(4));
t48 = cos(qJ(4));
t13 = (-pkin(4) * t48 - pkin(10) * t46 - pkin(3)) * qJD(3) - t18;
t45 = sin(qJ(5));
t62 = cos(qJ(5));
t19 = t64 * t47 + t49 * t56;
t17 = qJD(3) * pkin(9) + t19;
t22 = t44 * t35 - t41 * t55;
t10 = t48 * t17 + t46 * t22;
t8 = qJD(4) * pkin(10) + t10;
t4 = t45 * t13 + t62 * t8;
t59 = qJD(3) * t46;
t58 = t48 * qJD(3);
t57 = qJD(3) * qJD(4);
t3 = t62 * t13 - t45 * t8;
t9 = -t46 * t17 + t48 * t22;
t7 = -qJD(4) * pkin(4) - t9;
t51 = qJD(1) ^ 2;
t36 = -qJD(5) + t58;
t34 = t36 ^ 2 / 0.2e1;
t31 = t45 * qJD(4) + t62 * t59;
t29 = -t62 * qJD(4) + t45 * t59;
t26 = t31 ^ 2 / 0.2e1;
t25 = t29 ^ 2 / 0.2e1;
t24 = t31 * t36;
t23 = t29 * t36;
t20 = t31 * t29;
t16 = -qJD(3) * pkin(3) - t18;
t5 = t29 * pkin(5) + qJD(6) + t7;
t2 = -t29 * qJ(6) + t4;
t1 = -t36 * pkin(5) - t31 * qJ(6) + t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t51 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35 ^ 2 / 0.2e1 + (t40 ^ 2 / 0.2e1 + t43 ^ 2 / 0.2e1) * t51 * t42 ^ 2, 0, 0, 0, 0, 0, t63, t18 * qJD(3), -t19 * qJD(3), 0, t19 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t46 ^ 2 * t63, t46 * t50 * t48, t46 * t57, t48 ^ 2 * t63, t48 * t57, qJD(4) ^ 2 / 0.2e1, t9 * qJD(4) - t16 * t58, -t10 * qJD(4) + t16 * t59 (t10 * t48 - t46 * t9) * qJD(3), t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1, t26, -t20, -t24, t25, t23, t34, t7 * t29 - t3 * t36, t7 * t31 + t4 * t36, -t4 * t29 - t3 * t31, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t26, -t20, -t24, t25, t23, t34, -t1 * t36 + t5 * t29, t2 * t36 + t5 * t31, -t1 * t31 - t2 * t29, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t6;
