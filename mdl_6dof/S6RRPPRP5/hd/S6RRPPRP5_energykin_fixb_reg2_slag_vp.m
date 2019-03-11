% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPRP5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP5_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:43:56
% EndTime: 2019-03-09 08:43:56
% DurationCPUTime: 0.16s
% Computational Cost: add. (512->60), mult. (1119->126), div. (0->0), fcn. (649->6), ass. (0->51)
t48 = qJD(1) ^ 2;
t63 = t48 / 0.2e1;
t45 = sin(qJ(5));
t62 = cos(qJ(5));
t47 = cos(qJ(2));
t46 = sin(qJ(2));
t50 = -qJ(3) * t46 - pkin(1);
t58 = -pkin(2) - qJ(4);
t21 = (t58 * t47 + t50) * qJD(1);
t56 = t46 * qJD(1);
t54 = pkin(7) * t56 + qJD(3);
t22 = pkin(3) * t56 + t58 * qJD(2) + t54;
t43 = sin(pkin(9));
t44 = cos(pkin(9));
t10 = -t43 * t21 + t44 * t22;
t57 = qJD(1) * t47;
t28 = t44 * qJD(2) - t43 * t57;
t7 = pkin(4) * t56 - t28 * pkin(8) + t10;
t11 = t44 * t21 + t43 * t22;
t26 = t43 * qJD(2) + t44 * t57;
t9 = -t26 * pkin(8) + t11;
t4 = t45 * t7 + t62 * t9;
t14 = t62 * t26 + t45 * t28;
t16 = -t45 * t26 + t62 * t28;
t61 = t16 * t14;
t33 = qJD(5) + t56;
t60 = t33 * t14;
t59 = t47 * t48;
t30 = -pkin(7) * t57 - qJD(2) * qJ(3);
t55 = t14 ^ 2 / 0.2e1;
t53 = qJD(1) * qJD(2);
t52 = t46 * t53;
t51 = t47 * t53;
t24 = pkin(3) * t57 + qJD(4) - t30;
t3 = -t45 * t9 + t62 * t7;
t17 = t26 * pkin(4) + t24;
t42 = t47 ^ 2;
t41 = t46 ^ 2;
t39 = qJD(2) ^ 2 / 0.2e1;
t35 = t42 * t63;
t34 = t41 * t63;
t32 = t46 * t59;
t31 = t33 ^ 2 / 0.2e1;
t29 = -qJD(2) * pkin(2) + t54;
t25 = (-pkin(2) * t47 + t50) * qJD(1);
t13 = t16 ^ 2 / 0.2e1;
t12 = t16 * t33;
t5 = t14 * pkin(5) - t16 * qJ(6) + t17;
t2 = t33 * qJ(6) + t4;
t1 = -t33 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t63, 0, 0, 0, 0, t34, t32, t52, t35, t51, t39, pkin(1) * t59 - pkin(7) * t52, -t48 * pkin(1) * t46 - pkin(7) * t51 (t41 + t42) * t48 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t42 / 0.2e1 + t41 / 0.2e1) * pkin(7) ^ 2) * t48, t39, -t52, -t51, t34, t32, t35 (t29 * t46 - t30 * t47) * qJD(1), t29 * qJD(2) + t25 * t57, -t30 * qJD(2) - t25 * t56, t25 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1, t28 ^ 2 / 0.2e1, -t28 * t26, t28 * t56, t26 ^ 2 / 0.2e1, -t26 * t56, t34, t10 * t56 + t24 * t26, -t11 * t56 + t24 * t28, -t10 * t28 - t11 * t26, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1, t13, -t61, t12, t55, -t60, t31, t17 * t14 + t3 * t33, t17 * t16 - t4 * t33, -t4 * t14 - t3 * t16, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t13, t12, t61, t31, t60, t55, -t1 * t33 + t5 * t14, t1 * t16 - t2 * t14, -t5 * t16 + t2 * t33, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
