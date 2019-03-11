% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPRP4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP4_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:43:57
% EndTime: 2019-03-08 21:43:57
% DurationCPUTime: 0.15s
% Computational Cost: add. (283->55), mult. (660->113), div. (0->0), fcn. (399->8), ass. (0->52)
t48 = qJD(2) ^ 2;
t62 = t48 / 0.2e1;
t61 = -pkin(3) - pkin(9);
t46 = cos(qJ(3));
t43 = sin(qJ(3));
t51 = -qJ(4) * t43 - pkin(2);
t47 = cos(qJ(2));
t40 = sin(pkin(6));
t60 = qJD(1) * t40;
t55 = t47 * t60;
t11 = -t55 + (t61 * t46 + t51) * qJD(2);
t42 = sin(qJ(5));
t45 = cos(qJ(5));
t44 = sin(qJ(2));
t27 = qJD(2) * pkin(8) + t44 * t60;
t41 = cos(pkin(6));
t59 = qJD(1) * t41;
t15 = -t43 * t27 + t46 * t59;
t50 = qJD(4) - t15;
t57 = t43 * qJD(2);
t8 = pkin(4) * t57 + t61 * qJD(3) + t50;
t4 = t45 * t11 + t42 * t8;
t16 = t46 * t27 + t43 * t59;
t58 = qJD(2) * t46;
t56 = qJD(2) * qJD(3);
t13 = -qJD(3) * qJ(4) - t16;
t3 = -t42 * t11 + t45 * t8;
t54 = qJD(2) * t60;
t53 = t43 * t56;
t52 = t46 * t56;
t9 = pkin(4) * t58 - t13;
t49 = qJD(1) ^ 2;
t38 = qJD(3) ^ 2 / 0.2e1;
t35 = t46 ^ 2 * t62;
t34 = t43 ^ 2 * t62;
t33 = qJD(5) + t57;
t32 = t43 * t48 * t46;
t29 = t33 ^ 2 / 0.2e1;
t28 = -qJD(2) * pkin(2) - t55;
t26 = t45 * qJD(3) - t42 * t58;
t24 = t42 * qJD(3) + t45 * t58;
t21 = t26 ^ 2 / 0.2e1;
t20 = t24 ^ 2 / 0.2e1;
t19 = t26 * t33;
t18 = t24 * t33;
t17 = -t55 + (-pkin(3) * t46 + t51) * qJD(2);
t14 = t26 * t24;
t12 = -qJD(3) * pkin(3) + t50;
t5 = t24 * pkin(5) + qJD(6) + t9;
t2 = -t24 * qJ(6) + t4;
t1 = t33 * pkin(5) - t26 * qJ(6) + t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t49 / 0.2e1, 0, 0, 0, 0, 0, t62, t47 * t54, -t44 * t54, 0 (t41 ^ 2 / 0.2e1 + (t44 ^ 2 / 0.2e1 + t47 ^ 2 / 0.2e1) * t40 ^ 2) * t49, t34, t32, t53, t35, t52, t38, t15 * qJD(3) - t28 * t58, -t16 * qJD(3) + t28 * t57 (-t15 * t43 + t16 * t46) * qJD(2), t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1, t38, -t53, -t52, t34, t32, t35 (t12 * t43 - t13 * t46) * qJD(2), t12 * qJD(3) + t17 * t58, -t13 * qJD(3) - t17 * t57, t17 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1, t21, -t14, t19, t20, -t18, t29, t9 * t24 + t3 * t33, t9 * t26 - t4 * t33, -t4 * t24 - t3 * t26, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t21, -t14, t19, t20, -t18, t29, t1 * t33 + t5 * t24, -t2 * t33 + t5 * t26, -t1 * t26 - t2 * t24, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t6;
