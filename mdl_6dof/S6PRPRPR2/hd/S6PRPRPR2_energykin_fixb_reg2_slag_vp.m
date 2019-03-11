% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRPR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR2_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:32:52
% EndTime: 2019-03-08 19:32:52
% DurationCPUTime: 0.14s
% Computational Cost: add. (416->52), mult. (975->129), div. (0->0), fcn. (695->12), ass. (0->46)
t51 = qJD(2) ^ 2;
t40 = t51 / 0.2e1;
t50 = cos(qJ(2));
t43 = sin(pkin(6));
t58 = qJD(1) * t43;
t30 = qJD(2) * pkin(2) + t50 * t58;
t42 = sin(pkin(11));
t44 = cos(pkin(11));
t48 = sin(qJ(2));
t54 = t48 * t58;
t24 = t42 * t30 + t44 * t54;
t22 = qJD(2) * pkin(8) + t24;
t45 = cos(pkin(6));
t34 = t45 * qJD(1) + qJD(3);
t47 = sin(qJ(4));
t49 = cos(qJ(4));
t14 = t49 * t22 + t47 * t34;
t10 = qJD(4) * qJ(5) + t14;
t23 = t44 * t30 - t42 * t54;
t15 = (-pkin(4) * t49 - qJ(5) * t47 - pkin(3)) * qJD(2) - t23;
t41 = sin(pkin(12));
t59 = cos(pkin(12));
t6 = t59 * t10 + t41 * t15;
t60 = cos(qJ(6));
t57 = qJD(2) * t47;
t56 = t49 * qJD(2);
t55 = qJD(2) * qJD(4);
t53 = qJD(2) * t58;
t5 = -t41 * t10 + t59 * t15;
t13 = -t47 * t22 + t49 * t34;
t9 = -qJD(4) * pkin(4) + qJD(5) - t13;
t52 = qJD(1) ^ 2;
t46 = sin(qJ(6));
t37 = t49 ^ 2 * t40;
t35 = -qJD(6) + t56;
t29 = t41 * qJD(4) + t59 * t57;
t27 = -t59 * qJD(4) + t41 * t57;
t21 = -qJD(2) * pkin(3) - t23;
t18 = -t46 * t27 + t60 * t29;
t16 = t60 * t27 + t46 * t29;
t7 = t27 * pkin(5) + t9;
t4 = -t27 * pkin(9) + t6;
t3 = -pkin(5) * t56 - t29 * pkin(9) + t5;
t2 = t46 * t3 + t60 * t4;
t1 = t60 * t3 - t46 * t4;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t52 / 0.2e1, 0, 0, 0, 0, 0, t40, t50 * t53, -t48 * t53, 0 (t45 ^ 2 / 0.2e1 + (t48 ^ 2 / 0.2e1 + t50 ^ 2 / 0.2e1) * t43 ^ 2) * t52, 0, 0, 0, 0, 0, t40, t23 * qJD(2), -t24 * qJD(2), 0, t24 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1, t47 ^ 2 * t40, t47 * t51 * t49, t47 * t55, t37, t49 * t55, qJD(4) ^ 2 / 0.2e1, t13 * qJD(4) - t21 * t56, -t14 * qJD(4) + t21 * t57 (-t13 * t47 + t14 * t49) * qJD(2), t14 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t29 ^ 2 / 0.2e1, -t29 * t27, -t29 * t56, t27 ^ 2 / 0.2e1, t27 * t56, t37, t9 * t27 - t5 * t56, t9 * t29 + t6 * t56, -t6 * t27 - t5 * t29, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, -t18 * t35, t16 ^ 2 / 0.2e1, t16 * t35, t35 ^ 2 / 0.2e1, -t1 * t35 + t7 * t16, t7 * t18 + t2 * t35, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
