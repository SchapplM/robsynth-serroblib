% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPRP4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP4_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:39:40
% EndTime: 2019-03-09 08:39:40
% DurationCPUTime: 0.17s
% Computational Cost: add. (491->59), mult. (1143->123), div. (0->0), fcn. (710->6), ass. (0->50)
t47 = qJD(1) ^ 2;
t64 = t47 / 0.2e1;
t45 = sin(qJ(2));
t46 = cos(qJ(2));
t25 = (-pkin(2) * t46 - qJ(3) * t45 - pkin(1)) * qJD(1);
t56 = t46 * qJD(1);
t33 = pkin(7) * t56 + qJD(2) * qJ(3);
t42 = sin(pkin(9));
t58 = cos(pkin(9));
t20 = t42 * t25 + t58 * t33;
t14 = -qJ(4) * t56 + t20;
t57 = qJD(1) * t45;
t27 = -t58 * qJD(2) + t42 * t57;
t10 = t27 * pkin(8) + t14;
t44 = sin(qJ(5));
t63 = cos(qJ(5));
t19 = t58 * t25 - t42 * t33;
t13 = pkin(3) * t56 + qJD(4) - t19;
t29 = t42 * qJD(2) + t58 * t57;
t7 = pkin(4) * t56 - t29 * pkin(8) + t13;
t5 = t63 * t10 + t44 * t7;
t16 = -t63 * t27 + t44 * t29;
t18 = t44 * t27 + t63 * t29;
t62 = t18 * t16;
t61 = t29 * t27;
t35 = qJD(5) + t56;
t60 = t35 * t16;
t59 = t46 * t47;
t55 = t16 ^ 2 / 0.2e1;
t54 = t27 ^ 2 / 0.2e1;
t53 = qJD(1) * qJD(2);
t52 = t27 * t56;
t51 = t29 * t56;
t32 = -qJD(2) * pkin(2) + pkin(7) * t57 + qJD(3);
t50 = t45 * t53;
t49 = t46 * t53;
t12 = t27 * pkin(3) - t29 * qJ(4) + t32;
t4 = -t44 * t10 + t63 * t7;
t9 = -t27 * pkin(4) - t12;
t41 = t46 ^ 2;
t40 = t45 ^ 2;
t36 = t41 * t64;
t34 = t35 ^ 2 / 0.2e1;
t23 = t29 ^ 2 / 0.2e1;
t15 = t18 ^ 2 / 0.2e1;
t11 = t18 * t35;
t3 = t16 * pkin(5) - t18 * qJ(6) + t9;
t2 = t35 * qJ(6) + t5;
t1 = -t35 * pkin(5) + qJD(6) - t4;
t6 = [0, 0, 0, 0, 0, t64, 0, 0, 0, 0, t40 * t64, t45 * t59, t50, t36, t49, qJD(2) ^ 2 / 0.2e1, pkin(1) * t59 - pkin(7) * t50, -t47 * pkin(1) * t45 - pkin(7) * t49 (t40 + t41) * t47 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t41 / 0.2e1 + t40 / 0.2e1) * pkin(7) ^ 2) * t47, t23, -t61, -t51, t54, t52, t36, -t19 * t56 + t32 * t27, t20 * t56 + t32 * t29, -t19 * t29 - t20 * t27, t20 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1, t23, -t51, t61, t36, -t52, t54, t12 * t27 + t13 * t56, t13 * t29 - t14 * t27, -t12 * t29 - t14 * t56, t14 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1, t15, -t62, t11, t55, -t60, t34, t9 * t16 + t4 * t35, t9 * t18 - t5 * t35, -t5 * t16 - t4 * t18, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t15, t11, t62, t34, t60, t55, -t1 * t35 + t3 * t16, t1 * t18 - t2 * t16, -t3 * t18 + t2 * t35, t2 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
