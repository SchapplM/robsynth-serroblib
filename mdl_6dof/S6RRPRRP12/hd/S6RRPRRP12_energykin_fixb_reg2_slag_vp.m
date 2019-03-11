% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRP12_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP12_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:53:31
% EndTime: 2019-03-09 12:53:31
% DurationCPUTime: 0.16s
% Computational Cost: add. (536->60), mult. (1118->128), div. (0->0), fcn. (649->6), ass. (0->52)
t49 = qJD(1) ^ 2;
t64 = t49 / 0.2e1;
t63 = -pkin(2) - pkin(8);
t44 = sin(qJ(5));
t62 = cos(qJ(5));
t48 = cos(qJ(2));
t46 = sin(qJ(2));
t51 = -qJ(3) * t46 - pkin(1);
t21 = (t63 * t48 + t51) * qJD(1);
t57 = t46 * qJD(1);
t55 = pkin(7) * t57 + qJD(3);
t22 = pkin(3) * t57 + t63 * qJD(2) + t55;
t45 = sin(qJ(4));
t47 = cos(qJ(4));
t10 = -t45 * t21 + t47 * t22;
t58 = qJD(1) * t48;
t28 = t47 * qJD(2) - t45 * t58;
t34 = qJD(4) + t57;
t7 = t34 * pkin(4) - t28 * pkin(9) + t10;
t11 = t47 * t21 + t45 * t22;
t26 = t45 * qJD(2) + t47 * t58;
t9 = -t26 * pkin(9) + t11;
t4 = t44 * t7 + t62 * t9;
t14 = t62 * t26 + t44 * t28;
t16 = -t44 * t26 + t62 * t28;
t61 = t16 * t14;
t32 = qJD(5) + t34;
t60 = t32 * t14;
t59 = t48 * t49;
t31 = -pkin(7) * t58 - qJD(2) * qJ(3);
t56 = t14 ^ 2 / 0.2e1;
t54 = qJD(1) * qJD(2);
t24 = pkin(3) * t58 - t31;
t53 = t46 * t54;
t52 = t48 * t54;
t17 = t26 * pkin(4) + t24;
t3 = -t44 * t9 + t62 * t7;
t43 = t48 ^ 2;
t42 = t46 ^ 2;
t40 = qJD(2) ^ 2 / 0.2e1;
t36 = t43 * t64;
t35 = t42 * t64;
t33 = t46 * t59;
t30 = t32 ^ 2 / 0.2e1;
t29 = -qJD(2) * pkin(2) + t55;
t25 = (-pkin(2) * t48 + t51) * qJD(1);
t13 = t16 ^ 2 / 0.2e1;
t12 = t32 * t16;
t5 = t14 * pkin(5) - t16 * qJ(6) + t17;
t2 = t32 * qJ(6) + t4;
t1 = -t32 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t64, 0, 0, 0, 0, t35, t33, t53, t36, t52, t40, pkin(1) * t59 - pkin(7) * t53, -t49 * pkin(1) * t46 - pkin(7) * t52 (t42 + t43) * t49 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t43 / 0.2e1 + t42 / 0.2e1) * pkin(7) ^ 2) * t49, t40, -t53, -t52, t35, t33, t36 (t29 * t46 - t31 * t48) * qJD(1), t29 * qJD(2) + t25 * t58, -t31 * qJD(2) - t25 * t57, t25 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1, t28 ^ 2 / 0.2e1, -t28 * t26, t28 * t34, t26 ^ 2 / 0.2e1, -t26 * t34, t34 ^ 2 / 0.2e1, t10 * t34 + t24 * t26, -t11 * t34 + t24 * t28, -t10 * t28 - t11 * t26, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1, t13, -t61, t12, t56, -t60, t30, t17 * t14 + t3 * t32, t17 * t16 - t4 * t32, -t4 * t14 - t3 * t16, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t13, t12, t61, t30, t60, t56, -t1 * t32 + t5 * t14, t1 * t16 - t2 * t14, -t5 * t16 + t2 * t32, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
