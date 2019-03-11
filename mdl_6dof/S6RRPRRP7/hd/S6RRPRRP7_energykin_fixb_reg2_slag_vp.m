% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRRP7
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
% Datum: 2019-03-09 12:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRP7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP7_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:20:05
% EndTime: 2019-03-09 12:20:05
% DurationCPUTime: 0.18s
% Computational Cost: add. (491->60), mult. (1035->129), div. (0->0), fcn. (605->6), ass. (0->50)
t49 = sin(qJ(4));
t50 = sin(qJ(2));
t51 = cos(qJ(4));
t52 = cos(qJ(2));
t25 = (t49 * t50 + t51 * t52) * qJD(1);
t53 = qJD(1) ^ 2;
t67 = t53 / 0.2e1;
t62 = qJD(1) * t50;
t59 = pkin(7) * t62 + qJD(3);
t20 = -pkin(8) * t62 + (-pkin(2) - pkin(3)) * qJD(2) + t59;
t61 = qJD(1) * t52;
t31 = pkin(7) * t61 + qJD(2) * qJ(3);
t28 = -pkin(8) * t61 + t31;
t13 = t49 * t20 + t51 * t28;
t41 = qJD(2) - qJD(4);
t10 = -t41 * pkin(9) + t13;
t48 = sin(qJ(5));
t66 = cos(qJ(5));
t29 = -qJD(1) * pkin(1) - pkin(2) * t61 - qJ(3) * t62;
t19 = pkin(3) * t61 - t29;
t27 = (-t49 * t52 + t50 * t51) * qJD(1);
t7 = t25 * pkin(4) - t27 * pkin(9) + t19;
t4 = t66 * t10 + t48 * t7;
t15 = t48 * t27 + t66 * t41;
t17 = t66 * t27 - t48 * t41;
t65 = t17 * t15;
t24 = qJD(5) + t25;
t64 = t24 * t15;
t63 = t52 * t53;
t60 = t15 ^ 2 / 0.2e1;
t58 = qJD(1) * qJD(2);
t57 = t50 * t63;
t35 = t50 * t58;
t56 = t52 * t58;
t12 = t51 * t20 - t49 * t28;
t9 = t41 * pkin(4) - t12;
t3 = -t48 * t10 + t66 * t7;
t46 = t52 ^ 2;
t45 = t50 ^ 2;
t43 = qJD(2) ^ 2 / 0.2e1;
t34 = t46 * t67;
t33 = t45 * t67;
t30 = -qJD(2) * pkin(2) + t59;
t21 = t24 ^ 2 / 0.2e1;
t14 = t17 ^ 2 / 0.2e1;
t11 = t17 * t24;
t5 = t15 * pkin(5) - t17 * qJ(6) + t9;
t2 = t24 * qJ(6) + t4;
t1 = -t24 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t67, 0, 0, 0, 0, t33, t57, t35, t34, t56, t43, pkin(1) * t63 - pkin(7) * t35, -t53 * pkin(1) * t50 - pkin(7) * t56 (t45 + t46) * t53 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t46 / 0.2e1 + t45 / 0.2e1) * pkin(7) ^ 2) * t53, t33, t35, -t57, t43, -t56, t34, -t30 * qJD(2) - t29 * t61 (t30 * t50 + t31 * t52) * qJD(1), t31 * qJD(2) - t29 * t62, t31 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1, t27 ^ 2 / 0.2e1, -t27 * t25, -t27 * t41, t25 ^ 2 / 0.2e1, t25 * t41, t41 ^ 2 / 0.2e1, -t12 * t41 + t19 * t25, t13 * t41 + t19 * t27, -t12 * t27 - t13 * t25, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t14, -t65, t11, t60, -t64, t21, t9 * t15 + t3 * t24, t9 * t17 - t4 * t24, -t4 * t15 - t3 * t17, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t14, t11, t65, t21, t64, t60, -t1 * t24 + t5 * t15, t1 * t17 - t2 * t15, -t5 * t17 + t2 * t24, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
