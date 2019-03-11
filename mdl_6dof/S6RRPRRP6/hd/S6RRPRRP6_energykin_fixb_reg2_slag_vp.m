% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRP6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP6_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:14:08
% EndTime: 2019-03-09 12:14:08
% DurationCPUTime: 0.22s
% Computational Cost: add. (1055->67), mult. (2940->148), div. (0->0), fcn. (2299->10), ass. (0->56)
t59 = cos(qJ(2));
t67 = cos(pkin(6)) * qJD(1);
t64 = pkin(1) * t67;
t49 = t59 * t64;
t50 = qJD(2) + t67;
t58 = sin(qJ(2));
t53 = sin(pkin(6));
t68 = qJD(1) * t53;
t63 = t58 * t68;
t32 = t50 * pkin(2) + t49 + (-pkin(8) - qJ(3)) * t63;
t62 = t59 * t68;
t42 = pkin(8) * t62 + t58 * t64;
t35 = qJ(3) * t62 + t42;
t52 = sin(pkin(11));
t54 = cos(pkin(11));
t18 = t54 * t32 - t52 * t35;
t16 = -t50 * pkin(3) - t18;
t40 = (t52 * t59 + t54 * t58) * t68;
t57 = sin(qJ(4));
t73 = cos(qJ(4));
t29 = t57 * t40 - t73 * t50;
t31 = t73 * t40 + t57 * t50;
t10 = t29 * pkin(4) - t31 * pkin(10) + t16;
t56 = sin(qJ(5));
t72 = cos(qJ(5));
t19 = t52 * t32 + t54 * t35;
t17 = t50 * pkin(9) + t19;
t38 = t52 * t63 - t54 * t62;
t43 = qJD(3) + (-pkin(2) * t59 - pkin(1)) * t68;
t25 = t38 * pkin(3) - t40 * pkin(9) + t43;
t12 = t73 * t17 + t57 * t25;
t37 = qJD(4) + t38;
t8 = t37 * pkin(10) + t12;
t4 = t56 * t10 + t72 * t8;
t21 = t56 * t31 - t72 * t37;
t23 = t72 * t31 + t56 * t37;
t71 = t23 * t21;
t28 = qJD(5) + t29;
t70 = t28 * t21;
t60 = qJD(1) ^ 2;
t69 = t53 ^ 2 * t60;
t66 = t21 ^ 2 / 0.2e1;
t65 = t59 * t69;
t61 = t69 / 0.2e1;
t11 = -t57 * t17 + t73 * t25;
t3 = t72 * t10 - t56 * t8;
t7 = -t37 * pkin(4) - t11;
t46 = t50 ^ 2 / 0.2e1;
t41 = -pkin(8) * t63 + t49;
t26 = t28 ^ 2 / 0.2e1;
t20 = t23 ^ 2 / 0.2e1;
t13 = t23 * t28;
t5 = t21 * pkin(5) - t23 * qJ(6) + t7;
t2 = t28 * qJ(6) + t4;
t1 = -t28 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t60 / 0.2e1, 0, 0, 0, 0, t58 ^ 2 * t61, t58 * t65, t50 * t63, t59 ^ 2 * t61, t50 * t62, t46, pkin(1) * t65 + t41 * t50, -pkin(1) * t58 * t69 - t42 * t50 (-t41 * t58 + t42 * t59) * t68, t42 ^ 2 / 0.2e1 + t41 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t61, t40 ^ 2 / 0.2e1, -t40 * t38, t40 * t50, t38 ^ 2 / 0.2e1, -t38 * t50, t46, t18 * t50 + t43 * t38, -t19 * t50 + t43 * t40, -t18 * t40 - t19 * t38, t19 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1 + t43 ^ 2 / 0.2e1, t31 ^ 2 / 0.2e1, -t31 * t29, t31 * t37, t29 ^ 2 / 0.2e1, -t29 * t37, t37 ^ 2 / 0.2e1, t11 * t37 + t16 * t29, -t12 * t37 + t16 * t31, -t11 * t31 - t12 * t29, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1, t20, -t71, t13, t66, -t70, t26, t7 * t21 + t3 * t28, t7 * t23 - t4 * t28, -t4 * t21 - t3 * t23, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t20, t13, t71, t26, t70, t66, -t1 * t28 + t5 * t21, t1 * t23 - t2 * t21, t2 * t28 - t5 * t23, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
