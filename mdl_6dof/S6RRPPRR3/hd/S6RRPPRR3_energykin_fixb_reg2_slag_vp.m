% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPRR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:58:36
% EndTime: 2019-03-09 08:58:36
% DurationCPUTime: 0.26s
% Computational Cost: add. (1380->73), mult. (3882->165), div. (0->0), fcn. (3100->12), ass. (0->57)
t63 = cos(qJ(2));
t71 = cos(pkin(6)) * qJD(1);
t68 = pkin(1) * t71;
t52 = t63 * t68;
t53 = qJD(2) + t71;
t62 = sin(qJ(2));
t57 = sin(pkin(6));
t72 = qJD(1) * t57;
t67 = t62 * t72;
t35 = t53 * pkin(2) + t52 + (-pkin(8) - qJ(3)) * t67;
t66 = t63 * t72;
t45 = pkin(8) * t66 + t62 * t68;
t38 = qJ(3) * t66 + t45;
t56 = sin(pkin(11));
t58 = cos(pkin(11));
t26 = t56 * t35 + t58 * t38;
t24 = t53 * qJ(4) + t26;
t41 = t56 * t67 - t58 * t66;
t43 = (t56 * t63 + t58 * t62) * t72;
t46 = qJD(3) + (-pkin(2) * t63 - pkin(1)) * t72;
t29 = t41 * pkin(3) - t43 * qJ(4) + t46;
t55 = sin(pkin(12));
t73 = cos(pkin(12));
t13 = t73 * t24 + t55 * t29;
t32 = t55 * t43 - t73 * t53;
t11 = -t32 * pkin(9) + t13;
t61 = sin(qJ(5));
t76 = cos(qJ(5));
t12 = -t55 * t24 + t73 * t29;
t34 = t73 * t43 + t55 * t53;
t9 = t41 * pkin(4) - t34 * pkin(9) + t12;
t6 = t76 * t11 + t61 * t9;
t75 = cos(qJ(6));
t64 = qJD(1) ^ 2;
t74 = t57 ^ 2 * t64;
t70 = t41 ^ 2 / 0.2e1;
t69 = t63 * t74;
t65 = t74 / 0.2e1;
t20 = t76 * t32 + t61 * t34;
t25 = t58 * t35 - t56 * t38;
t5 = -t61 * t11 + t76 * t9;
t23 = -t53 * pkin(3) + qJD(4) - t25;
t14 = t32 * pkin(4) + t23;
t60 = sin(qJ(6));
t49 = t53 ^ 2 / 0.2e1;
t44 = -pkin(8) * t67 + t52;
t40 = qJD(5) + t41;
t22 = -t61 * t32 + t76 * t34;
t18 = qJD(6) + t20;
t17 = t75 * t22 + t60 * t40;
t15 = t60 * t22 - t75 * t40;
t7 = t20 * pkin(5) - t22 * pkin(10) + t14;
t4 = t40 * pkin(10) + t6;
t3 = -t40 * pkin(5) - t5;
t2 = t75 * t4 + t60 * t7;
t1 = -t60 * t4 + t75 * t7;
t8 = [0, 0, 0, 0, 0, t64 / 0.2e1, 0, 0, 0, 0, t62 ^ 2 * t65, t62 * t69, t53 * t67, t63 ^ 2 * t65, t53 * t66, t49, pkin(1) * t69 + t44 * t53, -pkin(1) * t62 * t74 - t45 * t53 (-t44 * t62 + t45 * t63) * t72, t45 ^ 2 / 0.2e1 + t44 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t65, t43 ^ 2 / 0.2e1, -t43 * t41, t43 * t53, t70, -t41 * t53, t49, t25 * t53 + t46 * t41, -t26 * t53 + t46 * t43, -t25 * t43 - t26 * t41, t26 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1 + t46 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, t34 * t41, t32 ^ 2 / 0.2e1, -t32 * t41, t70, t12 * t41 + t23 * t32, -t13 * t41 + t23 * t34, -t12 * t34 - t13 * t32, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * t40, t20 ^ 2 / 0.2e1, -t20 * t40, t40 ^ 2 / 0.2e1, t14 * t20 + t5 * t40, t14 * t22 - t6 * t40, -t6 * t20 - t5 * t22, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t17 ^ 2 / 0.2e1, -t17 * t15, t17 * t18, t15 ^ 2 / 0.2e1, -t15 * t18, t18 ^ 2 / 0.2e1, t1 * t18 + t3 * t15, t3 * t17 - t2 * t18, -t1 * t17 - t2 * t15, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
