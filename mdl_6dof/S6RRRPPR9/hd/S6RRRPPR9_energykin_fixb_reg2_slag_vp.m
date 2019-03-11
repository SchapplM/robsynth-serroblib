% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPPR9_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR9_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:17:53
% EndTime: 2019-03-09 16:17:53
% DurationCPUTime: 0.24s
% Computational Cost: add. (1046->67), mult. (2474->148), div. (0->0), fcn. (1896->10), ass. (0->57)
t77 = -pkin(4) - pkin(5);
t76 = cos(qJ(3));
t75 = cos(qJ(6));
t69 = cos(pkin(6)) * qJD(1);
t51 = qJD(2) + t69;
t57 = sin(qJ(3));
t58 = sin(qJ(2));
t54 = sin(pkin(6));
t70 = qJD(1) * t54;
t65 = t58 * t70;
t41 = t57 * t51 + t76 * t65;
t59 = cos(qJ(2));
t64 = t59 * t70;
t46 = -qJD(3) + t64;
t53 = sin(pkin(11));
t71 = cos(pkin(11));
t25 = t53 * t41 + t71 * t46;
t27 = t71 * t41 - t53 * t46;
t74 = t27 * t25;
t39 = -t76 * t51 + t57 * t65;
t73 = t39 * t25;
t60 = qJD(1) ^ 2;
t72 = t54 ^ 2 * t60;
t66 = pkin(1) * t69;
t42 = -pkin(8) * t65 + t59 * t66;
t32 = -t51 * pkin(2) - t42;
t15 = t39 * pkin(3) - t41 * qJ(4) + t32;
t43 = pkin(8) * t64 + t58 * t66;
t33 = t51 * pkin(9) + t43;
t34 = (-pkin(2) * t59 - pkin(9) * t58 - pkin(1)) * t70;
t22 = t76 * t33 + t57 * t34;
t19 = -t46 * qJ(4) + t22;
t10 = t53 * t15 + t71 * t19;
t21 = -t57 * t33 + t76 * t34;
t68 = t25 ^ 2 / 0.2e1;
t35 = t39 ^ 2 / 0.2e1;
t67 = t59 * t72;
t7 = t39 * qJ(5) + t10;
t63 = t72 / 0.2e1;
t9 = t71 * t15 - t53 * t19;
t18 = t46 * pkin(3) + qJD(4) - t21;
t62 = qJD(5) - t9;
t61 = t27 * qJ(5) - t18;
t56 = sin(qJ(6));
t37 = -qJD(6) + t39;
t24 = t27 ^ 2 / 0.2e1;
t20 = t27 * t39;
t13 = t56 * t25 + t75 * t27;
t11 = -t75 * t25 + t56 * t27;
t8 = t25 * pkin(4) - t61;
t6 = -t39 * pkin(4) + t62;
t5 = t77 * t25 + t61;
t4 = t25 * pkin(10) + t7;
t3 = -t27 * pkin(10) + t77 * t39 + t62;
t2 = t56 * t3 + t75 * t4;
t1 = t75 * t3 - t56 * t4;
t12 = [0, 0, 0, 0, 0, t60 / 0.2e1, 0, 0, 0, 0, t58 ^ 2 * t63, t58 * t67, t51 * t65, t59 ^ 2 * t63, t51 * t64, t51 ^ 2 / 0.2e1, pkin(1) * t67 + t42 * t51, -pkin(1) * t58 * t72 - t43 * t51 (-t42 * t58 + t43 * t59) * t70, t43 ^ 2 / 0.2e1 + t42 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t63, t41 ^ 2 / 0.2e1, -t41 * t39, -t41 * t46, t35, t39 * t46, t46 ^ 2 / 0.2e1, -t21 * t46 + t32 * t39, t22 * t46 + t32 * t41, -t21 * t41 - t22 * t39, t22 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1, t24, -t74, t20, t68, -t73, t35, t18 * t25 + t9 * t39, -t10 * t39 + t18 * t27, -t10 * t25 - t9 * t27, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t24, t20, t74, t35, t73, t68, t8 * t25 - t6 * t39, -t7 * t25 + t6 * t27, -t8 * t27 + t7 * t39, t7 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t13 * t11, -t13 * t37, t11 ^ 2 / 0.2e1, t11 * t37, t37 ^ 2 / 0.2e1, -t1 * t37 + t5 * t11, t5 * t13 + t2 * t37, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t12;
