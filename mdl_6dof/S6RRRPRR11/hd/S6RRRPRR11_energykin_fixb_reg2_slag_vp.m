% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRR11_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR11_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:34:13
% EndTime: 2019-03-09 19:34:13
% DurationCPUTime: 0.22s
% Computational Cost: add. (982->69), mult. (2282->150), div. (0->0), fcn. (1719->10), ass. (0->56)
t60 = sin(qJ(2));
t61 = cos(qJ(2));
t55 = sin(pkin(6));
t71 = qJD(1) * t55;
t65 = t61 * t71;
t70 = cos(pkin(6)) * qJD(1);
t67 = pkin(1) * t70;
t40 = pkin(8) * t65 + t60 * t67;
t53 = qJD(2) + t70;
t30 = t53 * pkin(9) + t40;
t32 = (-pkin(2) * t61 - pkin(9) * t60 - pkin(1)) * t71;
t59 = sin(qJ(3));
t77 = cos(qJ(3));
t20 = t77 * t30 + t59 * t32;
t46 = -qJD(3) + t65;
t15 = -t46 * qJ(4) + t20;
t66 = t60 * t71;
t36 = -t77 * t53 + t59 * t66;
t12 = t36 * pkin(10) + t15;
t58 = sin(qJ(5));
t76 = cos(qJ(5));
t38 = t59 * t53 + t77 * t66;
t19 = -t59 * t30 + t77 * t32;
t63 = qJD(4) - t19;
t9 = -t38 * pkin(10) + (pkin(3) + pkin(4)) * t46 + t63;
t6 = t76 * t12 + t58 * t9;
t75 = cos(qJ(6));
t74 = t38 * t36;
t73 = t46 * t36;
t62 = qJD(1) ^ 2;
t72 = t55 ^ 2 * t62;
t39 = -pkin(8) * t66 + t61 * t67;
t69 = t36 ^ 2 / 0.2e1;
t68 = t61 * t72;
t29 = -t53 * pkin(2) - t39;
t64 = t72 / 0.2e1;
t22 = -t76 * t36 + t58 * t38;
t13 = t36 * pkin(3) - t38 * qJ(4) + t29;
t5 = -t58 * t12 + t76 * t9;
t10 = -t36 * pkin(4) - t13;
t57 = sin(qJ(6));
t44 = qJD(5) + t46;
t42 = t46 ^ 2 / 0.2e1;
t33 = t38 ^ 2 / 0.2e1;
t25 = t38 * t46;
t24 = t58 * t36 + t76 * t38;
t21 = qJD(6) + t22;
t18 = t75 * t24 + t57 * t44;
t16 = t57 * t24 - t75 * t44;
t14 = t46 * pkin(3) + t63;
t7 = t22 * pkin(5) - t24 * pkin(11) + t10;
t4 = t44 * pkin(11) + t6;
t3 = -t44 * pkin(5) - t5;
t2 = t75 * t4 + t57 * t7;
t1 = -t57 * t4 + t75 * t7;
t8 = [0, 0, 0, 0, 0, t62 / 0.2e1, 0, 0, 0, 0, t60 ^ 2 * t64, t60 * t68, t53 * t66, t61 ^ 2 * t64, t53 * t65, t53 ^ 2 / 0.2e1, pkin(1) * t68 + t39 * t53, -pkin(1) * t60 * t72 - t40 * t53 (-t39 * t60 + t40 * t61) * t71, t40 ^ 2 / 0.2e1 + t39 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t64, t33, -t74, -t25, t69, t73, t42, -t19 * t46 + t29 * t36, t20 * t46 + t29 * t38, -t19 * t38 - t20 * t36, t20 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1, t33, -t25, t74, t42, -t73, t69, t13 * t36 + t14 * t46, t14 * t38 - t15 * t36, -t13 * t38 - t15 * t46, t15 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t24 ^ 2 / 0.2e1, -t24 * t22, t24 * t44, t22 ^ 2 / 0.2e1, -t22 * t44, t44 ^ 2 / 0.2e1, t10 * t22 + t5 * t44, t10 * t24 - t6 * t44, -t6 * t22 - t5 * t24, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t21, t16 ^ 2 / 0.2e1, -t16 * t21, t21 ^ 2 / 0.2e1, t1 * t21 + t3 * t16, t3 * t18 - t2 * t21, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
