% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRP9_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP9_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:11:38
% EndTime: 2019-03-10 02:11:38
% DurationCPUTime: 0.25s
% Computational Cost: add. (1311->68), mult. (3019->150), div. (0->0), fcn. (2366->10), ass. (0->55)
t65 = cos(qJ(2));
t62 = sin(qJ(2));
t57 = sin(pkin(6));
t73 = qJD(1) * t57;
t69 = t62 * t73;
t72 = cos(pkin(6)) * qJD(1);
t70 = pkin(1) * t72;
t46 = -pkin(8) * t69 + t65 * t70;
t55 = qJD(2) + t72;
t37 = -t55 * pkin(2) - t46;
t61 = sin(qJ(3));
t64 = cos(qJ(3));
t43 = -t64 * t55 + t61 * t69;
t45 = t61 * t55 + t64 * t69;
t24 = t43 * pkin(3) - t45 * pkin(10) + t37;
t68 = t65 * t73;
t47 = pkin(8) * t68 + t62 * t70;
t38 = t55 * pkin(9) + t47;
t41 = (-pkin(2) * t65 - pkin(9) * t62 - pkin(1)) * t73;
t29 = t64 * t38 + t61 * t41;
t50 = -qJD(3) + t68;
t27 = -t50 * pkin(10) + t29;
t60 = sin(qJ(4));
t63 = cos(qJ(4));
t13 = t60 * t24 + t63 * t27;
t31 = t60 * t45 + t63 * t50;
t10 = -t31 * pkin(11) + t13;
t59 = sin(qJ(5));
t12 = t63 * t24 - t60 * t27;
t33 = t63 * t45 - t60 * t50;
t42 = qJD(4) + t43;
t7 = t42 * pkin(4) - t33 * pkin(11) + t12;
t75 = cos(qJ(5));
t4 = t75 * t10 + t59 * t7;
t66 = qJD(1) ^ 2;
t74 = t57 ^ 2 * t66;
t71 = t65 * t74;
t67 = t74 / 0.2e1;
t3 = -t59 * t10 + t75 * t7;
t28 = -t61 * t38 + t64 * t41;
t26 = t50 * pkin(3) - t28;
t16 = t31 * pkin(4) + t26;
t40 = qJD(5) + t42;
t39 = t40 ^ 2 / 0.2e1;
t21 = -t59 * t31 + t75 * t33;
t19 = t75 * t31 + t59 * t33;
t18 = t21 ^ 2 / 0.2e1;
t17 = t19 ^ 2 / 0.2e1;
t15 = t21 * t40;
t14 = t19 * t40;
t11 = t21 * t19;
t8 = t19 * pkin(5) + qJD(6) + t16;
t2 = -t19 * qJ(6) + t4;
t1 = t40 * pkin(5) - t21 * qJ(6) + t3;
t5 = [0, 0, 0, 0, 0, t66 / 0.2e1, 0, 0, 0, 0, t62 ^ 2 * t67, t62 * t71, t55 * t69, t65 ^ 2 * t67, t55 * t68, t55 ^ 2 / 0.2e1, pkin(1) * t71 + t46 * t55, -pkin(1) * t62 * t74 - t47 * t55 (-t46 * t62 + t47 * t65) * t73, t47 ^ 2 / 0.2e1 + t46 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t67, t45 ^ 2 / 0.2e1, -t45 * t43, -t45 * t50, t43 ^ 2 / 0.2e1, t43 * t50, t50 ^ 2 / 0.2e1, -t28 * t50 + t37 * t43, t29 * t50 + t37 * t45, -t28 * t45 - t29 * t43, t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t37 ^ 2 / 0.2e1, t33 ^ 2 / 0.2e1, -t33 * t31, t33 * t42, t31 ^ 2 / 0.2e1, -t31 * t42, t42 ^ 2 / 0.2e1, t12 * t42 + t26 * t31, -t13 * t42 + t26 * t33, -t12 * t33 - t13 * t31, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1, t18, -t11, t15, t17, -t14, t39, t16 * t19 + t3 * t40, t16 * t21 - t4 * t40, -t4 * t19 - t3 * t21, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1, t18, -t11, t15, t17, -t14, t39, t1 * t40 + t8 * t19, -t2 * t40 + t8 * t21, -t1 * t21 - t2 * t19, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1;];
T_reg  = t5;
