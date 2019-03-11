% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRR8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR8_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:52:20
% EndTime: 2019-03-09 18:52:20
% DurationCPUTime: 0.26s
% Computational Cost: add. (1673->72), mult. (3927->167), div. (0->0), fcn. (3136->12), ass. (0->56)
t65 = sin(qJ(2));
t66 = cos(qJ(2));
t59 = sin(pkin(6));
t74 = qJD(1) * t59;
t69 = t66 * t74;
t73 = cos(pkin(6)) * qJD(1);
t71 = pkin(1) * t73;
t47 = pkin(8) * t69 + t65 * t71;
t56 = qJD(2) + t73;
t40 = pkin(9) * t56 + t47;
t42 = (-pkin(2) * t66 - pkin(9) * t65 - pkin(1)) * t74;
t64 = sin(qJ(3));
t78 = cos(qJ(3));
t28 = -t40 * t64 + t42 * t78;
t70 = t65 * t74;
t45 = t56 * t64 + t70 * t78;
t51 = -qJD(3) + t69;
t20 = -pkin(3) * t51 - qJ(4) * t45 + t28;
t29 = t40 * t78 + t42 * t64;
t43 = -t56 * t78 + t64 * t70;
t23 = -qJ(4) * t43 + t29;
t58 = sin(pkin(12));
t60 = cos(pkin(12));
t12 = t20 * t58 + t23 * t60;
t10 = -pkin(10) * t51 + t12;
t46 = -pkin(8) * t70 + t66 * t71;
t39 = -pkin(2) * t56 - t46;
t32 = pkin(3) * t43 + qJD(4) + t39;
t33 = t43 * t60 + t45 * t58;
t35 = -t43 * t58 + t45 * t60;
t15 = pkin(4) * t33 - pkin(10) * t35 + t32;
t63 = sin(qJ(5));
t77 = cos(qJ(5));
t6 = t10 * t77 + t15 * t63;
t76 = cos(qJ(6));
t67 = qJD(1) ^ 2;
t75 = t59 ^ 2 * t67;
t72 = t66 * t75;
t68 = t75 / 0.2e1;
t5 = -t10 * t63 + t15 * t77;
t11 = t20 * t60 - t23 * t58;
t31 = qJD(5) + t33;
t9 = pkin(4) * t51 - t11;
t62 = sin(qJ(6));
t49 = t51 ^ 2 / 0.2e1;
t30 = qJD(6) + t31;
t27 = t35 * t77 - t51 * t63;
t25 = t35 * t63 + t51 * t77;
t18 = -t25 * t62 + t27 * t76;
t16 = t25 * t76 + t27 * t62;
t7 = pkin(5) * t25 + t9;
t4 = -pkin(11) * t25 + t6;
t3 = pkin(5) * t31 - pkin(11) * t27 + t5;
t2 = t3 * t62 + t4 * t76;
t1 = t3 * t76 - t4 * t62;
t8 = [0, 0, 0, 0, 0, t67 / 0.2e1, 0, 0, 0, 0, t65 ^ 2 * t68, t65 * t72, t56 * t70, t66 ^ 2 * t68, t56 * t69, t56 ^ 2 / 0.2e1, pkin(1) * t72 + t46 * t56, -pkin(1) * t65 * t75 - t47 * t56 (-t46 * t65 + t47 * t66) * t74, t47 ^ 2 / 0.2e1 + t46 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t68, t45 ^ 2 / 0.2e1, -t45 * t43, -t45 * t51, t43 ^ 2 / 0.2e1, t43 * t51, t49, -t28 * t51 + t39 * t43, t29 * t51 + t39 * t45, -t28 * t45 - t29 * t43, t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t39 ^ 2 / 0.2e1, t35 ^ 2 / 0.2e1, -t35 * t33, -t35 * t51, t33 ^ 2 / 0.2e1, t33 * t51, t49, -t11 * t51 + t32 * t33, t12 * t51 + t32 * t35, -t11 * t35 - t12 * t33, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t31, t25 ^ 2 / 0.2e1, -t25 * t31, t31 ^ 2 / 0.2e1, t25 * t9 + t31 * t5, t27 * t9 - t31 * t6, -t25 * t6 - t27 * t5, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t30, t16 ^ 2 / 0.2e1, -t16 * t30, t30 ^ 2 / 0.2e1, t1 * t30 + t16 * t7, t18 * t7 - t2 * t30, -t1 * t18 - t16 * t2, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
