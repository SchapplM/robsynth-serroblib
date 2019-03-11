% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 04:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRR6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR6_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:21:33
% EndTime: 2019-03-10 04:21:34
% DurationCPUTime: 0.26s
% Computational Cost: add. (1726->72), mult. (3927->170), div. (0->0), fcn. (3136->12), ass. (0->56)
t65 = sin(qJ(2));
t67 = cos(qJ(2));
t59 = sin(pkin(6));
t75 = qJD(1) * t59;
t70 = t67 * t75;
t74 = cos(pkin(6)) * qJD(1);
t72 = pkin(1) * t74;
t47 = pkin(8) * t70 + t65 * t72;
t57 = qJD(2) + t74;
t40 = pkin(9) * t57 + t47;
t42 = (-pkin(2) * t67 - pkin(9) * t65 - pkin(1)) * t75;
t64 = sin(qJ(3));
t79 = cos(qJ(3));
t28 = -t40 * t64 + t42 * t79;
t71 = t65 * t75;
t45 = t57 * t64 + t71 * t79;
t52 = -qJD(3) + t70;
t20 = -pkin(3) * t52 - pkin(10) * t45 + t28;
t29 = t40 * t79 + t42 * t64;
t43 = -t57 * t79 + t64 * t71;
t23 = -pkin(10) * t43 + t29;
t63 = sin(qJ(4));
t66 = cos(qJ(4));
t12 = t20 * t63 + t23 * t66;
t49 = -qJD(4) + t52;
t10 = -pkin(11) * t49 + t12;
t32 = t43 * t66 + t45 * t63;
t34 = -t43 * t63 + t45 * t66;
t46 = -pkin(8) * t71 + t67 * t72;
t39 = -pkin(2) * t57 - t46;
t35 = pkin(3) * t43 + t39;
t15 = pkin(4) * t32 - pkin(11) * t34 + t35;
t62 = sin(qJ(5));
t78 = cos(qJ(5));
t6 = t10 * t78 + t15 * t62;
t77 = cos(qJ(6));
t68 = qJD(1) ^ 2;
t76 = t59 ^ 2 * t68;
t73 = t67 * t76;
t69 = t76 / 0.2e1;
t5 = -t10 * t62 + t15 * t78;
t11 = t20 * t66 - t23 * t63;
t31 = qJD(5) + t32;
t9 = pkin(4) * t49 - t11;
t61 = sin(qJ(6));
t30 = qJD(6) + t31;
t27 = t34 * t78 - t49 * t62;
t25 = t34 * t62 + t49 * t78;
t18 = -t25 * t61 + t27 * t77;
t16 = t25 * t77 + t27 * t61;
t7 = pkin(5) * t25 + t9;
t4 = -pkin(12) * t25 + t6;
t3 = pkin(5) * t31 - pkin(12) * t27 + t5;
t2 = t3 * t61 + t4 * t77;
t1 = t3 * t77 - t4 * t61;
t8 = [0, 0, 0, 0, 0, t68 / 0.2e1, 0, 0, 0, 0, t65 ^ 2 * t69, t65 * t73, t57 * t71, t67 ^ 2 * t69, t57 * t70, t57 ^ 2 / 0.2e1, pkin(1) * t73 + t46 * t57, -pkin(1) * t65 * t76 - t47 * t57 (-t46 * t65 + t47 * t67) * t75, t47 ^ 2 / 0.2e1 + t46 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t69, t45 ^ 2 / 0.2e1, -t45 * t43, -t45 * t52, t43 ^ 2 / 0.2e1, t43 * t52, t52 ^ 2 / 0.2e1, -t28 * t52 + t39 * t43, t29 * t52 + t39 * t45, -t28 * t45 - t29 * t43, t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t39 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, -t34 * t49, t32 ^ 2 / 0.2e1, t32 * t49, t49 ^ 2 / 0.2e1, -t11 * t49 + t32 * t35, t12 * t49 + t34 * t35, -t11 * t34 - t12 * t32, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t35 ^ 2 / 0.2e1, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t31, t25 ^ 2 / 0.2e1, -t25 * t31, t31 ^ 2 / 0.2e1, t25 * t9 + t31 * t5, t27 * t9 - t31 * t6, -t25 * t6 - t27 * t5, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t30, t16 ^ 2 / 0.2e1, -t16 * t30, t30 ^ 2 / 0.2e1, t1 * t30 + t16 * t7, t18 * t7 - t2 * t30, -t1 * t18 - t16 * t2, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
