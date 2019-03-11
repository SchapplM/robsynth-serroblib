% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRR11_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR11_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_energykin_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:46:03
% EndTime: 2019-03-09 07:46:03
% DurationCPUTime: 0.37s
% Computational Cost: add. (2223->82), mult. (7110->197), div. (0->0), fcn. (5932->14), ass. (0->66)
t63 = sin(pkin(13));
t66 = cos(pkin(13));
t65 = sin(pkin(6));
t83 = qJD(1) * t65;
t77 = qJ(2) * t83;
t68 = cos(pkin(6));
t82 = qJD(1) * t68;
t80 = pkin(1) * t82;
t54 = t63 * t80 + t66 * t77;
t64 = sin(pkin(7));
t67 = cos(pkin(7));
t85 = t65 * t66;
t41 = (t64 * t68 + t67 * t85) * qJD(1) * pkin(9) + t54;
t59 = t66 * t80;
t87 = t63 * t65;
t43 = t59 + (pkin(2) * t68 + (-pkin(9) * t67 - qJ(2)) * t87) * qJD(1);
t49 = qJD(2) + (-pkin(9) * t63 * t64 - pkin(2) * t66 - pkin(1)) * t83;
t72 = sin(qJ(3));
t74 = cos(qJ(3));
t25 = -t72 * t41 + (t43 * t67 + t49 * t64) * t74;
t75 = qJD(1) ^ 2;
t91 = t75 / 0.2e1;
t31 = -t64 * t43 + t67 * t49;
t79 = t66 * t83;
t44 = t72 * t63 * t83 + (-t64 * t82 - t67 * t79) * t74;
t84 = t67 * t72;
t86 = t64 * t72;
t46 = (t68 * t86 + (t63 * t74 + t66 * t84) * t65) * qJD(1);
t20 = t44 * pkin(3) - t46 * pkin(10) + t31;
t26 = t74 * t41 + t43 * t84 + t49 * t86;
t51 = t64 * t79 - t67 * t82 - qJD(3);
t24 = -t51 * pkin(10) + t26;
t71 = sin(qJ(4));
t73 = cos(qJ(4));
t12 = t71 * t20 + t73 * t24;
t42 = qJD(4) + t44;
t10 = t42 * pkin(11) + t12;
t23 = t51 * pkin(3) - t25;
t34 = t71 * t46 + t73 * t51;
t36 = t73 * t46 - t71 * t51;
t15 = t34 * pkin(4) - t36 * pkin(11) + t23;
t70 = sin(qJ(5));
t90 = cos(qJ(5));
t6 = t90 * t10 + t70 * t15;
t89 = cos(qJ(6));
t88 = t65 ^ 2 * t75;
t81 = t65 * t68 * t75;
t78 = t88 / 0.2e1;
t5 = -t70 * t10 + t90 * t15;
t11 = t73 * t20 - t71 * t24;
t33 = qJD(5) + t34;
t9 = -t42 * pkin(4) - t11;
t69 = sin(qJ(6));
t60 = -pkin(1) * t83 + qJD(2);
t53 = -t63 * t77 + t59;
t32 = qJD(6) + t33;
t30 = t90 * t36 + t70 * t42;
t28 = t70 * t36 - t90 * t42;
t18 = -t69 * t28 + t89 * t30;
t16 = t89 * t28 + t69 * t30;
t7 = t28 * pkin(5) + t9;
t4 = -t28 * pkin(12) + t6;
t3 = t33 * pkin(5) - t30 * pkin(12) + t5;
t2 = t69 * t3 + t89 * t4;
t1 = t89 * t3 - t69 * t4;
t8 = [0, 0, 0, 0, 0, t91, 0, 0, 0, 0, t63 ^ 2 * t78, t63 * t66 * t88, t63 * t81, t66 ^ 2 * t78, t66 * t81, t68 ^ 2 * t91 (t53 * t68 - t60 * t85) * qJD(1) (-t54 * t68 + t60 * t87) * qJD(1) (-t53 * t63 + t54 * t66) * t83, t54 ^ 2 / 0.2e1 + t53 ^ 2 / 0.2e1 + t60 ^ 2 / 0.2e1, t46 ^ 2 / 0.2e1, -t46 * t44, -t46 * t51, t44 ^ 2 / 0.2e1, t44 * t51, t51 ^ 2 / 0.2e1, -t25 * t51 + t31 * t44, t26 * t51 + t31 * t46, -t25 * t46 - t26 * t44, t26 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1, t36 ^ 2 / 0.2e1, -t36 * t34, t36 * t42, t34 ^ 2 / 0.2e1, -t34 * t42, t42 ^ 2 / 0.2e1, t11 * t42 + t23 * t34, -t12 * t42 + t23 * t36, -t11 * t36 - t12 * t34, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1, t30 ^ 2 / 0.2e1, -t30 * t28, t30 * t33, t28 ^ 2 / 0.2e1, -t28 * t33, t33 ^ 2 / 0.2e1, t9 * t28 + t5 * t33, t9 * t30 - t6 * t33, -t6 * t28 - t5 * t30, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t32, t16 ^ 2 / 0.2e1, -t16 * t32, t32 ^ 2 / 0.2e1, t1 * t32 + t7 * t16, t7 * t18 - t2 * t32, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
