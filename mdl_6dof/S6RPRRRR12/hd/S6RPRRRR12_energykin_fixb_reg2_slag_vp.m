% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRR12_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR12_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_energykin_fixb_reg2_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:00:39
% EndTime: 2019-03-09 08:00:40
% DurationCPUTime: 0.56s
% Computational Cost: add. (3555->86), mult. (11840->204), div. (0->0), fcn. (10104->16), ass. (0->71)
t103 = cos(qJ(3));
t70 = sin(pkin(7));
t74 = cos(pkin(7));
t75 = cos(pkin(6));
t71 = sin(pkin(6));
t72 = cos(pkin(14));
t92 = t71 * t72;
t107 = qJD(1) * (t70 * t75 + t74 * t92);
t68 = sin(pkin(14));
t79 = sin(qJ(3));
t90 = qJD(1) * t71;
t48 = t79 * t68 * t90 - t103 * t107;
t57 = -qJD(3) - (-t70 * t92 + t74 * t75) * qJD(1);
t69 = sin(pkin(8));
t73 = cos(pkin(8));
t109 = t48 * t73 + t57 * t69;
t91 = t74 * t79;
t93 = t70 * t79;
t50 = (t75 * t93 + (t103 * t68 + t72 * t91) * t71) * qJD(1);
t104 = pkin(11) * t50;
t85 = qJ(2) * t90;
t87 = pkin(1) * qJD(1) * t75;
t60 = t68 * t87 + t72 * t85;
t46 = pkin(10) * t107 + t60;
t65 = t72 * t87;
t95 = t68 * t71;
t47 = t65 + (pkin(2) * t75 + (-pkin(10) * t74 - qJ(2)) * t95) * qJD(1);
t54 = qJD(2) + (-pkin(10) * t68 * t70 - pkin(2) * t72 - pkin(1)) * t90;
t31 = -t79 * t46 + (t47 * t74 + t54 * t70) * t103;
t23 = -t57 * pkin(3) - t73 * t104 + t31;
t38 = -t70 * t47 + t74 * t54;
t30 = t48 * pkin(3) - t69 * t104 + t38;
t108 = t23 * t73 + t30 * t69;
t32 = t103 * t46 + t47 * t91 + t54 * t93;
t22 = -pkin(11) * t109 + t32;
t78 = sin(qJ(4));
t80 = cos(qJ(4));
t13 = t108 * t80 - t78 * t22;
t81 = qJD(1) ^ 2;
t105 = t81 / 0.2e1;
t14 = t108 * t78 + t80 * t22;
t40 = -t69 * t48 + t73 * t57 - qJD(4);
t10 = -t40 * pkin(12) + t14;
t102 = cos(qJ(5));
t15 = -t69 * t23 + t73 * t30;
t35 = t109 * t80 + t78 * t50;
t37 = -t109 * t78 + t80 * t50;
t12 = t35 * pkin(4) - t37 * pkin(12) + t15;
t77 = sin(qJ(5));
t6 = t102 * t10 + t77 * t12;
t101 = cos(qJ(6));
t96 = t71 ^ 2 * t81;
t88 = t71 * t75 * t81;
t86 = t96 / 0.2e1;
t25 = t102 * t40 + t77 * t37;
t5 = -t77 * t10 + t102 * t12;
t9 = t40 * pkin(4) - t13;
t76 = sin(qJ(6));
t66 = -pkin(1) * t90 + qJD(2);
t59 = -t68 * t85 + t65;
t34 = qJD(5) + t35;
t27 = t102 * t37 - t77 * t40;
t24 = qJD(6) + t25;
t18 = t101 * t27 + t76 * t34;
t16 = -t101 * t34 + t76 * t27;
t7 = t25 * pkin(5) - t27 * pkin(13) + t9;
t4 = t34 * pkin(13) + t6;
t3 = -t34 * pkin(5) - t5;
t2 = t101 * t4 + t76 * t7;
t1 = t101 * t7 - t76 * t4;
t8 = [0, 0, 0, 0, 0, t105, 0, 0, 0, 0, t68 ^ 2 * t86, t68 * t72 * t96, t68 * t88, t72 ^ 2 * t86, t72 * t88, t75 ^ 2 * t105 (t59 * t75 - t66 * t92) * qJD(1) (-t60 * t75 + t66 * t95) * qJD(1) (-t59 * t68 + t60 * t72) * t90, t60 ^ 2 / 0.2e1 + t59 ^ 2 / 0.2e1 + t66 ^ 2 / 0.2e1, t50 ^ 2 / 0.2e1, -t50 * t48, -t57 * t50, t48 ^ 2 / 0.2e1, t57 * t48, t57 ^ 2 / 0.2e1, -t31 * t57 + t38 * t48, t32 * t57 + t38 * t50, -t31 * t50 - t32 * t48, t32 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1, t37 ^ 2 / 0.2e1, -t35 * t37, -t40 * t37, t35 ^ 2 / 0.2e1, t40 * t35, t40 ^ 2 / 0.2e1, -t13 * t40 + t15 * t35, t14 * t40 + t15 * t37, -t13 * t37 - t14 * t35, t14 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t34, t25 ^ 2 / 0.2e1, -t25 * t34, t34 ^ 2 / 0.2e1, t9 * t25 + t5 * t34, t9 * t27 - t6 * t34, -t6 * t25 - t5 * t27, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t24, t16 ^ 2 / 0.2e1, -t16 * t24, t24 ^ 2 / 0.2e1, t1 * t24 + t3 * t16, t3 * t18 - t2 * t24, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
