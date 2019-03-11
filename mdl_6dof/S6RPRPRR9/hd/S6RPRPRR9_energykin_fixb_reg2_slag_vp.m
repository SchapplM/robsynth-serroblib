% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRR9_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR9_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_energykin_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:05:04
% EndTime: 2019-03-09 04:05:04
% DurationCPUTime: 0.34s
% Computational Cost: add. (2118->82), mult. (7222->194), div. (0->0), fcn. (6038->14), ass. (0->68)
t73 = qJD(1) ^ 2;
t92 = t73 / 0.2e1;
t64 = sin(pkin(12));
t68 = cos(pkin(12));
t66 = sin(pkin(6));
t83 = qJD(1) * t66;
t76 = qJ(2) * t83;
t69 = cos(pkin(6));
t82 = qJD(1) * t69;
t80 = pkin(1) * t82;
t54 = t64 * t80 + t68 * t76;
t65 = sin(pkin(7));
t84 = cos(pkin(7));
t85 = t66 * t68;
t39 = (t65 * t69 + t84 * t85) * qJD(1) * pkin(9) + t54;
t59 = t68 * t80;
t87 = t64 * t66;
t41 = t59 + (pkin(2) * t69 + (-t84 * pkin(9) - qJ(2)) * t87) * qJD(1);
t48 = qJD(2) + (-pkin(9) * t64 * t65 - pkin(2) * t68 - pkin(1)) * t83;
t72 = sin(qJ(3));
t91 = cos(qJ(3));
t74 = t84 * t91;
t77 = t65 * t91;
t23 = -t72 * t39 + t41 * t74 + t48 * t77;
t75 = t72 * t84;
t86 = t65 * t72;
t44 = (t69 * t86 + (t91 * t64 + t68 * t75) * t66) * qJD(1);
t79 = t68 * t83;
t51 = t65 * t79 - t84 * t82 - qJD(3);
t19 = -t51 * pkin(3) - t44 * qJ(4) + t23;
t24 = t91 * t39 + t41 * t75 + t48 * t86;
t42 = t72 * t64 * t83 - t74 * t79 - t77 * t82;
t22 = -t42 * qJ(4) + t24;
t63 = sin(pkin(13));
t67 = cos(pkin(13));
t12 = t63 * t19 + t67 * t22;
t10 = -t51 * pkin(10) + t12;
t35 = -t65 * t41 + t84 * t48;
t25 = t42 * pkin(3) + qJD(4) + t35;
t32 = t67 * t42 + t63 * t44;
t34 = -t63 * t42 + t67 * t44;
t14 = t32 * pkin(4) - t34 * pkin(10) + t25;
t71 = sin(qJ(5));
t90 = cos(qJ(5));
t6 = t90 * t10 + t71 * t14;
t89 = cos(qJ(6));
t88 = t66 ^ 2 * t73;
t81 = t66 * t69 * t73;
t78 = t88 / 0.2e1;
t11 = t67 * t19 - t63 * t22;
t27 = t71 * t34 + t90 * t51;
t9 = t51 * pkin(4) - t11;
t5 = -t71 * t10 + t90 * t14;
t70 = sin(qJ(6));
t60 = -pkin(1) * t83 + qJD(2);
t53 = -t64 * t76 + t59;
t50 = t51 ^ 2 / 0.2e1;
t31 = qJD(5) + t32;
t29 = t90 * t34 - t71 * t51;
t26 = qJD(6) + t27;
t17 = t89 * t29 + t70 * t31;
t15 = t70 * t29 - t89 * t31;
t7 = t27 * pkin(5) - t29 * pkin(11) + t9;
t4 = t31 * pkin(11) + t6;
t3 = -t31 * pkin(5) - t5;
t2 = t89 * t4 + t70 * t7;
t1 = -t70 * t4 + t89 * t7;
t8 = [0, 0, 0, 0, 0, t92, 0, 0, 0, 0, t64 ^ 2 * t78, t64 * t68 * t88, t64 * t81, t68 ^ 2 * t78, t68 * t81, t69 ^ 2 * t92 (t53 * t69 - t60 * t85) * qJD(1) (-t54 * t69 + t60 * t87) * qJD(1) (-t53 * t64 + t54 * t68) * t83, t54 ^ 2 / 0.2e1 + t53 ^ 2 / 0.2e1 + t60 ^ 2 / 0.2e1, t44 ^ 2 / 0.2e1, -t44 * t42, -t44 * t51, t42 ^ 2 / 0.2e1, t42 * t51, t50, -t23 * t51 + t35 * t42, t24 * t51 + t35 * t44, -t23 * t44 - t24 * t42, t24 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + t35 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, -t34 * t51, t32 ^ 2 / 0.2e1, t32 * t51, t50, -t11 * t51 + t25 * t32, t12 * t51 + t25 * t34, -t11 * t34 - t12 * t32, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1, t29 ^ 2 / 0.2e1, -t29 * t27, t29 * t31, t27 ^ 2 / 0.2e1, -t27 * t31, t31 ^ 2 / 0.2e1, t9 * t27 + t5 * t31, t9 * t29 - t6 * t31, -t6 * t27 - t5 * t29, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t17 ^ 2 / 0.2e1, -t17 * t15, t17 * t26, t15 ^ 2 / 0.2e1, -t15 * t26, t26 ^ 2 / 0.2e1, t1 * t26 + t3 * t15, t3 * t17 - t2 * t26, -t1 * t17 - t2 * t15, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
