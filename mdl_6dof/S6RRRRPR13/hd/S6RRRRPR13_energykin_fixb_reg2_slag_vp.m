% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPR13_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR13_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:04:02
% EndTime: 2019-03-10 00:04:03
% DurationCPUTime: 0.24s
% Computational Cost: add. (1076->67), mult. (2474->150), div. (0->0), fcn. (1896->10), ass. (0->58)
t78 = -pkin(4) - pkin(5);
t77 = cos(qJ(3));
t76 = cos(qJ(4));
t75 = cos(qJ(6));
t70 = cos(pkin(6)) * qJD(1);
t52 = qJD(2) + t70;
t58 = sin(qJ(3));
t59 = sin(qJ(2));
t54 = sin(pkin(6));
t71 = qJD(1) * t54;
t66 = t59 * t71;
t42 = t58 * t52 + t77 * t66;
t60 = cos(qJ(2));
t65 = t60 * t71;
t47 = -qJD(3) + t65;
t57 = sin(qJ(4));
t25 = t57 * t42 + t76 * t47;
t27 = t76 * t42 - t57 * t47;
t74 = t27 * t25;
t40 = -t77 * t52 + t58 * t66;
t39 = qJD(4) + t40;
t73 = t39 * t25;
t61 = qJD(1) ^ 2;
t72 = t54 ^ 2 * t61;
t67 = pkin(1) * t70;
t43 = -pkin(8) * t66 + t60 * t67;
t32 = -t52 * pkin(2) - t43;
t15 = t40 * pkin(3) - t42 * pkin(10) + t32;
t44 = pkin(8) * t65 + t59 * t67;
t33 = t52 * pkin(9) + t44;
t38 = (-pkin(2) * t60 - pkin(9) * t59 - pkin(1)) * t71;
t22 = t77 * t33 + t58 * t38;
t19 = -t47 * pkin(10) + t22;
t10 = t57 * t15 + t76 * t19;
t21 = -t58 * t33 + t77 * t38;
t69 = t25 ^ 2 / 0.2e1;
t68 = t60 * t72;
t7 = t39 * qJ(5) + t10;
t64 = t72 / 0.2e1;
t18 = t47 * pkin(3) - t21;
t9 = t76 * t15 - t57 * t19;
t63 = qJD(5) - t9;
t62 = t27 * qJ(5) - t18;
t56 = sin(qJ(6));
t36 = -qJD(6) + t39;
t34 = t39 ^ 2 / 0.2e1;
t24 = t27 ^ 2 / 0.2e1;
t20 = t27 * t39;
t13 = t56 * t25 + t75 * t27;
t11 = -t75 * t25 + t56 * t27;
t8 = t25 * pkin(4) - t62;
t6 = -t39 * pkin(4) + t63;
t5 = t78 * t25 + t62;
t4 = t25 * pkin(11) + t7;
t3 = -t27 * pkin(11) + t78 * t39 + t63;
t2 = t56 * t3 + t75 * t4;
t1 = t75 * t3 - t56 * t4;
t12 = [0, 0, 0, 0, 0, t61 / 0.2e1, 0, 0, 0, 0, t59 ^ 2 * t64, t59 * t68, t52 * t66, t60 ^ 2 * t64, t52 * t65, t52 ^ 2 / 0.2e1, pkin(1) * t68 + t43 * t52, -pkin(1) * t59 * t72 - t44 * t52 (-t43 * t59 + t44 * t60) * t71, t44 ^ 2 / 0.2e1 + t43 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t64, t42 ^ 2 / 0.2e1, -t42 * t40, -t42 * t47, t40 ^ 2 / 0.2e1, t40 * t47, t47 ^ 2 / 0.2e1, -t21 * t47 + t32 * t40, t22 * t47 + t32 * t42, -t21 * t42 - t22 * t40, t22 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1, t24, -t74, t20, t69, -t73, t34, t18 * t25 + t9 * t39, -t10 * t39 + t18 * t27, -t10 * t25 - t9 * t27, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t24, t20, t74, t34, t73, t69, t8 * t25 - t6 * t39, -t7 * t25 + t6 * t27, -t8 * t27 + t7 * t39, t7 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t13 * t11, -t13 * t36, t11 ^ 2 / 0.2e1, t11 * t36, t36 ^ 2 / 0.2e1, -t1 * t36 + t5 * t11, t5 * t13 + t2 * t36, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t12;
