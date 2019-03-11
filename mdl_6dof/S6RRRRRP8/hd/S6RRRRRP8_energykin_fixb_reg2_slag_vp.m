% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRRP8
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
% Datum: 2019-03-10 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRP8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP8_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:55:48
% EndTime: 2019-03-10 01:55:48
% DurationCPUTime: 0.23s
% Computational Cost: add. (1289->66), mult. (2994->150), div. (0->0), fcn. (2344->10), ass. (0->55)
t68 = cos(pkin(6)) * qJD(1);
t52 = qJD(2) + t68;
t58 = sin(qJ(3));
t59 = sin(qJ(2));
t54 = sin(pkin(6));
t69 = qJD(1) * t54;
t64 = t59 * t69;
t75 = cos(qJ(3));
t38 = -t75 * t52 + t58 * t64;
t40 = t58 * t52 + t75 * t64;
t57 = sin(qJ(4));
t74 = cos(qJ(4));
t27 = t74 * t38 + t57 * t40;
t29 = -t57 * t38 + t74 * t40;
t60 = cos(qJ(2));
t65 = pkin(1) * t68;
t41 = -pkin(8) * t64 + t60 * t65;
t34 = -t52 * pkin(2) - t41;
t30 = t38 * pkin(3) + t34;
t12 = t27 * pkin(4) - t29 * pkin(11) + t30;
t56 = sin(qJ(5));
t73 = cos(qJ(5));
t63 = t60 * t69;
t42 = pkin(8) * t63 + t59 * t65;
t35 = t52 * pkin(9) + t42;
t37 = (-pkin(2) * t60 - pkin(9) * t59 - pkin(1)) * t69;
t23 = -t58 * t35 + t75 * t37;
t47 = -qJD(3) + t63;
t15 = -t47 * pkin(3) - t40 * pkin(10) + t23;
t24 = t75 * t35 + t58 * t37;
t18 = -t38 * pkin(10) + t24;
t10 = t57 * t15 + t74 * t18;
t44 = -qJD(4) + t47;
t8 = -t44 * pkin(11) + t10;
t4 = t56 * t12 + t73 * t8;
t20 = t56 * t29 + t73 * t44;
t22 = t73 * t29 - t56 * t44;
t72 = t22 * t20;
t26 = qJD(5) + t27;
t71 = t26 * t20;
t61 = qJD(1) ^ 2;
t70 = t54 ^ 2 * t61;
t67 = t20 ^ 2 / 0.2e1;
t66 = t60 * t70;
t62 = t70 / 0.2e1;
t9 = t74 * t15 - t57 * t18;
t3 = t73 * t12 - t56 * t8;
t7 = t44 * pkin(4) - t9;
t25 = t26 ^ 2 / 0.2e1;
t19 = t22 ^ 2 / 0.2e1;
t13 = t22 * t26;
t5 = t20 * pkin(5) - t22 * qJ(6) + t7;
t2 = t26 * qJ(6) + t4;
t1 = -t26 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t61 / 0.2e1, 0, 0, 0, 0, t59 ^ 2 * t62, t59 * t66, t52 * t64, t60 ^ 2 * t62, t52 * t63, t52 ^ 2 / 0.2e1, pkin(1) * t66 + t41 * t52, -pkin(1) * t59 * t70 - t42 * t52 (-t41 * t59 + t42 * t60) * t69, t42 ^ 2 / 0.2e1 + t41 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t62, t40 ^ 2 / 0.2e1, -t40 * t38, -t40 * t47, t38 ^ 2 / 0.2e1, t38 * t47, t47 ^ 2 / 0.2e1, -t23 * t47 + t34 * t38, t24 * t47 + t34 * t40, -t23 * t40 - t24 * t38, t24 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1, t29 ^ 2 / 0.2e1, -t29 * t27, -t29 * t44, t27 ^ 2 / 0.2e1, t27 * t44, t44 ^ 2 / 0.2e1, t30 * t27 - t9 * t44, t10 * t44 + t30 * t29, -t10 * t27 - t9 * t29, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1, t19, -t72, t13, t67, -t71, t25, t7 * t20 + t3 * t26, t7 * t22 - t4 * t26, -t4 * t20 - t3 * t22, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t19, t13, t72, t25, t71, t67, -t1 * t26 + t5 * t20, t1 * t22 - t2 * t20, t2 * t26 - t5 * t22, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
