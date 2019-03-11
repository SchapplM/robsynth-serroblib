% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRP4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP4_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:45:55
% EndTime: 2019-03-09 16:45:55
% DurationCPUTime: 0.20s
% Computational Cost: add. (537->59), mult. (1195->123), div. (0->0), fcn. (776->6), ass. (0->55)
t41 = sin(qJ(3));
t42 = sin(qJ(2));
t43 = cos(qJ(3));
t44 = cos(qJ(2));
t27 = (t41 * t44 + t42 * t43) * qJD(1);
t45 = qJD(1) ^ 2;
t67 = t45 / 0.2e1;
t66 = pkin(3) + pkin(9);
t65 = -pkin(8) - pkin(7);
t40 = sin(qJ(5));
t64 = cos(qJ(5));
t37 = qJD(2) + qJD(3);
t57 = qJD(1) * t42;
t31 = qJD(2) * pkin(2) + t65 * t57;
t56 = qJD(1) * t44;
t32 = t65 * t56;
t15 = t43 * t31 + t41 * t32;
t49 = qJD(4) - t15;
t8 = t27 * pkin(4) - t66 * t37 + t49;
t25 = t41 * t57 - t43 * t56;
t33 = (-pkin(2) * t44 - pkin(1)) * qJD(1);
t47 = -t27 * qJ(4) + t33;
t9 = t66 * t25 + t47;
t4 = t40 * t8 + t64 * t9;
t18 = -t64 * t25 + t40 * t37;
t20 = t40 * t25 + t64 * t37;
t63 = t20 * t18;
t24 = qJD(5) + t27;
t62 = t24 * t18;
t61 = t27 * t25;
t60 = t27 * t37;
t59 = t37 * t25;
t58 = t44 * t45;
t16 = t41 * t31 - t43 * t32;
t55 = t18 ^ 2 / 0.2e1;
t54 = t25 ^ 2 / 0.2e1;
t53 = t27 ^ 2 / 0.2e1;
t52 = qJD(1) * qJD(2);
t14 = -t37 * qJ(4) - t16;
t51 = t42 * t52;
t50 = t44 * t52;
t10 = -t25 * pkin(4) - t14;
t3 = -t40 * t9 + t64 * t8;
t39 = t44 ^ 2;
t38 = t42 ^ 2;
t35 = t37 ^ 2 / 0.2e1;
t21 = t24 ^ 2 / 0.2e1;
t17 = t20 ^ 2 / 0.2e1;
t13 = -t37 * pkin(3) + t49;
t12 = t25 * pkin(3) + t47;
t11 = t20 * t24;
t5 = t18 * pkin(5) - t20 * qJ(6) + t10;
t2 = t24 * qJ(6) + t4;
t1 = -t24 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t67, 0, 0, 0, 0, t38 * t67, t42 * t58, t51, t39 * t67, t50, qJD(2) ^ 2 / 0.2e1, pkin(1) * t58 - pkin(7) * t51, -t45 * pkin(1) * t42 - pkin(7) * t50 (t38 + t39) * t45 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t39 / 0.2e1 + t38 / 0.2e1) * pkin(7) ^ 2) * t45, t53, -t61, t60, t54, -t59, t35, t15 * t37 + t33 * t25, -t16 * t37 + t33 * t27, -t15 * t27 - t16 * t25, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1, t35, -t60, t59, t53, -t61, t54, t13 * t27 + t14 * t25, -t12 * t25 + t13 * t37, -t12 * t27 - t14 * t37, t12 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1, t17, -t63, t11, t55, -t62, t21, t10 * t18 + t3 * t24, t10 * t20 - t4 * t24, -t4 * t18 - t3 * t20, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1, t17, t11, t63, t21, t62, t55, -t1 * t24 + t5 * t18, t1 * t20 - t2 * t18, t2 * t24 - t5 * t20, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
