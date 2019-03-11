% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRPRP8
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
% Datum: 2019-03-09 17:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRP8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP8_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:19:38
% EndTime: 2019-03-09 17:19:39
% DurationCPUTime: 0.16s
% Computational Cost: add. (546->62), mult. (1162->124), div. (0->0), fcn. (731->6), ass. (0->52)
t54 = qJD(1) ^ 2;
t68 = t54 / 0.2e1;
t52 = sin(qJ(2));
t53 = cos(qJ(2));
t28 = (-pkin(2) * t53 - pkin(8) * t52 - pkin(1)) * qJD(1);
t61 = t53 * qJD(1);
t39 = pkin(7) * t61 + qJD(2) * pkin(8);
t51 = sin(qJ(3));
t67 = cos(qJ(3));
t24 = t51 * t28 + t67 * t39;
t43 = -qJD(3) + t61;
t16 = -t43 * qJ(4) + t24;
t62 = qJD(1) * t52;
t32 = -t67 * qJD(2) + t51 * t62;
t11 = t32 * pkin(9) + t16;
t50 = sin(qJ(5));
t66 = cos(qJ(5));
t34 = t51 * qJD(2) + t67 * t62;
t23 = t67 * t28 - t51 * t39;
t56 = qJD(4) - t23;
t9 = -t34 * pkin(9) + (pkin(3) + pkin(4)) * t43 + t56;
t4 = t66 * t11 + t50 * t9;
t65 = t34 * t32;
t64 = t43 * t32;
t63 = t53 * t54;
t38 = -qJD(2) * pkin(2) + pkin(7) * t62;
t60 = t32 ^ 2 / 0.2e1;
t59 = qJD(1) * qJD(2);
t3 = -t50 * t11 + t66 * t9;
t58 = t52 * t59;
t57 = t53 * t59;
t17 = t32 * pkin(3) - t34 * qJ(4) + t38;
t12 = -t32 * pkin(4) - t17;
t48 = t53 ^ 2;
t47 = t52 ^ 2;
t42 = qJD(5) + t43;
t40 = t43 ^ 2 / 0.2e1;
t37 = t42 ^ 2 / 0.2e1;
t29 = t34 ^ 2 / 0.2e1;
t25 = t34 * t43;
t22 = t50 * t32 + t66 * t34;
t20 = -t66 * t32 + t50 * t34;
t19 = t22 ^ 2 / 0.2e1;
t18 = t20 ^ 2 / 0.2e1;
t15 = t43 * pkin(3) + t56;
t14 = t42 * t22;
t13 = t42 * t20;
t8 = t22 * t20;
t5 = t20 * pkin(5) + qJD(6) + t12;
t2 = -t20 * qJ(6) + t4;
t1 = t42 * pkin(5) - t22 * qJ(6) + t3;
t6 = [0, 0, 0, 0, 0, t68, 0, 0, 0, 0, t47 * t68, t52 * t63, t58, t48 * t68, t57, qJD(2) ^ 2 / 0.2e1, pkin(1) * t63 - pkin(7) * t58, -t54 * pkin(1) * t52 - pkin(7) * t57 (t47 + t48) * t54 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t48 / 0.2e1 + t47 / 0.2e1) * pkin(7) ^ 2) * t54, t29, -t65, -t25, t60, t64, t40, -t23 * t43 + t38 * t32, t24 * t43 + t38 * t34, -t23 * t34 - t24 * t32, t24 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1, t29, -t25, t65, t40, -t64, t60, t15 * t43 + t17 * t32, t15 * t34 - t16 * t32, -t16 * t43 - t17 * t34, t16 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t19, -t8, t14, t18, -t13, t37, t12 * t20 + t3 * t42, t12 * t22 - t4 * t42, -t4 * t20 - t3 * t22, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1, t19, -t8, t14, t18, -t13, t37, t1 * t42 + t5 * t20, -t2 * t42 + t5 * t22, -t1 * t22 - t2 * t20, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t6;
