% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRP5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP5_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:23:02
% EndTime: 2019-03-10 01:23:02
% DurationCPUTime: 0.21s
% Computational Cost: add. (1009->65), mult. (2216->144), div. (0->0), fcn. (1584->8), ass. (0->52)
t61 = qJD(1) ^ 2;
t71 = t61 / 0.2e1;
t58 = sin(qJ(2));
t60 = cos(qJ(2));
t37 = (-pkin(2) * t60 - pkin(8) * t58 - pkin(1)) * qJD(1);
t66 = t60 * qJD(1);
t46 = pkin(7) * t66 + qJD(2) * pkin(8);
t57 = sin(qJ(3));
t70 = cos(qJ(3));
t31 = t70 * t37 - t57 * t46;
t67 = qJD(1) * t58;
t40 = t57 * qJD(2) + t70 * t67;
t49 = -qJD(3) + t66;
t24 = -t49 * pkin(3) - t40 * pkin(9) + t31;
t32 = t57 * t37 + t70 * t46;
t38 = -t70 * qJD(2) + t57 * t67;
t26 = -t38 * pkin(9) + t32;
t56 = sin(qJ(4));
t59 = cos(qJ(4));
t13 = t56 * t24 + t59 * t26;
t28 = t59 * t38 + t56 * t40;
t10 = -t28 * pkin(10) + t13;
t55 = sin(qJ(5));
t69 = cos(qJ(5));
t12 = t59 * t24 - t56 * t26;
t30 = -t56 * t38 + t59 * t40;
t47 = -qJD(4) + t49;
t7 = -t47 * pkin(4) - t30 * pkin(10) + t12;
t4 = t69 * t10 + t55 * t7;
t68 = t60 * t61;
t65 = qJD(1) * qJD(2);
t3 = -t55 * t10 + t69 * t7;
t64 = t58 * t65;
t63 = t60 * t65;
t45 = -qJD(2) * pkin(2) + pkin(7) * t67;
t33 = t38 * pkin(3) + t45;
t21 = t28 * pkin(4) + t33;
t54 = t60 ^ 2;
t53 = t58 ^ 2;
t43 = -qJD(5) + t47;
t42 = t43 ^ 2 / 0.2e1;
t20 = -t55 * t28 + t69 * t30;
t18 = t69 * t28 + t55 * t30;
t17 = t20 ^ 2 / 0.2e1;
t16 = t18 ^ 2 / 0.2e1;
t15 = t20 * t43;
t14 = t18 * t43;
t11 = t18 * pkin(5) + qJD(6) + t21;
t8 = t20 * t18;
t2 = -t18 * qJ(6) + t4;
t1 = -t43 * pkin(5) - t20 * qJ(6) + t3;
t5 = [0, 0, 0, 0, 0, t71, 0, 0, 0, 0, t53 * t71, t58 * t68, t64, t54 * t71, t63, qJD(2) ^ 2 / 0.2e1, pkin(1) * t68 - pkin(7) * t64, -t61 * pkin(1) * t58 - pkin(7) * t63 (t53 + t54) * t61 * pkin(7) (pkin(1) ^ 2 / 0.2e1 + (t54 / 0.2e1 + t53 / 0.2e1) * pkin(7) ^ 2) * t61, t40 ^ 2 / 0.2e1, -t40 * t38, -t40 * t49, t38 ^ 2 / 0.2e1, t38 * t49, t49 ^ 2 / 0.2e1, -t31 * t49 + t45 * t38, t32 * t49 + t45 * t40, -t31 * t40 - t32 * t38, t32 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1 + t45 ^ 2 / 0.2e1, t30 ^ 2 / 0.2e1, -t30 * t28, -t30 * t47, t28 ^ 2 / 0.2e1, t28 * t47, t47 ^ 2 / 0.2e1, -t12 * t47 + t33 * t28, t13 * t47 + t33 * t30, -t12 * t30 - t13 * t28, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1, t17, -t8, -t15, t16, t14, t42, t21 * t18 - t3 * t43, t21 * t20 + t4 * t43, -t4 * t18 - t3 * t20, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t17, -t8, -t15, t16, t14, t42, -t1 * t43 + t11 * t18, t11 * t20 + t2 * t43, -t1 * t20 - t2 * t18, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1;];
T_reg  = t5;
