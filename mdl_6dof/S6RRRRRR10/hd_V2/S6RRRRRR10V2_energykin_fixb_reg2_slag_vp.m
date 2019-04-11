% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRRR10V2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRR10V2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10V2_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:51:47
% EndTime: 2019-04-11 14:51:47
% DurationCPUTime: 0.19s
% Computational Cost: add. (673->54), mult. (1424->141), div. (0->0), fcn. (1075->10), ass. (0->45)
t46 = qJD(1) ^ 2;
t59 = t46 / 0.2e1;
t58 = cos(qJ(4));
t57 = cos(qJ(5));
t56 = cos(qJ(6));
t44 = cos(qJ(2));
t55 = t44 * t46;
t37 = qJD(2) + qJD(3);
t41 = sin(qJ(3));
t54 = pkin(2) * qJD(2);
t50 = t41 * t54;
t32 = t37 * pkin(5) + t50;
t40 = sin(qJ(4));
t42 = sin(qJ(2));
t43 = cos(qJ(3));
t27 = (t41 * t42 - t43 * t44) * qJD(1);
t29 = (t41 * t44 + t42 * t43) * qJD(1);
t34 = (-pkin(2) * t44 - pkin(1)) * qJD(1);
t48 = t27 * pkin(3) - t29 * pkin(5) + t34;
t14 = t58 * t32 + t40 * t48;
t49 = t43 * t54;
t33 = -t37 * pkin(3) - t49;
t39 = sin(qJ(5));
t10 = t57 * t14 + t39 * t33;
t8 = t39 * t14 - t57 * t33;
t53 = t8 ^ 2 / 0.2e1;
t12 = t40 * t32 - t58 * t48;
t52 = t12 ^ 2 / 0.2e1;
t51 = qJD(1) * qJD(2);
t24 = t58 * t29 + t40 * t37;
t26 = qJD(4) + t27;
t16 = t39 * t24 - t57 * t26;
t22 = t40 * t29 - t58 * t37;
t45 = qJD(2) ^ 2;
t38 = sin(qJ(6));
t21 = qJD(5) + t22;
t18 = t57 * t24 + t39 * t26;
t15 = qJD(6) + t16;
t7 = t56 * t18 + t38 * t21;
t5 = t38 * t18 - t56 * t21;
t4 = t21 * pkin(6) + t10;
t3 = -t18 * pkin(6) + t12;
t2 = t38 * t3 + t56 * t4;
t1 = t56 * t3 - t38 * t4;
t6 = [0, 0, 0, 0, 0, t59, 0, 0, 0, 0, t42 ^ 2 * t59, t42 * t55, t42 * t51, t44 ^ 2 * t59, t44 * t51, t45 / 0.2e1, pkin(1) * t55, -t46 * pkin(1) * t42, 0, pkin(1) ^ 2 * t59, t29 ^ 2 / 0.2e1, -t29 * t27, t29 * t37, t27 ^ 2 / 0.2e1, -t27 * t37, t37 ^ 2 / 0.2e1, t34 * t27 + t37 * t49, t34 * t29 - t37 * t50 (-t27 * t41 - t29 * t43) * t54, t34 ^ 2 / 0.2e1 + (t41 ^ 2 / 0.2e1 + t43 ^ 2 / 0.2e1) * pkin(2) ^ 2 * t45, t24 ^ 2 / 0.2e1, -t24 * t22, t24 * t26, t22 ^ 2 / 0.2e1, -t22 * t26, t26 ^ 2 / 0.2e1, -t12 * t26 + t33 * t22, -t14 * t26 + t33 * t24, t12 * t24 - t14 * t22, t14 ^ 2 / 0.2e1 + t52 + t33 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t21, t16 ^ 2 / 0.2e1, -t16 * t21, t21 ^ 2 / 0.2e1, t12 * t16 - t8 * t21, -t10 * t21 + t12 * t18, -t10 * t16 + t8 * t18, t10 ^ 2 / 0.2e1 + t53 + t52, t7 ^ 2 / 0.2e1, -t7 * t5, t7 * t15, t5 ^ 2 / 0.2e1, -t5 * t15, t15 ^ 2 / 0.2e1, t1 * t15 + t8 * t5, -t2 * t15 + t8 * t7, -t1 * t7 - t2 * t5, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t53;];
T_reg  = t6;
