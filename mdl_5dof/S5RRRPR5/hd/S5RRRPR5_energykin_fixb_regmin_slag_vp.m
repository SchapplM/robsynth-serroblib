% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:12
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR5_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:10:42
% EndTime: 2021-01-15 23:10:42
% DurationCPUTime: 0.10s
% Computational Cost: add. (298->41), mult. (733->90), div. (0->0), fcn. (514->8), ass. (0->39)
t63 = qJD(1) ^ 2;
t74 = t63 / 0.2e1;
t73 = pkin(7) + pkin(6);
t72 = cos(qJ(3));
t62 = cos(qJ(2));
t71 = t62 * t63;
t59 = sin(qJ(3));
t60 = sin(qJ(2));
t48 = (t59 * t62 + t72 * t60) * qJD(1);
t55 = qJD(2) + qJD(3);
t69 = qJD(1) * t60;
t50 = qJD(2) * pkin(2) - t73 * t69;
t68 = qJD(1) * t62;
t51 = t73 * t68;
t64 = t72 * t50 - t59 * t51;
t35 = t55 * pkin(3) - t48 * qJ(4) + t64;
t47 = t59 * t69 - t72 * t68;
t70 = t59 * t50 + t72 * t51;
t37 = -t47 * qJ(4) + t70;
t56 = sin(pkin(9));
t57 = cos(pkin(9));
t32 = t56 * t35 + t57 * t37;
t67 = qJD(1) * qJD(2);
t66 = t60 * t67;
t65 = t62 * t67;
t41 = t57 * t47 + t56 * t48;
t52 = (-pkin(2) * t62 - pkin(1)) * qJD(1);
t31 = t57 * t35 - t56 * t37;
t43 = t47 * pkin(3) + qJD(4) + t52;
t61 = cos(qJ(5));
t58 = sin(qJ(5));
t42 = -t56 * t47 + t57 * t48;
t40 = qJD(5) + t41;
t39 = t61 * t42 + t58 * t55;
t38 = t58 * t42 - t61 * t55;
t33 = t41 * pkin(4) - t42 * pkin(8) + t43;
t30 = t55 * pkin(8) + t32;
t29 = -t55 * pkin(4) - t31;
t1 = [t74, 0, 0, t60 ^ 2 * t74, t60 * t71, t66, t65, qJD(2) ^ 2 / 0.2e1, pkin(1) * t71 - pkin(6) * t66, -t63 * pkin(1) * t60 - pkin(6) * t65, t48 ^ 2 / 0.2e1, -t48 * t47, t48 * t55, -t47 * t55, t55 ^ 2 / 0.2e1, t52 * t47 + t64 * t55, t52 * t48 - t70 * t55, t31 * t55 + t43 * t41, -t32 * t55 + t43 * t42, -t31 * t42 - t32 * t41, t32 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1 + t43 ^ 2 / 0.2e1, t39 ^ 2 / 0.2e1, -t39 * t38, t39 * t40, -t38 * t40, t40 ^ 2 / 0.2e1, (-t58 * t30 + t61 * t33) * t40 + t29 * t38, -(t61 * t30 + t58 * t33) * t40 + t29 * t39;];
T_reg = t1;
