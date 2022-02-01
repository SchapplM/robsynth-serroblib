% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:25:59
% EndTime: 2022-01-23 09:25:59
% DurationCPUTime: 0.20s
% Computational Cost: add. (225->43), mult. (589->97), div. (0->0), fcn. (389->8), ass. (0->38)
t49 = sin(pkin(8));
t46 = t49 ^ 2;
t56 = qJD(1) ^ 2;
t66 = t46 * t56;
t51 = cos(pkin(8));
t36 = qJD(2) + (-pkin(2) * t51 - pkin(6) * t49 - pkin(1)) * qJD(1);
t55 = cos(qJ(3));
t35 = t55 * t36;
t62 = t51 * qJD(1);
t42 = -qJD(3) + t62;
t53 = sin(qJ(3));
t27 = -t42 * pkin(3) + t35 + (-qJ(2) * t51 * t53 - qJ(4) * t49 * t55) * qJD(1);
t63 = qJD(1) * t49;
t59 = t53 * t63;
t61 = qJ(2) * qJD(1);
t58 = t51 * t61;
t65 = t53 * t36 + t55 * t58;
t30 = -qJ(4) * t59 + t65;
t48 = sin(pkin(9));
t50 = cos(pkin(9));
t22 = t48 * t27 + t50 * t30;
t64 = qJ(2) * t56;
t60 = t46 * t64;
t37 = pkin(3) * t59 + t49 * t61 + qJD(4);
t21 = t50 * t27 - t48 * t30;
t54 = cos(qJ(5));
t52 = sin(qJ(5));
t47 = t51 ^ 2;
t45 = -qJD(1) * pkin(1) + qJD(2);
t39 = -qJD(5) + t42;
t33 = (-t48 * t53 + t50 * t55) * t63;
t32 = (t48 * t55 + t50 * t53) * t63;
t29 = t32 * pkin(4) + t37;
t24 = -t52 * t32 + t54 * t33;
t23 = t54 * t32 + t52 * t33;
t20 = -t32 * pkin(7) + t22;
t19 = -t42 * pkin(4) - t33 * pkin(7) + t21;
t1 = [t56 / 0.2e1, 0, 0, -t45 * t62, (t46 + t47) * t64, t45 ^ 2 / 0.2e1 + (t47 / 0.2e1 + t46 / 0.2e1) * qJ(2) ^ 2 * t56, t55 ^ 2 * t66 / 0.2e1, -t55 * t53 * t66, -t55 * t42 * t63, t42 * t59, t42 ^ 2 / 0.2e1, t53 * t60 - (-t53 * t58 + t35) * t42, t65 * t42 + t55 * t60, -t21 * t42 + t37 * t32, t22 * t42 + t37 * t33, -t21 * t33 - t22 * t32, t22 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1 + t37 ^ 2 / 0.2e1, t24 ^ 2 / 0.2e1, -t24 * t23, -t24 * t39, t23 * t39, t39 ^ 2 / 0.2e1, -(t54 * t19 - t52 * t20) * t39 + t29 * t23, (t52 * t19 + t54 * t20) * t39 + t29 * t24;];
T_reg = t1;
