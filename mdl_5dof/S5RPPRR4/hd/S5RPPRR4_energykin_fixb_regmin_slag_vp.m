% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:16:57
% EndTime: 2022-01-23 09:16:58
% DurationCPUTime: 0.08s
% Computational Cost: add. (174->42), mult. (496->92), div. (0->0), fcn. (340->8), ass. (0->38)
t54 = sin(pkin(8));
t55 = cos(pkin(9));
t70 = t54 * t55;
t56 = cos(pkin(8));
t40 = qJD(2) + (-pkin(2) * t56 - qJ(3) * t54 - pkin(1)) * qJD(1);
t39 = t55 * t40;
t53 = sin(pkin(9));
t29 = t39 + (-pkin(6) * t70 + (-qJ(2) * t53 - pkin(3)) * t56) * qJD(1);
t66 = qJ(2) * qJD(1);
t64 = t56 * t66;
t34 = t53 * t40 + t55 * t64;
t68 = qJD(1) * t54;
t65 = t53 * t68;
t32 = -pkin(6) * t65 + t34;
t58 = sin(qJ(4));
t60 = cos(qJ(4));
t69 = t58 * t29 + t60 * t32;
t67 = t56 * qJD(1);
t46 = t54 * t66 + qJD(3);
t41 = pkin(3) * t65 + t46;
t63 = t60 * t29 - t58 * t32;
t47 = -qJD(4) + t67;
t61 = qJD(1) ^ 2;
t59 = cos(qJ(5));
t57 = sin(qJ(5));
t52 = t56 ^ 2;
t51 = t54 ^ 2;
t50 = -qJD(1) * pkin(1) + qJD(2);
t44 = -qJD(5) + t47;
t37 = (-t53 * t58 + t55 * t60) * t68;
t36 = (t53 * t60 + t55 * t58) * t68;
t33 = -t53 * t64 + t39;
t31 = t36 * pkin(4) + t41;
t26 = -t57 * t36 + t59 * t37;
t25 = t59 * t36 + t57 * t37;
t24 = -t36 * pkin(7) + t69;
t23 = -t47 * pkin(4) - t37 * pkin(7) + t63;
t1 = [t61 / 0.2e1, 0, 0, -t50 * t67, (t51 + t52) * t61 * qJ(2), t50 ^ 2 / 0.2e1 + (t52 / 0.2e1 + t51 / 0.2e1) * qJ(2) ^ 2 * t61, (t46 * t53 * t54 - t33 * t56) * qJD(1), (t34 * t56 + t46 * t70) * qJD(1), t34 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1 + t46 ^ 2 / 0.2e1, t37 ^ 2 / 0.2e1, -t37 * t36, -t37 * t47, t36 * t47, t47 ^ 2 / 0.2e1, t41 * t36 - t63 * t47, t41 * t37 + t69 * t47, t26 ^ 2 / 0.2e1, -t26 * t25, -t26 * t44, t25 * t44, t44 ^ 2 / 0.2e1, -(t59 * t23 - t57 * t24) * t44 + t31 * t25, (t57 * t23 + t59 * t24) * t44 + t31 * t26;];
T_reg = t1;
