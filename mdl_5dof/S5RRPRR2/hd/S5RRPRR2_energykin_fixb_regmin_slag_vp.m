% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:23:34
% EndTime: 2021-01-15 21:23:34
% DurationCPUTime: 0.10s
% Computational Cost: add. (260->41), mult. (690->90), div. (0->0), fcn. (504->8), ass. (0->39)
t55 = qJD(1) ^ 2;
t66 = t55 / 0.2e1;
t65 = cos(qJ(4));
t54 = cos(qJ(2));
t64 = t54 * t55;
t63 = pkin(6) + qJ(3);
t52 = sin(qJ(2));
t61 = qJD(1) * t52;
t42 = qJD(2) * pkin(2) - t63 * t61;
t60 = qJD(1) * t54;
t43 = t63 * t60;
t48 = sin(pkin(9));
t49 = cos(pkin(9));
t33 = t49 * t42 - t48 * t43;
t40 = (t48 * t54 + t49 * t52) * qJD(1);
t28 = qJD(2) * pkin(3) - t40 * pkin(7) + t33;
t34 = t48 * t42 + t49 * t43;
t39 = t48 * t61 - t49 * t60;
t29 = -t39 * pkin(7) + t34;
t51 = sin(qJ(4));
t62 = t51 * t28 + t65 * t29;
t59 = qJD(1) * qJD(2);
t47 = qJD(2) + qJD(4);
t58 = t52 * t59;
t57 = t54 * t59;
t56 = t65 * t28 - t51 * t29;
t44 = qJD(3) + (-pkin(2) * t54 - pkin(1)) * qJD(1);
t35 = t39 * pkin(3) + t44;
t53 = cos(qJ(5));
t50 = sin(qJ(5));
t46 = qJD(5) + t47;
t32 = -t51 * t39 + t65 * t40;
t31 = t65 * t39 + t51 * t40;
t24 = t31 * pkin(4) + t35;
t23 = -t50 * t31 + t53 * t32;
t22 = t53 * t31 + t50 * t32;
t21 = -t31 * pkin(8) + t62;
t20 = t47 * pkin(4) - t32 * pkin(8) + t56;
t1 = [t66, 0, 0, t52 ^ 2 * t66, t52 * t64, t58, t57, qJD(2) ^ 2 / 0.2e1, pkin(1) * t64 - pkin(6) * t58, -t55 * pkin(1) * t52 - pkin(6) * t57, t33 * qJD(2) + t44 * t39, -t34 * qJD(2) + t44 * t40, -t33 * t40 - t34 * t39, t34 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1 + t44 ^ 2 / 0.2e1, t32 ^ 2 / 0.2e1, -t32 * t31, t32 * t47, -t31 * t47, t47 ^ 2 / 0.2e1, t35 * t31 + t56 * t47, t35 * t32 - t62 * t47, t23 ^ 2 / 0.2e1, -t23 * t22, t23 * t46, -t22 * t46, t46 ^ 2 / 0.2e1, t24 * t22 + (t53 * t20 - t50 * t21) * t46, t24 * t23 - (t50 * t20 + t53 * t21) * t46;];
T_reg = t1;
