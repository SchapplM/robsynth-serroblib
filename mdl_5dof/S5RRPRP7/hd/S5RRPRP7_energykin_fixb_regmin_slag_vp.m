% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:42
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRP7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP7_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:41:22
% EndTime: 2021-01-15 20:41:22
% DurationCPUTime: 0.09s
% Computational Cost: add. (244->38), mult. (595->83), div. (0->0), fcn. (384->6), ass. (0->34)
t52 = qJD(1) ^ 2;
t62 = t52 / 0.2e1;
t51 = cos(qJ(2));
t61 = t51 * t52;
t60 = pkin(6) + qJ(3);
t46 = sin(pkin(8));
t47 = cos(pkin(8));
t57 = qJD(1) * t51;
t49 = sin(qJ(2));
t58 = qJD(1) * t49;
t37 = t46 * t58 - t47 * t57;
t38 = (t46 * t51 + t47 * t49) * qJD(1);
t43 = qJD(3) + (-pkin(2) * t51 - pkin(1)) * qJD(1);
t27 = t37 * pkin(3) - t38 * pkin(7) + t43;
t41 = qJD(2) * pkin(2) - t60 * t58;
t42 = t60 * t57;
t32 = t46 * t41 + t47 * t42;
t30 = qJD(2) * pkin(7) + t32;
t48 = sin(qJ(4));
t50 = cos(qJ(4));
t59 = t48 * t27 + t50 * t30;
t56 = qJD(1) * qJD(2);
t55 = t49 * t56;
t54 = t51 * t56;
t31 = t47 * t41 - t46 * t42;
t53 = t50 * t27 - t48 * t30;
t29 = -qJD(2) * pkin(3) - t31;
t36 = qJD(4) + t37;
t34 = t48 * qJD(2) + t50 * t38;
t33 = -t50 * qJD(2) + t48 * t38;
t25 = t33 * pkin(4) - t34 * qJ(5) + t29;
t24 = t36 * qJ(5) + t59;
t23 = -t36 * pkin(4) + qJD(5) - t53;
t1 = [t62, 0, 0, t49 ^ 2 * t62, t49 * t61, t55, t54, qJD(2) ^ 2 / 0.2e1, pkin(1) * t61 - pkin(6) * t55, -t52 * pkin(1) * t49 - pkin(6) * t54, t31 * qJD(2) + t43 * t37, -t32 * qJD(2) + t43 * t38, -t31 * t38 - t32 * t37, t32 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1 + t43 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t33, t34 * t36, -t33 * t36, t36 ^ 2 / 0.2e1, t29 * t33 + t53 * t36, t29 * t34 - t59 * t36, -t23 * t36 + t25 * t33, t23 * t34 - t24 * t33, t24 * t36 - t25 * t34, t24 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1;];
T_reg = t1;
