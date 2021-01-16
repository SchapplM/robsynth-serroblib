% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRPR1
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
% Datum: 2021-01-15 22:50
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:48:56
% EndTime: 2021-01-15 22:48:56
% DurationCPUTime: 0.10s
% Computational Cost: add. (306->41), mult. (773->90), div. (0->0), fcn. (554->8), ass. (0->39)
t54 = qJD(1) ^ 2;
t66 = t54 / 0.2e1;
t65 = pkin(7) + pkin(6);
t64 = cos(qJ(3));
t53 = cos(qJ(2));
t63 = t53 * t54;
t50 = sin(qJ(3));
t51 = sin(qJ(2));
t40 = (t50 * t53 + t64 * t51) * qJD(1);
t47 = qJD(2) + qJD(3);
t60 = qJD(1) * t51;
t42 = qJD(2) * pkin(2) - t65 * t60;
t59 = qJD(1) * t53;
t43 = t65 * t59;
t55 = t64 * t42 - t50 * t43;
t29 = t47 * pkin(3) - t40 * qJ(4) + t55;
t39 = t50 * t60 - t64 * t59;
t62 = t50 * t42 + t64 * t43;
t31 = -t39 * qJ(4) + t62;
t48 = sin(pkin(9));
t61 = cos(pkin(9));
t23 = t48 * t29 + t61 * t31;
t58 = qJD(1) * qJD(2);
t57 = t51 * t58;
t56 = t53 * t58;
t22 = t61 * t29 - t48 * t31;
t44 = (-pkin(2) * t53 - pkin(1)) * qJD(1);
t35 = t39 * pkin(3) + qJD(4) + t44;
t52 = cos(qJ(5));
t49 = sin(qJ(5));
t46 = qJD(5) + t47;
t34 = -t48 * t39 + t61 * t40;
t33 = t61 * t39 + t48 * t40;
t26 = t33 * pkin(4) + t35;
t25 = -t49 * t33 + t52 * t34;
t24 = t52 * t33 + t49 * t34;
t21 = -t33 * pkin(8) + t23;
t20 = t47 * pkin(4) - t34 * pkin(8) + t22;
t1 = [t66, 0, 0, t51 ^ 2 * t66, t51 * t63, t57, t56, qJD(2) ^ 2 / 0.2e1, pkin(1) * t63 - pkin(6) * t57, -t54 * pkin(1) * t51 - pkin(6) * t56, t40 ^ 2 / 0.2e1, -t40 * t39, t40 * t47, -t39 * t47, t47 ^ 2 / 0.2e1, t44 * t39 + t55 * t47, t44 * t40 - t62 * t47, t22 * t47 + t35 * t33, -t23 * t47 + t35 * t34, -t22 * t34 - t23 * t33, t23 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1 + t35 ^ 2 / 0.2e1, t25 ^ 2 / 0.2e1, -t25 * t24, t25 * t46, -t24 * t46, t46 ^ 2 / 0.2e1, t26 * t24 + (t52 * t20 - t49 * t21) * t46, t26 * t25 - (t49 * t20 + t52 * t21) * t46;];
T_reg = t1;
