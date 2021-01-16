% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:38
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPP6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP6_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:37:05
% EndTime: 2021-01-15 22:37:05
% DurationCPUTime: 0.09s
% Computational Cost: add. (307->39), mult. (674->84), div. (0->0), fcn. (421->6), ass. (0->33)
t52 = qJD(1) ^ 2;
t62 = t52 / 0.2e1;
t61 = cos(qJ(3));
t51 = cos(qJ(2));
t60 = t51 * t52;
t49 = sin(qJ(3));
t50 = sin(qJ(2));
t58 = qJD(1) * t50;
t39 = t49 * qJD(2) + t61 * t58;
t57 = t51 * qJD(1);
t43 = -qJD(3) + t57;
t37 = (-pkin(2) * t51 - pkin(7) * t50 - pkin(1)) * qJD(1);
t42 = pkin(6) * t57 + qJD(2) * pkin(7);
t53 = t61 * t37 - t49 * t42;
t28 = -t43 * pkin(3) - t39 * qJ(4) + t53;
t38 = -t61 * qJD(2) + t49 * t58;
t59 = t49 * t37 + t61 * t42;
t30 = -t38 * qJ(4) + t59;
t47 = sin(pkin(8));
t48 = cos(pkin(8));
t25 = t47 * t28 + t48 * t30;
t56 = qJD(1) * qJD(2);
t55 = t50 * t56;
t54 = t51 * t56;
t41 = -qJD(2) * pkin(2) + pkin(6) * t58;
t24 = t48 * t28 - t47 * t30;
t33 = t38 * pkin(3) + qJD(4) + t41;
t32 = -t47 * t38 + t48 * t39;
t31 = t48 * t38 + t47 * t39;
t26 = t31 * pkin(4) - t32 * qJ(5) + t33;
t23 = -t43 * qJ(5) + t25;
t22 = t43 * pkin(4) + qJD(5) - t24;
t1 = [t62, 0, 0, t50 ^ 2 * t62, t50 * t60, t55, t54, qJD(2) ^ 2 / 0.2e1, pkin(1) * t60 - pkin(6) * t55, -t52 * pkin(1) * t50 - pkin(6) * t54, t39 ^ 2 / 0.2e1, -t39 * t38, -t39 * t43, t38 * t43, t43 ^ 2 / 0.2e1, t41 * t38 - t53 * t43, t41 * t39 + t59 * t43, -t24 * t43 + t33 * t31, t25 * t43 + t33 * t32, -t24 * t32 - t25 * t31, t25 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1, t22 * t43 + t26 * t31, t22 * t32 - t23 * t31, -t23 * t43 - t26 * t32, t23 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1;];
T_reg = t1;
