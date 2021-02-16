% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRPP4
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
% Datum: 2021-01-15 22:26
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:24:47
% EndTime: 2021-01-15 22:24:47
% DurationCPUTime: 0.09s
% Computational Cost: add. (307->38), mult. (738->83), div. (0->0), fcn. (485->6), ass. (0->34)
t47 = qJD(1) ^ 2;
t58 = t47 / 0.2e1;
t57 = pkin(7) + pkin(6);
t56 = cos(qJ(3));
t46 = cos(qJ(2));
t55 = t46 * t47;
t44 = sin(qJ(3));
t45 = sin(qJ(2));
t35 = (t44 * t46 + t56 * t45) * qJD(1);
t41 = qJD(2) + qJD(3);
t53 = qJD(1) * t45;
t37 = qJD(2) * pkin(2) - t57 * t53;
t52 = qJD(1) * t46;
t38 = t57 * t52;
t48 = t56 * t37 - t44 * t38;
t25 = t41 * pkin(3) - t35 * qJ(4) + t48;
t34 = t44 * t53 - t56 * t52;
t54 = t44 * t37 + t56 * t38;
t27 = -t34 * qJ(4) + t54;
t42 = sin(pkin(8));
t43 = cos(pkin(8));
t22 = t42 * t25 + t43 * t27;
t51 = qJD(1) * qJD(2);
t50 = t45 * t51;
t49 = t46 * t51;
t39 = (-pkin(2) * t46 - pkin(1)) * qJD(1);
t21 = t43 * t25 - t42 * t27;
t30 = t34 * pkin(3) + qJD(4) + t39;
t29 = -t42 * t34 + t43 * t35;
t28 = t43 * t34 + t42 * t35;
t23 = t28 * pkin(4) - t29 * qJ(5) + t30;
t20 = t41 * qJ(5) + t22;
t19 = -t41 * pkin(4) + qJD(5) - t21;
t1 = [t58, 0, 0, t45 ^ 2 * t58, t45 * t55, t50, t49, qJD(2) ^ 2 / 0.2e1, pkin(1) * t55 - pkin(6) * t50, -t47 * pkin(1) * t45 - pkin(6) * t49, t35 ^ 2 / 0.2e1, -t35 * t34, t35 * t41, -t34 * t41, t41 ^ 2 / 0.2e1, t39 * t34 + t48 * t41, t39 * t35 - t54 * t41, t21 * t41 + t30 * t28, -t22 * t41 + t30 * t29, -t21 * t29 - t22 * t28, t22 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1, -t19 * t41 + t23 * t28, t19 * t29 - t20 * t28, t20 * t41 - t23 * t29, t20 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1;];
T_reg = t1;
