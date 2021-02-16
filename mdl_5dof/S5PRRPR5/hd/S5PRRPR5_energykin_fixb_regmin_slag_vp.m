% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:05
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRPR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR5_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:04:41
% EndTime: 2021-01-15 16:04:41
% DurationCPUTime: 0.08s
% Computational Cost: add. (153->36), mult. (391->82), div. (0->0), fcn. (278->10), ass. (0->36)
t54 = qJD(2) ^ 2;
t64 = t54 / 0.2e1;
t50 = sin(qJ(2));
t62 = qJD(1) * sin(pkin(5));
t38 = qJD(2) * pkin(7) + t50 * t62;
t52 = cos(qJ(3));
t61 = qJD(1) * cos(pkin(5));
t42 = t52 * t61;
t49 = sin(qJ(3));
t58 = qJ(4) * qJD(2);
t29 = qJD(3) * pkin(3) + t42 + (-t38 - t58) * t49;
t63 = t52 * t38 + t49 * t61;
t30 = t52 * t58 + t63;
t44 = sin(pkin(10));
t46 = cos(pkin(10));
t25 = t44 * t29 + t46 * t30;
t60 = qJD(2) * t49;
t59 = qJD(2) * t52;
t57 = qJD(2) * qJD(3);
t53 = cos(qJ(2));
t56 = t53 * t62;
t55 = qJD(2) * t62;
t35 = t44 * t60 - t46 * t59;
t24 = t46 * t29 - t44 * t30;
t33 = -t56 + qJD(4) + (-pkin(3) * t52 - pkin(2)) * qJD(2);
t51 = cos(qJ(5));
t48 = sin(qJ(5));
t39 = -qJD(2) * pkin(2) - t56;
t36 = (t44 * t52 + t46 * t49) * qJD(2);
t34 = qJD(5) + t35;
t32 = t48 * qJD(3) + t51 * t36;
t31 = -t51 * qJD(3) + t48 * t36;
t26 = t35 * pkin(4) - t36 * pkin(8) + t33;
t23 = qJD(3) * pkin(8) + t25;
t22 = -qJD(3) * pkin(4) - t24;
t1 = [qJD(1) ^ 2 / 0.2e1, t64, t53 * t55, -t50 * t55, t49 ^ 2 * t64, t49 * t54 * t52, t49 * t57, t52 * t57, qJD(3) ^ 2 / 0.2e1, (-t49 * t38 + t42) * qJD(3) - t39 * t59, -t63 * qJD(3) + t39 * t60, t24 * qJD(3) + t33 * t35, -t25 * qJD(3) + t33 * t36, -t24 * t36 - t25 * t35, t25 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1, t32 ^ 2 / 0.2e1, -t32 * t31, t32 * t34, -t31 * t34, t34 ^ 2 / 0.2e1, (-t48 * t23 + t51 * t26) * t34 + t22 * t31, -(t51 * t23 + t48 * t26) * t34 + t22 * t32;];
T_reg = t1;
