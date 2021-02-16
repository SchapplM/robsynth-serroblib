% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:53
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:52:43
% EndTime: 2021-01-15 23:52:43
% DurationCPUTime: 0.09s
% Computational Cost: add. (263->37), mult. (649->82), div. (0->0), fcn. (452->6), ass. (0->35)
t45 = qJD(1) ^ 2;
t58 = t45 / 0.2e1;
t57 = pkin(7) + pkin(6);
t56 = cos(qJ(4));
t44 = cos(qJ(2));
t55 = t44 * t45;
t41 = sin(qJ(3));
t42 = sin(qJ(2));
t43 = cos(qJ(3));
t32 = (t41 * t44 + t42 * t43) * qJD(1);
t39 = qJD(2) + qJD(3);
t52 = qJD(1) * t42;
t34 = qJD(2) * pkin(2) - t57 * t52;
t51 = qJD(1) * t44;
t35 = t57 * t51;
t46 = t43 * t34 - t41 * t35;
t22 = t39 * pkin(3) - t32 * pkin(8) + t46;
t31 = t41 * t52 - t43 * t51;
t53 = t41 * t34 + t43 * t35;
t24 = -t31 * pkin(8) + t53;
t40 = sin(qJ(4));
t54 = t40 * t22 + t56 * t24;
t50 = qJD(1) * qJD(2);
t49 = t42 * t50;
t48 = t44 * t50;
t47 = t56 * t22 - t40 * t24;
t36 = (-pkin(2) * t44 - pkin(1)) * qJD(1);
t27 = t31 * pkin(3) + t36;
t38 = qJD(4) + t39;
t26 = -t40 * t31 + t56 * t32;
t25 = t56 * t31 + t40 * t32;
t19 = t25 * pkin(4) + qJD(5) + t27;
t18 = -t25 * qJ(5) + t54;
t17 = t38 * pkin(4) - t26 * qJ(5) + t47;
t1 = [t58, 0, 0, t42 ^ 2 * t58, t42 * t55, t49, t48, qJD(2) ^ 2 / 0.2e1, pkin(1) * t55 - pkin(6) * t49, -t45 * pkin(1) * t42 - pkin(6) * t48, t32 ^ 2 / 0.2e1, -t32 * t31, t32 * t39, -t31 * t39, t39 ^ 2 / 0.2e1, t36 * t31 + t46 * t39, t36 * t32 - t53 * t39, t26 ^ 2 / 0.2e1, -t26 * t25, t26 * t38, -t25 * t38, t38 ^ 2 / 0.2e1, t27 * t25 + t47 * t38, t27 * t26 - t54 * t38, t17 * t38 + t19 * t25, -t18 * t38 + t19 * t26, -t17 * t26 - t18 * t25, t18 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1;];
T_reg = t1;
