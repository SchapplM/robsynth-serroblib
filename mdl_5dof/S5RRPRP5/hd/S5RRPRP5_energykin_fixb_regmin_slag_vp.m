% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRP5
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
% Datum: 2021-01-15 20:19
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRP5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:18:17
% EndTime: 2021-01-15 20:18:17
% DurationCPUTime: 0.08s
% Computational Cost: add. (259->38), mult. (655->83), div. (0->0), fcn. (435->6), ass. (0->34)
t43 = qJD(1) ^ 2;
t53 = t43 / 0.2e1;
t42 = cos(qJ(2));
t52 = t42 * t43;
t51 = pkin(6) + qJ(3);
t40 = sin(qJ(2));
t49 = qJD(1) * t40;
t32 = qJD(2) * pkin(2) - t51 * t49;
t48 = qJD(1) * t42;
t33 = t51 * t48;
t37 = sin(pkin(8));
t38 = cos(pkin(8));
t23 = t38 * t32 - t37 * t33;
t30 = (t37 * t42 + t38 * t40) * qJD(1);
t19 = qJD(2) * pkin(3) - t30 * pkin(7) + t23;
t24 = t37 * t32 + t38 * t33;
t29 = t37 * t49 - t38 * t48;
t20 = -t29 * pkin(7) + t24;
t39 = sin(qJ(4));
t41 = cos(qJ(4));
t50 = t39 * t19 + t41 * t20;
t47 = qJD(1) * qJD(2);
t46 = t40 * t47;
t45 = t42 * t47;
t44 = t41 * t19 - t39 * t20;
t34 = qJD(3) + (-pkin(2) * t42 - pkin(1)) * qJD(1);
t25 = t29 * pkin(3) + t34;
t36 = qJD(2) + qJD(4);
t22 = -t39 * t29 + t41 * t30;
t21 = t41 * t29 + t39 * t30;
t16 = t21 * pkin(4) - t22 * qJ(5) + t25;
t15 = t36 * qJ(5) + t50;
t14 = -t36 * pkin(4) + qJD(5) - t44;
t1 = [t53, 0, 0, t40 ^ 2 * t53, t40 * t52, t46, t45, qJD(2) ^ 2 / 0.2e1, pkin(1) * t52 - pkin(6) * t46, -t43 * pkin(1) * t40 - pkin(6) * t45, t23 * qJD(2) + t34 * t29, -t24 * qJD(2) + t34 * t30, -t23 * t30 - t24 * t29, t24 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t21, t22 * t36, -t21 * t36, t36 ^ 2 / 0.2e1, t25 * t21 + t44 * t36, t25 * t22 - t50 * t36, -t14 * t36 + t16 * t21, t14 * t22 - t15 * t21, t15 * t36 - t16 * t22, t15 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1;];
T_reg = t1;
