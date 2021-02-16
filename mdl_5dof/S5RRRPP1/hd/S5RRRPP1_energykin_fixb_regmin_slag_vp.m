% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRPP1
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
% T_reg [1x21]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:14:46
% EndTime: 2021-01-15 22:14:46
% DurationCPUTime: 0.06s
% Computational Cost: add. (237->30), mult. (345->70), div. (0->0), fcn. (179->6), ass. (0->28)
t29 = qJD(1) + qJD(2);
t28 = t29 ^ 2;
t44 = t28 / 0.2e1;
t32 = sin(qJ(3));
t43 = t29 * t32;
t34 = cos(qJ(3));
t42 = t29 * t34;
t41 = pkin(1) * qJD(1);
t38 = sin(qJ(2)) * t41;
t25 = t29 * pkin(7) + t38;
t36 = qJ(4) * t29 + t25;
t20 = qJD(3) * pkin(3) - t36 * t32;
t21 = t36 * t34;
t30 = sin(pkin(8));
t31 = cos(pkin(8));
t17 = t30 * t20 + t31 * t21;
t40 = qJD(3) * t32;
t39 = qJD(3) * t34;
t37 = cos(qJ(2)) * t41;
t16 = t31 * t20 - t30 * t21;
t22 = -t37 + qJD(4) + (-pkin(3) * t34 - pkin(2)) * t29;
t26 = -t29 * pkin(2) - t37;
t24 = (t30 * t34 + t31 * t32) * t29;
t23 = t30 * t43 - t31 * t42;
t15 = qJD(3) * qJ(5) + t17;
t14 = -qJD(3) * pkin(4) + qJD(5) - t16;
t13 = t23 * pkin(4) - t24 * qJ(5) + t22;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t44, t29 * t37, -t29 * t38, t32 ^ 2 * t44, t32 * t28 * t34, t29 * t40, t29 * t39, qJD(3) ^ 2 / 0.2e1, -t25 * t40 - t26 * t42, -t25 * t39 + t26 * t43, t16 * qJD(3) + t22 * t23, -t17 * qJD(3) + t22 * t24, -t16 * t24 - t17 * t23, t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, -t14 * qJD(3) + t13 * t23, t14 * t24 - t15 * t23, t15 * qJD(3) - t13 * t24, t15 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1;];
T_reg = t1;
