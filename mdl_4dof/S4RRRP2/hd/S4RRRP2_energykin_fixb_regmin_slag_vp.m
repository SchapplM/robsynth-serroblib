% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% T_reg [1x17]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:04
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRRP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:04:28
% EndTime: 2021-01-15 11:04:28
% DurationCPUTime: 0.04s
% Computational Cost: add. (71->16), mult. (113->46), div. (0->0), fcn. (43->4), ass. (0->19)
t9 = qJD(1) + qJD(2);
t8 = t9 ^ 2;
t22 = t8 / 0.2e1;
t10 = sin(qJ(3));
t21 = t10 * t9;
t12 = cos(qJ(3));
t20 = t12 * t9;
t19 = pkin(1) * qJD(1);
t18 = qJD(3) * t10;
t17 = qJD(3) * t12;
t16 = sin(qJ(2)) * t19;
t15 = cos(qJ(2)) * t19;
t6 = t9 * pkin(6) + t16;
t14 = qJ(4) * t9 + t6;
t7 = -t9 * pkin(2) - t15;
t5 = -t15 + qJD(4) + (-pkin(3) * t12 - pkin(2)) * t9;
t4 = t14 * t12;
t3 = qJD(3) * pkin(3) - t14 * t10;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t22, t9 * t15, -t9 * t16, t10 ^ 2 * t22, t10 * t8 * t12, t9 * t18, t9 * t17, qJD(3) ^ 2 / 0.2e1, -t6 * t18 - t7 * t20, -t6 * t17 + t7 * t21, t3 * qJD(3) - t5 * t20, -t4 * qJD(3) + t5 * t21, (-t10 * t3 + t12 * t4) * t9, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t1;
