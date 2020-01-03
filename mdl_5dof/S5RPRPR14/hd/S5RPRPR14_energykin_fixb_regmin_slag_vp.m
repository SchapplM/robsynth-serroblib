% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR14_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:35:15
% EndTime: 2019-12-31 18:35:15
% DurationCPUTime: 0.07s
% Computational Cost: add. (126->32), mult. (270->71), div. (0->0), fcn. (150->6), ass. (0->28)
t118 = sin(pkin(8));
t119 = cos(pkin(8));
t121 = sin(qJ(3));
t123 = cos(qJ(3));
t110 = (t118 * t123 + t119 * t121) * qJD(1);
t124 = qJD(1) ^ 2;
t130 = t124 / 0.2e1;
t129 = t124 * qJ(2);
t113 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t126 = -qJ(4) * qJD(1) + t113;
t107 = qJD(3) * pkin(3) + t126 * t123;
t108 = t126 * t121;
t101 = t118 * t107 + t119 * t108;
t128 = qJD(3) * t113;
t127 = qJD(1) * qJD(3);
t112 = qJD(4) + (pkin(3) * t121 + qJ(2)) * qJD(1);
t100 = t119 * t107 - t118 * t108;
t122 = cos(qJ(5));
t120 = sin(qJ(5));
t115 = -qJD(1) * pkin(1) + qJD(2);
t111 = (-t118 * t121 + t119 * t123) * qJD(1);
t109 = qJD(5) + t110;
t104 = t120 * qJD(3) + t122 * t111;
t103 = -t122 * qJD(3) + t120 * t111;
t102 = t110 * pkin(4) - t111 * pkin(7) + t112;
t99 = qJD(3) * pkin(7) + t101;
t98 = -qJD(3) * pkin(4) - t100;
t1 = [t130, 0, 0, t115 * qJD(1), t129, qJ(2) ^ 2 * t130 + t115 ^ 2 / 0.2e1, t123 ^ 2 * t130, -t123 * t124 * t121, t123 * t127, -t121 * t127, qJD(3) ^ 2 / 0.2e1, t121 * t129 + t123 * t128, -t121 * t128 + t123 * t129, -t100 * t111 - t101 * t110, t101 ^ 2 / 0.2e1 + t100 ^ 2 / 0.2e1 + t112 ^ 2 / 0.2e1, t104 ^ 2 / 0.2e1, -t104 * t103, t104 * t109, -t103 * t109, t109 ^ 2 / 0.2e1, (t122 * t102 - t120 * t99) * t109 + t98 * t103, -(t120 * t102 + t122 * t99) * t109 + t98 * t104;];
T_reg = t1;
