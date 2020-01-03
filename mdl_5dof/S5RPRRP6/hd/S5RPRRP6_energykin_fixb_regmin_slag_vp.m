% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% T_reg [1x20]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP6_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:43:15
% EndTime: 2019-12-31 18:43:15
% DurationCPUTime: 0.11s
% Computational Cost: add. (95->29), mult. (234->70), div. (0->0), fcn. (123->6), ass. (0->28)
t117 = qJD(1) ^ 2;
t128 = t117 / 0.2e1;
t127 = cos(qJ(4));
t112 = sin(pkin(8));
t106 = (pkin(1) * t112 + pkin(6)) * qJD(1);
t115 = sin(qJ(3));
t116 = cos(qJ(3));
t125 = t115 * qJD(2) + t116 * t106;
t100 = qJD(3) * pkin(7) + t125;
t113 = cos(pkin(8));
t121 = -pkin(1) * t113 - pkin(2);
t101 = (-pkin(3) * t116 - pkin(7) * t115 + t121) * qJD(1);
t114 = sin(qJ(4));
t126 = t127 * t100 + t114 * t101;
t124 = qJD(1) * t115;
t123 = t116 * qJD(1);
t122 = qJD(1) * qJD(3);
t120 = -t114 * t100 + t127 * t101;
t119 = t116 * qJD(2) - t115 * t106;
t99 = -qJD(3) * pkin(3) - t119;
t108 = -qJD(4) + t123;
t107 = t121 * qJD(1);
t105 = t114 * qJD(3) + t127 * t124;
t104 = -t127 * qJD(3) + t114 * t124;
t95 = t104 * pkin(4) + qJD(5) + t99;
t94 = -t104 * qJ(5) + t126;
t93 = -t108 * pkin(4) - t105 * qJ(5) + t120;
t1 = [t128, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t112 ^ 2 / 0.2e1 + t113 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t117, t115 ^ 2 * t128, t115 * t117 * t116, t115 * t122, t116 * t122, qJD(3) ^ 2 / 0.2e1, t119 * qJD(3) - t107 * t123, -t125 * qJD(3) + t107 * t124, t105 ^ 2 / 0.2e1, -t105 * t104, -t105 * t108, t104 * t108, t108 ^ 2 / 0.2e1, t99 * t104 - t120 * t108, t99 * t105 + t126 * t108, -t94 * t104 - t93 * t105, t94 ^ 2 / 0.2e1 + t93 ^ 2 / 0.2e1 + t95 ^ 2 / 0.2e1;];
T_reg = t1;
