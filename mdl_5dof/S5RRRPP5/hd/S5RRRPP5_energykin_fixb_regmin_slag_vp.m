% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPP5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP5_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:58:27
% EndTime: 2019-12-31 20:58:27
% DurationCPUTime: 0.13s
% Computational Cost: add. (183->35), mult. (420->75), div. (0->0), fcn. (245->4), ass. (0->31)
t109 = cos(qJ(2));
t103 = (pkin(2) * t109 + pkin(1)) * qJD(1);
t107 = sin(qJ(3));
t108 = sin(qJ(2));
t122 = cos(qJ(3));
t98 = (t107 * t109 + t122 * t108) * qJD(1);
t126 = -qJ(4) * t98 - t103;
t125 = -pkin(3) - pkin(4);
t124 = -pkin(7) - pkin(6);
t110 = qJD(1) ^ 2;
t123 = t110 / 0.2e1;
t118 = qJD(1) * t108;
t101 = qJD(2) * pkin(2) + t124 * t118;
t117 = qJD(1) * t109;
t102 = t124 * t117;
t120 = t107 * t101 - t122 * t102;
t119 = t109 * t110;
t116 = qJD(1) * qJD(2);
t106 = qJD(2) + qJD(3);
t95 = t106 * qJ(4) + t120;
t114 = t108 * t116;
t113 = t109 * t116;
t112 = t122 * t101 + t107 * t102;
t111 = qJD(4) - t112;
t97 = t107 * t118 - t122 * t117;
t94 = -t106 * pkin(3) + t111;
t93 = pkin(3) * t97 + t126;
t92 = qJ(5) * t97 + t95;
t91 = t125 * t97 + qJD(5) - t126;
t90 = -t98 * qJ(5) + t125 * t106 + t111;
t1 = [t123, 0, 0, t108 ^ 2 * t123, t108 * t119, t114, t113, qJD(2) ^ 2 / 0.2e1, pkin(1) * t119 - pkin(6) * t114, -pkin(1) * t108 * t110 - pkin(6) * t113, t98 ^ 2 / 0.2e1, -t98 * t97, t98 * t106, -t97 * t106, t106 ^ 2 / 0.2e1, -t103 * t97 + t112 * t106, -t103 * t98 - t120 * t106, -t106 * t94 + t93 * t97, t94 * t98 - t95 * t97, t106 * t95 - t93 * t98, t95 ^ 2 / 0.2e1 + t93 ^ 2 / 0.2e1 + t94 ^ 2 / 0.2e1, -t106 * t90 - t91 * t97, t106 * t92 + t91 * t98, -t90 * t98 + t92 * t97, t92 ^ 2 / 0.2e1 + t90 ^ 2 / 0.2e1 + t91 ^ 2 / 0.2e1;];
T_reg = t1;
