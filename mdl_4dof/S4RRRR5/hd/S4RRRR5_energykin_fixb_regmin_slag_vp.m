% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRRR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR5_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:28:13
% EndTime: 2019-12-31 17:28:13
% DurationCPUTime: 0.06s
% Computational Cost: add. (95->28), mult. (240->67), div. (0->0), fcn. (154->6), ass. (0->29)
t114 = qJD(1) ^ 2;
t124 = t114 / 0.2e1;
t123 = cos(qJ(3));
t113 = cos(qJ(2));
t119 = t113 * qJD(1);
t102 = pkin(5) * t119 + qJD(2) * pkin(6);
t110 = sin(qJ(3));
t111 = sin(qJ(2));
t97 = (-pkin(2) * t113 - pkin(6) * t111 - pkin(1)) * qJD(1);
t122 = t123 * t102 + t110 * t97;
t121 = t113 * t114;
t120 = qJD(1) * t111;
t118 = qJD(1) * qJD(2);
t117 = t111 * t118;
t116 = t113 * t118;
t101 = -qJD(2) * pkin(2) + pkin(5) * t120;
t115 = -t110 * t102 + t123 * t97;
t105 = -qJD(3) + t119;
t112 = cos(qJ(4));
t109 = sin(qJ(4));
t103 = -qJD(4) + t105;
t99 = t110 * qJD(2) + t123 * t120;
t98 = -t123 * qJD(2) + t110 * t120;
t93 = t98 * pkin(3) + t101;
t92 = -t109 * t98 + t112 * t99;
t91 = t109 * t99 + t112 * t98;
t90 = -t98 * pkin(7) + t122;
t89 = -t105 * pkin(3) - t99 * pkin(7) + t115;
t1 = [t124, 0, 0, t111 ^ 2 * t124, t111 * t121, t117, t116, qJD(2) ^ 2 / 0.2e1, pkin(1) * t121 - pkin(5) * t117, -t114 * pkin(1) * t111 - pkin(5) * t116, t99 ^ 2 / 0.2e1, -t99 * t98, -t99 * t105, t98 * t105, t105 ^ 2 / 0.2e1, t101 * t98 - t115 * t105, t101 * t99 + t122 * t105, t92 ^ 2 / 0.2e1, -t92 * t91, -t92 * t103, t91 * t103, t103 ^ 2 / 0.2e1, -(-t109 * t90 + t112 * t89) * t103 + t93 * t91, (t109 * t89 + t112 * t90) * t103 + t93 * t92;];
T_reg = t1;
