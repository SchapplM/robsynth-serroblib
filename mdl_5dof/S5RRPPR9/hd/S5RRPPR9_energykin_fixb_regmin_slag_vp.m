% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPPR9_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR9_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:41:47
% EndTime: 2019-12-31 19:41:47
% DurationCPUTime: 0.11s
% Computational Cost: add. (101->35), mult. (235->79), div. (0->0), fcn. (98->4), ass. (0->28)
t127 = -pkin(2) - pkin(3);
t117 = qJD(1) ^ 2;
t126 = t117 / 0.2e1;
t116 = cos(qJ(2));
t125 = t116 * t117;
t124 = qJD(1) * t116;
t104 = pkin(6) * t124 + qJD(2) * qJ(3);
t114 = sin(qJ(2));
t123 = t114 * qJD(1);
t122 = pkin(6) * t123 + qJD(3);
t121 = qJD(1) * qJD(2);
t100 = -qJD(1) * pkin(1) - pkin(2) * t124 - qJ(3) * t123;
t120 = t114 * t121;
t119 = t116 * t121;
t96 = pkin(3) * t124 + qJD(4) - t100;
t99 = qJ(4) * t124 - t104;
t118 = -qJ(4) * t123 + t122;
t115 = cos(qJ(5));
t113 = sin(qJ(5));
t105 = qJD(5) + t123;
t103 = -qJD(2) * pkin(2) + t122;
t102 = -t113 * qJD(2) - t115 * t124;
t101 = -t115 * qJD(2) + t113 * t124;
t98 = qJD(2) * pkin(4) - t99;
t97 = t127 * qJD(2) + t118;
t95 = (-pkin(7) + t127) * qJD(2) + t118;
t94 = (pkin(4) * t114 + pkin(7) * t116) * qJD(1) + t96;
t1 = [t126, 0, 0, t114 ^ 2 * t126, t114 * t125, t120, t119, qJD(2) ^ 2 / 0.2e1, pkin(1) * t125 - pkin(6) * t120, -t117 * pkin(1) * t114 - pkin(6) * t119, -t103 * qJD(2) - t100 * t124, (t103 * t114 + t104 * t116) * qJD(1), t104 * qJD(2) - t100 * t123, t104 ^ 2 / 0.2e1 + t100 ^ 2 / 0.2e1 + t103 ^ 2 / 0.2e1, -t99 * qJD(2) + t96 * t123, t97 * qJD(2) - t96 * t124, (-t114 * t97 + t116 * t99) * qJD(1), t97 ^ 2 / 0.2e1 + t99 ^ 2 / 0.2e1 + t96 ^ 2 / 0.2e1, t102 ^ 2 / 0.2e1, t102 * t101, t102 * t105, t101 * t105, t105 ^ 2 / 0.2e1, (-t113 * t95 + t115 * t94) * t105 - t98 * t101, -(t113 * t94 + t115 * t95) * t105 + t98 * t102;];
T_reg = t1;
