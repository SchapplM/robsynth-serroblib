% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% T_reg [1x21]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR7_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:52
% EndTime: 2019-12-31 17:59:52
% DurationCPUTime: 0.05s
% Computational Cost: add. (67->29), mult. (156->64), div. (0->0), fcn. (69->6), ass. (0->27)
t99 = qJD(1) ^ 2;
t109 = t99 / 0.2e1;
t94 = cos(pkin(8));
t103 = -pkin(1) * t94 - pkin(2);
t83 = qJD(3) + (-pkin(6) + t103) * qJD(1);
t96 = sin(qJ(4));
t98 = cos(qJ(4));
t108 = t98 * qJD(2) + t96 * t83;
t107 = qJD(1) * t98;
t93 = sin(pkin(8));
t102 = -pkin(1) * t93 - qJ(3);
t87 = t102 * qJD(1);
t106 = t87 * qJD(1);
t105 = t96 * qJD(1);
t104 = qJD(1) * qJD(4);
t101 = -t96 * qJD(2) + t98 * t83;
t97 = cos(qJ(5));
t95 = sin(qJ(5));
t92 = qJD(2) ^ 2 / 0.2e1;
t89 = qJD(5) + t105;
t86 = t103 * qJD(1) + qJD(3);
t85 = t95 * qJD(4) + t97 * t107;
t84 = -t97 * qJD(4) + t95 * t107;
t81 = (pkin(4) * t96 - pkin(7) * t98 - t102) * qJD(1);
t80 = qJD(4) * pkin(7) + t108;
t79 = -qJD(4) * pkin(4) - t101;
t1 = [t109, 0, 0, t92 + (t93 ^ 2 / 0.2e1 + t94 ^ 2 / 0.2e1) * t99 * pkin(1) ^ 2, t86 * qJD(1), -t106, t92 + t87 ^ 2 / 0.2e1 + t86 ^ 2 / 0.2e1, t98 ^ 2 * t109, -t98 * t99 * t96, t98 * t104, -t96 * t104, qJD(4) ^ 2 / 0.2e1, t101 * qJD(4) - t87 * t105, -t108 * qJD(4) - t98 * t106, t85 ^ 2 / 0.2e1, -t85 * t84, t85 * t89, -t84 * t89, t89 ^ 2 / 0.2e1, (-t95 * t80 + t97 * t81) * t89 + t79 * t84, -(t97 * t80 + t95 * t81) * t89 + t79 * t85;];
T_reg = t1;
