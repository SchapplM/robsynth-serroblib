% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% T_reg [1x21]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPRR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR7_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:54:04
% EndTime: 2019-12-31 16:54:04
% DurationCPUTime: 0.06s
% Computational Cost: add. (82->29), mult. (243->65), div. (0->0), fcn. (156->6), ass. (0->27)
t109 = pkin(5) + qJ(2);
t100 = sin(qJ(3));
t102 = cos(qJ(3));
t97 = sin(pkin(7));
t107 = qJD(1) * t97;
t89 = t109 * t107;
t98 = cos(pkin(7));
t106 = qJD(1) * t98;
t90 = t109 * t106;
t108 = -t100 * t89 + t102 * t90;
t87 = t100 * t107 - t102 * t106;
t105 = -t100 * t90 - t102 * t89;
t91 = qJD(2) + (-pkin(2) * t98 - pkin(1)) * qJD(1);
t103 = qJD(1) ^ 2;
t101 = cos(qJ(4));
t99 = sin(qJ(4));
t96 = t98 ^ 2;
t95 = t97 ^ 2;
t93 = -qJD(1) * pkin(1) + qJD(2);
t88 = (t100 * t98 + t102 * t97) * qJD(1);
t84 = qJD(4) + t87;
t83 = t99 * qJD(3) + t101 * t88;
t82 = -t101 * qJD(3) + t99 * t88;
t81 = qJD(3) * pkin(6) + t108;
t80 = -qJD(3) * pkin(3) - t105;
t79 = t87 * pkin(3) - t88 * pkin(6) + t91;
t1 = [t103 / 0.2e1, 0, 0, -t93 * t106, t93 * t107, (t95 + t96) * t103 * qJ(2), t93 ^ 2 / 0.2e1 + (t96 / 0.2e1 + t95 / 0.2e1) * qJ(2) ^ 2 * t103, t88 ^ 2 / 0.2e1, -t88 * t87, t88 * qJD(3), -t87 * qJD(3), qJD(3) ^ 2 / 0.2e1, t105 * qJD(3) + t91 * t87, -t108 * qJD(3) + t91 * t88, t83 ^ 2 / 0.2e1, -t83 * t82, t83 * t84, -t82 * t84, t84 ^ 2 / 0.2e1, (t101 * t79 - t99 * t81) * t84 + t80 * t82, -(t101 * t81 + t99 * t79) * t84 + t80 * t83;];
T_reg = t1;
