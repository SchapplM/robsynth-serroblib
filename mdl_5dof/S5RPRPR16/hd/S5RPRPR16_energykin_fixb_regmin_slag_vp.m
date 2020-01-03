% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRPR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR16_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR16_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:39:32
% EndTime: 2019-12-31 18:39:32
% DurationCPUTime: 0.11s
% Computational Cost: add. (85->34), mult. (176->72), div. (0->0), fcn. (69->4), ass. (0->26)
t104 = qJD(1) ^ 2;
t113 = t104 / 0.2e1;
t101 = sin(qJ(3));
t109 = qJD(1) * t101;
t112 = pkin(3) * t109 + qJD(1) * qJ(2);
t94 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t111 = qJD(3) * t94;
t110 = t104 * qJ(2);
t108 = qJD(3) * qJ(4);
t103 = cos(qJ(3));
t107 = t103 * qJD(1);
t106 = qJD(1) * qJD(3);
t105 = pkin(4) * qJD(1) - t94;
t102 = cos(qJ(5));
t100 = sin(qJ(5));
t98 = -pkin(1) * qJD(1) + qJD(2);
t96 = qJD(5) + t107;
t93 = qJD(3) * t102 + t100 * t109;
t92 = qJD(3) * t100 - t102 * t109;
t91 = -t101 * t94 - t108;
t90 = -qJ(4) * t107 + t112;
t89 = -qJD(3) * pkin(3) - t103 * t94 + qJD(4);
t88 = -t101 * t105 + t108;
t87 = (pkin(7) * t101 - qJ(4) * t103) * qJD(1) + t112;
t86 = qJD(4) + t105 * t103 + (-pkin(3) - pkin(7)) * qJD(3);
t1 = [t113, 0, 0, t98 * qJD(1), t110, qJ(2) ^ 2 * t113 + t98 ^ 2 / 0.2e1, t103 ^ 2 * t113, -t103 * t104 * t101, t103 * t106, -t101 * t106, qJD(3) ^ 2 / 0.2e1, t101 * t110 + t103 * t111, -t101 * t111 + t103 * t110, (t101 * t91 + t103 * t89) * qJD(1), qJD(3) * t89 - t109 * t90, -qJD(3) * t91 - t107 * t90, t90 ^ 2 / 0.2e1 + t91 ^ 2 / 0.2e1 + t89 ^ 2 / 0.2e1, t93 ^ 2 / 0.2e1, -t93 * t92, t93 * t96, -t92 * t96, t96 ^ 2 / 0.2e1, (-t100 * t87 + t102 * t86) * t96 + t88 * t92, -(t100 * t86 + t102 * t87) * t96 + t88 * t93;];
T_reg = t1;
