% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% T_reg [1x19]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:23
% EndTime: 2019-12-31 17:49:23
% DurationCPUTime: 0.11s
% Computational Cost: add. (130->34), mult. (315->73), div. (0->0), fcn. (179->6), ass. (0->26)
t102 = sin(qJ(4));
t103 = cos(qJ(4));
t99 = sin(pkin(7));
t94 = (pkin(1) * t99 + qJ(3)) * qJD(1);
t100 = cos(pkin(8));
t97 = t100 * qJD(2);
t98 = sin(pkin(8));
t85 = t97 + (-pkin(6) * qJD(1) - t94) * t98;
t108 = qJD(1) * t100;
t88 = t98 * qJD(2) + t100 * t94;
t86 = pkin(6) * t108 + t88;
t110 = t102 * t85 + t103 * t86;
t109 = qJD(1) * t98;
t101 = cos(pkin(7));
t107 = -pkin(1) * t101 - pkin(2);
t106 = -t102 * t86 + t103 * t85;
t89 = qJD(3) + (-pkin(3) * t100 + t107) * qJD(1);
t104 = qJD(1) ^ 2;
t93 = t107 * qJD(1) + qJD(3);
t91 = (t100 * t102 + t103 * t98) * qJD(1);
t90 = t102 * t109 - t103 * t108;
t87 = -t98 * t94 + t97;
t82 = t90 * pkin(4) - t91 * qJ(5) + t89;
t81 = qJD(4) * qJ(5) + t110;
t80 = -qJD(4) * pkin(4) + qJD(5) - t106;
t1 = [t104 / 0.2e1, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t99 ^ 2 / 0.2e1 + t101 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t104, -t93 * t108, t93 * t109, (t100 * t88 - t87 * t98) * qJD(1), t88 ^ 2 / 0.2e1 + t87 ^ 2 / 0.2e1 + t93 ^ 2 / 0.2e1, t91 ^ 2 / 0.2e1, -t91 * t90, t91 * qJD(4), -t90 * qJD(4), qJD(4) ^ 2 / 0.2e1, t106 * qJD(4) + t89 * t90, -t110 * qJD(4) + t89 * t91, -t80 * qJD(4) + t82 * t90, t80 * t91 - t81 * t90, t81 * qJD(4) - t82 * t91, t81 ^ 2 / 0.2e1 + t82 ^ 2 / 0.2e1 + t80 ^ 2 / 0.2e1;];
T_reg = t1;
