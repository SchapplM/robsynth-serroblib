% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP3_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:53
% EndTime: 2019-12-31 18:12:53
% DurationCPUTime: 0.07s
% Computational Cost: add. (158->37), mult. (410->74), div. (0->0), fcn. (242->4), ass. (0->28)
t115 = pkin(3) + qJ(5);
t114 = pkin(6) + qJ(2);
t104 = sin(qJ(3));
t105 = cos(qJ(3));
t102 = sin(pkin(7));
t112 = qJD(1) * t102;
t95 = t114 * t112;
t103 = cos(pkin(7));
t111 = qJD(1) * t103;
t96 = t114 * t111;
t113 = -t104 * t95 + t105 * t96;
t110 = -t104 * t96 - t105 * t95;
t89 = -qJD(3) * qJ(4) - t113;
t109 = qJD(4) - t110;
t97 = qJD(2) + (-pkin(2) * t103 - pkin(1)) * qJD(1);
t94 = (t102 * t105 + t103 * t104) * qJD(1);
t108 = -t94 * qJ(4) + t97;
t106 = qJD(1) ^ 2;
t101 = t103 ^ 2;
t100 = t102 ^ 2;
t99 = -qJD(1) * pkin(1) + qJD(2);
t93 = t104 * t112 - t105 * t111;
t88 = -qJD(3) * pkin(3) + t109;
t87 = t93 * pkin(3) + t108;
t86 = -t93 * pkin(4) + qJD(5) - t89;
t85 = t94 * pkin(4) - t115 * qJD(3) + t109;
t84 = t115 * t93 + t108;
t1 = [t106 / 0.2e1, 0, 0, -t99 * t111, t99 * t112, (t100 + t101) * t106 * qJ(2), t99 ^ 2 / 0.2e1 + (t101 / 0.2e1 + t100 / 0.2e1) * qJ(2) ^ 2 * t106, t94 ^ 2 / 0.2e1, -t94 * t93, t94 * qJD(3), -t93 * qJD(3), qJD(3) ^ 2 / 0.2e1, t110 * qJD(3) + t97 * t93, -t113 * qJD(3) + t97 * t94, t88 * t94 + t89 * t93, t88 * qJD(3) - t87 * t93, -t89 * qJD(3) - t87 * t94, t87 ^ 2 / 0.2e1 + t89 ^ 2 / 0.2e1 + t88 ^ 2 / 0.2e1, t85 * t94 - t86 * t93, t86 * qJD(3) - t84 * t94, -t85 * qJD(3) + t84 * t93, t84 ^ 2 / 0.2e1 + t85 ^ 2 / 0.2e1 + t86 ^ 2 / 0.2e1;];
T_reg = t1;
