% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% 
% Output:
% T_reg [1x19]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:49
% EndTime: 2019-12-31 18:14:49
% DurationCPUTime: 0.10s
% Computational Cost: add. (134->29), mult. (262->65), div. (0->0), fcn. (121->4), ass. (0->24)
t105 = qJD(1) ^ 2;
t111 = t105 / 0.2e1;
t101 = sin(pkin(7));
t102 = cos(pkin(7));
t104 = cos(qJ(3));
t96 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t106 = -qJ(4) * qJD(1) + t96;
t91 = qJD(3) * pkin(3) + t106 * t104;
t103 = sin(qJ(3));
t92 = t106 * t103;
t87 = t101 * t91 + t102 * t92;
t110 = qJD(3) * t96;
t109 = t105 * qJ(2);
t108 = qJD(1) * t103;
t107 = qJD(1) * qJD(3);
t95 = pkin(3) * t108 + qJD(1) * qJ(2) + qJD(4);
t86 = -t101 * t92 + t102 * t91;
t99 = -pkin(1) * qJD(1) + qJD(2);
t94 = t102 * t104 * qJD(1) - t101 * t108;
t93 = (t101 * t104 + t102 * t103) * qJD(1);
t88 = t93 * pkin(4) - t94 * qJ(5) + t95;
t85 = qJD(3) * qJ(5) + t87;
t84 = -qJD(3) * pkin(4) + qJD(5) - t86;
t1 = [t111, 0, 0, t99 * qJD(1), t109, qJ(2) ^ 2 * t111 + t99 ^ 2 / 0.2e1, t104 ^ 2 * t111, -t104 * t105 * t103, t104 * t107, -t103 * t107, qJD(3) ^ 2 / 0.2e1, t103 * t109 + t104 * t110, -t103 * t110 + t104 * t109, -t86 * t94 - t87 * t93, t87 ^ 2 / 0.2e1 + t86 ^ 2 / 0.2e1 + t95 ^ 2 / 0.2e1, -t84 * qJD(3) + t88 * t93, t84 * t94 - t85 * t93, t85 * qJD(3) - t88 * t94, t85 ^ 2 / 0.2e1 + t88 ^ 2 / 0.2e1 + t84 ^ 2 / 0.2e1;];
T_reg = t1;
