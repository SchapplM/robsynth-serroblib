% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPPRP1
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
% T_reg [1x17]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:25:48
% EndTime: 2020-01-03 11:25:48
% DurationCPUTime: 0.06s
% Computational Cost: add. (90->30), mult. (221->70), div. (0->0), fcn. (112->6), ass. (0->27)
t106 = sin(qJ(4));
t107 = cos(qJ(4));
t102 = sin(pkin(8));
t104 = cos(pkin(8));
t105 = cos(pkin(7));
t113 = -pkin(1) * t105 - pkin(2);
t91 = qJD(3) + (-pkin(3) * t104 - pkin(6) * t102 + t113) * qJD(1);
t103 = sin(pkin(7));
t97 = (pkin(1) * t103 + qJ(3)) * qJD(1);
t95 = t102 * qJD(2) + t104 * t97;
t118 = t106 * t91 + t107 * t95;
t108 = qJD(1) ^ 2;
t117 = t102 ^ 2 * t108;
t116 = qJD(1) * t102;
t115 = qJD(1) * t106;
t114 = t104 * qJD(1);
t112 = t102 * t115;
t111 = t107 * t116;
t110 = -t106 * t95 + t107 * t91;
t100 = t104 * qJD(2);
t98 = -qJD(4) + t114;
t96 = t113 * qJD(1) + qJD(3);
t93 = t102 * t97 - t100;
t88 = qJD(5) - t100 + (pkin(4) * t115 + t97) * t102;
t87 = -qJ(5) * t112 + t118;
t86 = -t98 * pkin(4) - qJ(5) * t111 + t110;
t1 = [t108 / 0.2e1, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t103 ^ 2 / 0.2e1 + t105 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t108, -t96 * t114, t96 * t116, (t102 * t93 + t104 * t95) * qJD(1), t95 ^ 2 / 0.2e1 + t93 ^ 2 / 0.2e1 + t96 ^ 2 / 0.2e1, t107 ^ 2 * t117 / 0.2e1, -t107 * t106 * t117, -t98 * t111, t98 * t112, t98 ^ 2 / 0.2e1, -t110 * t98 + t93 * t112, t93 * t111 + t118 * t98, (-t106 * t87 - t107 * t86) * t116, t87 ^ 2 / 0.2e1 + t86 ^ 2 / 0.2e1 + t88 ^ 2 / 0.2e1;];
T_reg = t1;
