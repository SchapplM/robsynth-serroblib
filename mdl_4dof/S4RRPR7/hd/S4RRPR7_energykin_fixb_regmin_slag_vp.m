% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% T_reg [1x19]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRPR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:06:43
% EndTime: 2019-12-31 17:06:43
% DurationCPUTime: 0.06s
% Computational Cost: add. (85->26), mult. (241->63), div. (0->0), fcn. (150->6), ass. (0->29)
t112 = qJD(1) ^ 2;
t120 = t112 / 0.2e1;
t119 = pkin(5) + qJ(3);
t109 = sin(qJ(2));
t117 = qJD(1) * t109;
t101 = qJD(2) * pkin(2) - t119 * t117;
t111 = cos(qJ(2));
t116 = qJD(1) * t111;
t102 = t119 * t116;
t106 = sin(pkin(7));
t107 = cos(pkin(7));
t93 = t106 * t101 + t107 * t102;
t118 = t111 * t112;
t115 = qJD(1) * qJD(2);
t114 = t109 * t115;
t113 = t111 * t115;
t98 = -t106 * t117 + t107 * t116;
t92 = t107 * t101 - t106 * t102;
t103 = qJD(3) + (-pkin(2) * t111 - pkin(1)) * qJD(1);
t110 = cos(qJ(4));
t108 = sin(qJ(4));
t99 = (t106 * t111 + t107 * t109) * qJD(1);
t97 = qJD(4) - t98;
t95 = t108 * qJD(2) + t110 * t99;
t94 = -t110 * qJD(2) + t108 * t99;
t91 = qJD(2) * pkin(6) + t93;
t90 = -qJD(2) * pkin(3) - t92;
t89 = -t98 * pkin(3) - t99 * pkin(6) + t103;
t1 = [t120, 0, 0, t109 ^ 2 * t120, t109 * t118, t114, t113, qJD(2) ^ 2 / 0.2e1, pkin(1) * t118 - pkin(5) * t114, -t112 * pkin(1) * t109 - pkin(5) * t113, -t92 * t99 + t93 * t98, t93 ^ 2 / 0.2e1 + t92 ^ 2 / 0.2e1 + t103 ^ 2 / 0.2e1, t95 ^ 2 / 0.2e1, -t95 * t94, t95 * t97, -t94 * t97, t97 ^ 2 / 0.2e1, (-t108 * t91 + t110 * t89) * t97 + t90 * t94, -(t108 * t89 + t110 * t91) * t97 + t90 * t95;];
T_reg = t1;
