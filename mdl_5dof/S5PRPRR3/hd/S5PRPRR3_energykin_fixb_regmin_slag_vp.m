% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% T_reg [1x19]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPRR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR3_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:47:31
% EndTime: 2019-12-05 15:47:31
% DurationCPUTime: 0.06s
% Computational Cost: add. (69->27), mult. (172->67), div. (0->0), fcn. (107->8), ass. (0->29)
t119 = qJD(2) ^ 2;
t126 = t119 / 0.2e1;
t114 = sin(qJ(4));
t117 = cos(qJ(4));
t118 = cos(qJ(2));
t104 = qJD(2) * pkin(2) + t118 * qJD(1);
t111 = sin(pkin(9));
t112 = cos(pkin(9));
t115 = sin(qJ(2));
t124 = qJD(1) * t115;
t100 = t111 * t104 + t112 * t124;
t98 = qJD(2) * pkin(6) + t100;
t125 = t114 * qJD(3) + t117 * t98;
t123 = qJD(2) * t114;
t122 = qJD(2) * t117;
t121 = qJD(1) * qJD(2);
t120 = qJD(2) * qJD(4);
t99 = t112 * t104 - t111 * t124;
t116 = cos(qJ(5));
t113 = sin(qJ(5));
t110 = qJD(4) + qJD(5);
t109 = t117 * qJD(3);
t102 = (t113 * t117 + t114 * t116) * qJD(2);
t101 = t113 * t123 - t116 * t122;
t97 = -qJD(2) * pkin(3) - t99;
t95 = (-pkin(4) * t117 - pkin(3)) * qJD(2) - t99;
t94 = pkin(7) * t122 + t125;
t93 = qJD(4) * pkin(4) + t109 + (-pkin(7) * qJD(2) - t98) * t114;
t1 = [qJD(1) ^ 2 / 0.2e1, t126, t118 * t121, -t115 * t121, t100 ^ 2 / 0.2e1 + t99 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t114 ^ 2 * t126, t114 * t119 * t117, t114 * t120, t117 * t120, qJD(4) ^ 2 / 0.2e1, (-t114 * t98 + t109) * qJD(4) - t97 * t122, -t125 * qJD(4) + t97 * t123, t102 ^ 2 / 0.2e1, -t102 * t101, t102 * t110, -t101 * t110, t110 ^ 2 / 0.2e1, (-t113 * t94 + t116 * t93) * t110 + t95 * t101, -(t113 * t93 + t116 * t94) * t110 + t95 * t102;];
T_reg = t1;
