% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:42:57
% EndTime: 2019-12-05 18:42:57
% DurationCPUTime: 0.07s
% Computational Cost: add. (208->31), mult. (311->73), div. (0->0), fcn. (191->8), ass. (0->31)
t117 = qJD(1) + qJD(2);
t115 = t117 ^ 2;
t133 = t115 / 0.2e1;
t132 = pkin(1) * qJD(1);
t127 = cos(qJ(2)) * t132;
t131 = (-t117 * pkin(2) - t127) * t117;
t121 = sin(qJ(3));
t128 = sin(qJ(2)) * t132;
t113 = t117 * pkin(7) + t128;
t126 = qJ(4) * t117 + t113;
t107 = qJD(3) * pkin(3) - t126 * t121;
t124 = cos(qJ(3));
t109 = t126 * t124;
t118 = sin(pkin(9));
t119 = cos(pkin(9));
t100 = t118 * t107 + t119 * t109;
t130 = qJD(3) * t121;
t129 = qJD(3) * t124;
t99 = t119 * t107 - t118 * t109;
t110 = -t127 + qJD(4) + (-pkin(3) * t124 - pkin(2)) * t117;
t123 = cos(qJ(5));
t120 = sin(qJ(5));
t116 = qJD(3) + qJD(5);
t112 = (t118 * t124 + t119 * t121) * t117;
t111 = (-t118 * t121 + t119 * t124) * t117;
t103 = -t111 * pkin(4) + t110;
t102 = t120 * t111 + t123 * t112;
t101 = -t123 * t111 + t120 * t112;
t98 = t111 * pkin(8) + t100;
t97 = qJD(3) * pkin(4) - t112 * pkin(8) + t99;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t133, t117 * t127, -t117 * t128, t121 ^ 2 * t133, t121 * t115 * t124, t117 * t130, t117 * t129, qJD(3) ^ 2 / 0.2e1, -t113 * t130 - t124 * t131, -t113 * t129 + t121 * t131, t100 * t111 - t99 * t112, t100 ^ 2 / 0.2e1 + t99 ^ 2 / 0.2e1 + t110 ^ 2 / 0.2e1, t102 ^ 2 / 0.2e1, -t102 * t101, t102 * t116, -t101 * t116, t116 ^ 2 / 0.2e1, t103 * t101 + (-t120 * t98 + t123 * t97) * t116, t103 * t102 - (t120 * t97 + t123 * t98) * t116;];
T_reg = t1;
