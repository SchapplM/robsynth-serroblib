% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
% 
% Output:
% T_reg [1x31]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_energykin_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:56:28
% EndTime: 2019-12-05 18:56:28
% DurationCPUTime: 0.12s
% Computational Cost: add. (179->31), mult. (406->80), div. (0->0), fcn. (317->8), ass. (0->31)
t128 = qJD(1) ^ 2;
t135 = t128 / 0.2e1;
t134 = cos(qJ(4));
t127 = cos(qJ(2));
t133 = qJD(1) * t127;
t123 = sin(qJ(3));
t132 = qJD(2) * t123;
t126 = cos(qJ(3));
t131 = qJD(2) * t126;
t130 = qJD(1) * qJD(2);
t124 = sin(qJ(2));
t114 = t123 * t124 * qJD(1) - t126 * t133;
t115 = (t123 * t127 + t124 * t126) * qJD(1);
t107 = -pkin(1) * t133 + t114 * pkin(2) - t115 * pkin(5);
t120 = qJD(2) + qJD(3);
t116 = pkin(1) * t132 + t120 * pkin(5);
t122 = sin(qJ(4));
t129 = t134 * t107 - t122 * t116;
t117 = -pkin(1) * t131 - t120 * pkin(2);
t113 = qJD(4) + t114;
t125 = cos(qJ(5));
t121 = sin(qJ(5));
t112 = qJD(5) + t113;
t111 = t134 * t115 + t122 * t120;
t110 = t122 * t115 - t134 * t120;
t108 = t110 * pkin(3) + t117;
t105 = t122 * t107 + t134 * t116;
t104 = -t121 * t110 + t125 * t111;
t103 = t125 * t110 + t121 * t111;
t102 = t113 * pkin(3) + t129;
t1 = [t135, 0, 0, t124 ^ 2 * t135, t124 * t128 * t127, t124 * t130, t127 * t130, qJD(2) ^ 2 / 0.2e1, 0, 0, t115 ^ 2 / 0.2e1, -t115 * t114, t115 * t120, -t114 * t120, t120 ^ 2 / 0.2e1, (-t114 * t133 + t120 * t131) * pkin(1), (-t115 * t133 - t120 * t132) * pkin(1), t111 ^ 2 / 0.2e1, -t111 * t110, t111 * t113, -t110 * t113, t113 ^ 2 / 0.2e1, t117 * t110 + t129 * t113, -t105 * t113 + t117 * t111, t104 ^ 2 / 0.2e1, -t104 * t103, t104 * t112, -t103 * t112, t112 ^ 2 / 0.2e1, (t125 * t102 - t121 * t105) * t112 + t108 * t103, -(t121 * t102 + t125 * t105) * t112 + t108 * t104;];
T_reg = t1;
