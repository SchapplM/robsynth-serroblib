% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP9_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP9_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:49:10
% EndTime: 2019-12-31 18:49:10
% DurationCPUTime: 0.13s
% Computational Cost: add. (228->39), mult. (605->81), div. (0->0), fcn. (419->6), ass. (0->32)
t134 = cos(qJ(3));
t133 = pkin(6) + qJ(2);
t120 = sin(pkin(8));
t121 = cos(pkin(8));
t123 = sin(qJ(3));
t111 = (t134 * t120 + t121 * t123) * qJD(1);
t130 = qJD(1) * t120;
t112 = t133 * t130;
t129 = qJD(1) * t121;
t113 = t133 * t129;
t128 = -t134 * t112 - t123 * t113;
t101 = qJD(3) * pkin(3) - t111 * pkin(7) + t128;
t110 = t123 * t130 - t134 * t129;
t131 = -t123 * t112 + t134 * t113;
t102 = -t110 * pkin(7) + t131;
t122 = sin(qJ(4));
t124 = cos(qJ(4));
t132 = t122 * t101 + t124 * t102;
t127 = t124 * t101 - t122 * t102;
t114 = qJD(2) + (-pkin(2) * t121 - pkin(1)) * qJD(1);
t105 = t110 * pkin(3) + t114;
t125 = qJD(1) ^ 2;
t119 = qJD(3) + qJD(4);
t118 = t121 ^ 2;
t117 = t120 ^ 2;
t116 = -qJD(1) * pkin(1) + qJD(2);
t104 = -t122 * t110 + t124 * t111;
t103 = t124 * t110 + t122 * t111;
t98 = t103 * pkin(4) - t104 * qJ(5) + t105;
t97 = t119 * qJ(5) + t132;
t96 = -t119 * pkin(4) + qJD(5) - t127;
t1 = [t125 / 0.2e1, 0, 0, -t116 * t129, t116 * t130, (t117 + t118) * t125 * qJ(2), t116 ^ 2 / 0.2e1 + (t118 / 0.2e1 + t117 / 0.2e1) * qJ(2) ^ 2 * t125, t111 ^ 2 / 0.2e1, -t111 * t110, t111 * qJD(3), -t110 * qJD(3), qJD(3) ^ 2 / 0.2e1, t128 * qJD(3) + t114 * t110, -t131 * qJD(3) + t114 * t111, t104 ^ 2 / 0.2e1, -t104 * t103, t104 * t119, -t103 * t119, t119 ^ 2 / 0.2e1, t105 * t103 + t127 * t119, t105 * t104 - t132 * t119, t98 * t103 - t96 * t119, -t97 * t103 + t96 * t104, -t98 * t104 + t97 * t119, t97 ^ 2 / 0.2e1 + t98 ^ 2 / 0.2e1 + t96 ^ 2 / 0.2e1;];
T_reg = t1;
