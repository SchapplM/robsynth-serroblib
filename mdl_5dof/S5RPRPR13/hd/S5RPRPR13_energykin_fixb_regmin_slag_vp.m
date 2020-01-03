% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR13_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR13_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:32:52
% EndTime: 2019-12-31 18:32:52
% DurationCPUTime: 0.12s
% Computational Cost: add. (152->40), mult. (408->81), div. (0->0), fcn. (259->6), ass. (0->33)
t142 = pkin(3) + pkin(7);
t141 = pkin(6) + qJ(2);
t127 = sin(pkin(8));
t139 = qJD(1) * t127;
t120 = t141 * t139;
t128 = cos(pkin(8));
t138 = qJD(1) * t128;
t121 = t141 * t138;
t130 = sin(qJ(3));
t132 = cos(qJ(3));
t140 = -t130 * t120 + t132 * t121;
t137 = -t132 * t120 - t130 * t121;
t110 = -qJD(3) * qJ(4) - t140;
t136 = qJD(4) - t137;
t122 = qJD(2) + (-pkin(2) * t128 - pkin(1)) * qJD(1);
t119 = (t127 * t132 + t128 * t130) * qJD(1);
t135 = -t119 * qJ(4) + t122;
t133 = qJD(1) ^ 2;
t131 = cos(qJ(5));
t129 = sin(qJ(5));
t126 = t128 ^ 2;
t125 = t127 ^ 2;
t124 = -qJD(1) * pkin(1) + qJD(2);
t118 = t130 * t139 - t132 * t138;
t114 = qJD(5) + t119;
t112 = t131 * qJD(3) + t129 * t118;
t111 = t129 * qJD(3) - t131 * t118;
t109 = -qJD(3) * pkin(3) + t136;
t108 = t118 * pkin(3) + t135;
t107 = -t118 * pkin(4) - t110;
t106 = t119 * pkin(4) - t142 * qJD(3) + t136;
t105 = t142 * t118 + t135;
t1 = [t133 / 0.2e1, 0, 0, -t124 * t138, t124 * t139, (t125 + t126) * t133 * qJ(2), t124 ^ 2 / 0.2e1 + (t126 / 0.2e1 + t125 / 0.2e1) * qJ(2) ^ 2 * t133, t119 ^ 2 / 0.2e1, -t119 * t118, t119 * qJD(3), -t118 * qJD(3), qJD(3) ^ 2 / 0.2e1, t137 * qJD(3) + t122 * t118, -t140 * qJD(3) + t122 * t119, t109 * t119 + t110 * t118, t109 * qJD(3) - t108 * t118, -t110 * qJD(3) - t108 * t119, t108 ^ 2 / 0.2e1 + t110 ^ 2 / 0.2e1 + t109 ^ 2 / 0.2e1, t112 ^ 2 / 0.2e1, -t112 * t111, t112 * t114, -t111 * t114, t114 ^ 2 / 0.2e1, (-t129 * t105 + t131 * t106) * t114 + t107 * t111, -(t131 * t105 + t129 * t106) * t114 + t107 * t112;];
T_reg = t1;
