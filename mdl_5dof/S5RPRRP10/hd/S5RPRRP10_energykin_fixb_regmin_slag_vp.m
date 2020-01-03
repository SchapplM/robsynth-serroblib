% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRP10
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
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP10_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP10_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:52:03
% EndTime: 2019-12-31 18:52:03
% DurationCPUTime: 0.08s
% Computational Cost: add. (159->37), mult. (423->77), div. (0->0), fcn. (280->6), ass. (0->32)
t140 = cos(qJ(4));
t139 = pkin(6) + qJ(2);
t129 = sin(qJ(3));
t130 = cos(qJ(3));
t127 = cos(pkin(8));
t135 = qJD(1) * t127;
t126 = sin(pkin(8));
t136 = qJD(1) * t126;
t116 = t129 * t136 - t130 * t135;
t117 = (t126 * t130 + t127 * t129) * qJD(1);
t120 = qJD(2) + (-pkin(2) * t127 - pkin(1)) * qJD(1);
t106 = t116 * pkin(3) - t117 * pkin(7) + t120;
t118 = t139 * t136;
t119 = t139 * t135;
t137 = -t129 * t118 + t130 * t119;
t109 = qJD(3) * pkin(7) + t137;
t128 = sin(qJ(4));
t138 = t128 * t106 + t140 * t109;
t134 = t140 * t106 - t128 * t109;
t133 = -t130 * t118 - t129 * t119;
t108 = -qJD(3) * pkin(3) - t133;
t131 = qJD(1) ^ 2;
t125 = t127 ^ 2;
t124 = t126 ^ 2;
t122 = -qJD(1) * pkin(1) + qJD(2);
t112 = qJD(4) + t116;
t111 = t128 * qJD(3) + t140 * t117;
t110 = -t140 * qJD(3) + t128 * t117;
t103 = t110 * pkin(4) + qJD(5) + t108;
t102 = -t110 * qJ(5) + t138;
t101 = t112 * pkin(4) - t111 * qJ(5) + t134;
t1 = [t131 / 0.2e1, 0, 0, -t122 * t135, t122 * t136, (t124 + t125) * t131 * qJ(2), t122 ^ 2 / 0.2e1 + (t125 / 0.2e1 + t124 / 0.2e1) * qJ(2) ^ 2 * t131, t117 ^ 2 / 0.2e1, -t117 * t116, t117 * qJD(3), -t116 * qJD(3), qJD(3) ^ 2 / 0.2e1, t133 * qJD(3) + t120 * t116, -t137 * qJD(3) + t120 * t117, t111 ^ 2 / 0.2e1, -t111 * t110, t111 * t112, -t110 * t112, t112 ^ 2 / 0.2e1, t108 * t110 + t134 * t112, t108 * t111 - t138 * t112, -t101 * t111 - t102 * t110, t102 ^ 2 / 0.2e1 + t101 ^ 2 / 0.2e1 + t103 ^ 2 / 0.2e1;];
T_reg = t1;
