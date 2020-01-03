% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% T_reg [1x20]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:19:50
% EndTime: 2019-12-31 18:19:50
% DurationCPUTime: 0.08s
% Computational Cost: add. (128->34), mult. (319->80), div. (0->0), fcn. (193->8), ass. (0->32)
t137 = qJD(1) ^ 2;
t144 = t137 / 0.2e1;
t130 = sin(pkin(8));
t123 = (pkin(1) * t130 + pkin(6)) * qJD(1);
t136 = cos(qJ(3));
t128 = t136 * qJD(2);
t134 = sin(qJ(3));
t114 = qJD(3) * pkin(3) + t128 + (-qJ(4) * qJD(1) - t123) * t134;
t141 = qJD(1) * t136;
t143 = t134 * qJD(2) + t136 * t123;
t117 = qJ(4) * t141 + t143;
t129 = sin(pkin(9));
t131 = cos(pkin(9));
t110 = t129 * t114 + t131 * t117;
t142 = qJD(1) * t134;
t140 = qJD(1) * qJD(3);
t132 = cos(pkin(8));
t139 = -pkin(1) * t132 - pkin(2);
t120 = -t129 * t142 + t131 * t141;
t109 = t131 * t114 - t129 * t117;
t119 = qJD(4) + (-pkin(3) * t136 + t139) * qJD(1);
t135 = cos(qJ(5));
t133 = sin(qJ(5));
t124 = t139 * qJD(1);
t121 = (t129 * t136 + t131 * t134) * qJD(1);
t118 = qJD(5) - t120;
t116 = t133 * qJD(3) + t135 * t121;
t115 = -t135 * qJD(3) + t133 * t121;
t111 = -t120 * pkin(4) - t121 * pkin(7) + t119;
t108 = qJD(3) * pkin(7) + t110;
t107 = -qJD(3) * pkin(4) - t109;
t1 = [t144, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t130 ^ 2 / 0.2e1 + t132 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t137, t134 ^ 2 * t144, t134 * t137 * t136, t134 * t140, t136 * t140, qJD(3) ^ 2 / 0.2e1, -t124 * t141 + (-t134 * t123 + t128) * qJD(3), -t143 * qJD(3) + t124 * t142, -t109 * t121 + t110 * t120, t110 ^ 2 / 0.2e1 + t109 ^ 2 / 0.2e1 + t119 ^ 2 / 0.2e1, t116 ^ 2 / 0.2e1, -t116 * t115, t116 * t118, -t115 * t118, t118 ^ 2 / 0.2e1, (-t133 * t108 + t135 * t111) * t118 + t107 * t115, -(t135 * t108 + t133 * t111) * t118 + t107 * t116;];
T_reg = t1;
