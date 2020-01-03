% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR7_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:04:08
% EndTime: 2019-12-31 19:04:09
% DurationCPUTime: 0.07s
% Computational Cost: add. (129->34), mult. (306->81), div. (0->0), fcn. (188->8), ass. (0->33)
t139 = qJD(1) ^ 2;
t150 = t139 / 0.2e1;
t149 = cos(qJ(4));
t132 = sin(pkin(9));
t124 = (pkin(1) * t132 + pkin(6)) * qJD(1);
t136 = sin(qJ(3));
t138 = cos(qJ(3));
t147 = t136 * qJD(2) + t138 * t124;
t117 = qJD(3) * pkin(7) + t147;
t133 = cos(pkin(9));
t143 = -pkin(1) * t133 - pkin(2);
t118 = (-pkin(3) * t138 - pkin(7) * t136 + t143) * qJD(1);
t135 = sin(qJ(4));
t148 = t149 * t117 + t135 * t118;
t146 = qJD(1) * t136;
t145 = t138 * qJD(1);
t144 = qJD(1) * qJD(3);
t142 = -t135 * t117 + t149 * t118;
t141 = t138 * qJD(2) - t136 * t124;
t128 = -qJD(4) + t145;
t116 = -qJD(3) * pkin(3) - t141;
t137 = cos(qJ(5));
t134 = sin(qJ(5));
t126 = -qJD(5) + t128;
t125 = t143 * qJD(1);
t123 = t135 * qJD(3) + t149 * t146;
t122 = -t149 * qJD(3) + t135 * t146;
t112 = -t134 * t122 + t137 * t123;
t111 = t137 * t122 + t134 * t123;
t110 = t122 * pkin(4) + t116;
t109 = -t122 * pkin(8) + t148;
t108 = -t128 * pkin(4) - t123 * pkin(8) + t142;
t1 = [t150, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t132 ^ 2 / 0.2e1 + t133 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t139, t136 ^ 2 * t150, t136 * t139 * t138, t136 * t144, t138 * t144, qJD(3) ^ 2 / 0.2e1, t141 * qJD(3) - t125 * t145, -t147 * qJD(3) + t125 * t146, t123 ^ 2 / 0.2e1, -t123 * t122, -t123 * t128, t122 * t128, t128 ^ 2 / 0.2e1, t116 * t122 - t142 * t128, t116 * t123 + t148 * t128, t112 ^ 2 / 0.2e1, -t112 * t111, -t112 * t126, t111 * t126, t126 ^ 2 / 0.2e1, -(t137 * t108 - t134 * t109) * t126 + t110 * t111, (t134 * t108 + t137 * t109) * t126 + t110 * t112;];
T_reg = t1;
