% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRPR12
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
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR12_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR12_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:30:14
% EndTime: 2019-12-31 18:30:14
% DurationCPUTime: 0.12s
% Computational Cost: add. (247->43), mult. (640->89), div. (0->0), fcn. (446->8), ass. (0->36)
t160 = pkin(6) + qJ(2);
t159 = cos(pkin(9));
t150 = sin(qJ(3));
t152 = cos(qJ(3));
t148 = cos(pkin(8));
t156 = qJD(1) * t148;
t147 = sin(pkin(8));
t157 = qJD(1) * t147;
t136 = t150 * t157 - t152 * t156;
t137 = (t147 * t152 + t148 * t150) * qJD(1);
t140 = qJD(2) + (-pkin(2) * t148 - pkin(1)) * qJD(1);
t125 = t136 * pkin(3) - t137 * qJ(4) + t140;
t138 = t160 * t157;
t139 = t160 * t156;
t158 = -t150 * t138 + t152 * t139;
t128 = qJD(3) * qJ(4) + t158;
t146 = sin(pkin(9));
t119 = t146 * t125 + t159 * t128;
t118 = t159 * t125 - t128 * t146;
t155 = -t138 * t152 - t150 * t139;
t127 = -qJD(3) * pkin(3) + qJD(4) - t155;
t153 = qJD(1) ^ 2;
t151 = cos(qJ(5));
t149 = sin(qJ(5));
t145 = t148 ^ 2;
t144 = t147 ^ 2;
t142 = -qJD(1) * pkin(1) + qJD(2);
t132 = qJD(5) + t136;
t131 = t146 * qJD(3) + t137 * t159;
t130 = -qJD(3) * t159 + t146 * t137;
t122 = -t130 * t149 + t131 * t151;
t121 = t151 * t130 + t131 * t149;
t120 = pkin(4) * t130 + t127;
t117 = -pkin(7) * t130 + t119;
t116 = pkin(4) * t136 - pkin(7) * t131 + t118;
t1 = [t153 / 0.2e1, 0, 0, -t142 * t156, t142 * t157, (t144 + t145) * t153 * qJ(2), t142 ^ 2 / 0.2e1 + (t145 / 0.2e1 + t144 / 0.2e1) * qJ(2) ^ 2 * t153, t137 ^ 2 / 0.2e1, -t137 * t136, t137 * qJD(3), -t136 * qJD(3), qJD(3) ^ 2 / 0.2e1, qJD(3) * t155 + t140 * t136, -qJD(3) * t158 + t140 * t137, t118 * t136 + t127 * t130, -t119 * t136 + t127 * t131, -t118 * t131 - t119 * t130, t119 ^ 2 / 0.2e1 + t118 ^ 2 / 0.2e1 + t127 ^ 2 / 0.2e1, t122 ^ 2 / 0.2e1, -t122 * t121, t122 * t132, -t121 * t132, t132 ^ 2 / 0.2e1, (t116 * t151 - t117 * t149) * t132 + t120 * t121, -(t116 * t149 + t117 * t151) * t132 + t120 * t122;];
T_reg = t1;
