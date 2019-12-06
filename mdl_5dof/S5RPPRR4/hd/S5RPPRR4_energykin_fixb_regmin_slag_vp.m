% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:44:54
% EndTime: 2019-12-05 17:44:54
% DurationCPUTime: 0.17s
% Computational Cost: add. (184->43), mult. (523->95), div. (0->0), fcn. (355->8), ass. (0->37)
t153 = sin(pkin(8));
t154 = cos(pkin(9));
t168 = t153 * t154;
t155 = cos(pkin(8));
t139 = qJD(2) + (-pkin(2) * t155 - qJ(3) * t153 - pkin(1)) * qJD(1);
t138 = t154 * t139;
t152 = sin(pkin(9));
t128 = t138 + (-pkin(6) * t168 + (-qJ(2) * t152 - pkin(3)) * t155) * qJD(1);
t165 = t155 * qJD(1);
t164 = qJ(2) * t165;
t133 = t152 * t139 + t154 * t164;
t166 = qJD(1) * t153;
t163 = t152 * t166;
t131 = -pkin(6) * t163 + t133;
t157 = sin(qJ(4));
t159 = cos(qJ(4));
t167 = t157 * t128 + t159 * t131;
t145 = qJ(2) * t166 + qJD(3);
t140 = pkin(3) * t163 + t145;
t162 = t159 * t128 - t157 * t131;
t146 = -qJD(4) + t165;
t160 = qJD(1) ^ 2;
t158 = cos(qJ(5));
t156 = sin(qJ(5));
t151 = t155 ^ 2;
t150 = t153 ^ 2;
t149 = -qJD(1) * pkin(1) + qJD(2);
t143 = -qJD(5) + t146;
t136 = (-t152 * t157 + t154 * t159) * t166;
t135 = (t152 * t159 + t154 * t157) * t166;
t132 = -t152 * t164 + t138;
t130 = t135 * pkin(4) + t140;
t125 = -t156 * t135 + t158 * t136;
t124 = t158 * t135 + t156 * t136;
t123 = -t135 * pkin(7) + t167;
t122 = -t146 * pkin(4) - t136 * pkin(7) + t162;
t1 = [t160 / 0.2e1, 0, 0, -t149 * t165, t149 * t166, (t150 + t151) * t160 * qJ(2), t149 ^ 2 / 0.2e1 + (t151 / 0.2e1 + t150 / 0.2e1) * qJ(2) ^ 2 * t160, (t145 * t152 * t153 - t132 * t155) * qJD(1), (t133 * t155 + t145 * t168) * qJD(1), (-t132 * t154 - t133 * t152) * t166, t133 ^ 2 / 0.2e1 + t132 ^ 2 / 0.2e1 + t145 ^ 2 / 0.2e1, t136 ^ 2 / 0.2e1, -t136 * t135, -t136 * t146, t135 * t146, t146 ^ 2 / 0.2e1, t140 * t135 - t162 * t146, t140 * t136 + t167 * t146, t125 ^ 2 / 0.2e1, -t125 * t124, -t125 * t143, t124 * t143, t143 ^ 2 / 0.2e1, -(t158 * t122 - t156 * t123) * t143 + t130 * t124, (t156 * t122 + t158 * t123) * t143 + t130 * t125;];
T_reg = t1;
