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
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:31:21
% EndTime: 2020-01-03 11:31:21
% DurationCPUTime: 0.16s
% Computational Cost: add. (184->43), mult. (523->95), div. (0->0), fcn. (355->8), ass. (0->37)
t152 = sin(pkin(8));
t153 = cos(pkin(9));
t167 = t152 * t153;
t154 = cos(pkin(8));
t138 = qJD(2) + (-pkin(2) * t154 - qJ(3) * t152 - pkin(1)) * qJD(1);
t137 = t153 * t138;
t151 = sin(pkin(9));
t127 = t137 + (-pkin(6) * t167 + (-qJ(2) * t151 - pkin(3)) * t154) * qJD(1);
t164 = t154 * qJD(1);
t163 = qJ(2) * t164;
t132 = t151 * t138 + t153 * t163;
t165 = qJD(1) * t152;
t162 = t151 * t165;
t130 = -pkin(6) * t162 + t132;
t156 = sin(qJ(4));
t158 = cos(qJ(4));
t166 = t156 * t127 + t158 * t130;
t144 = qJ(2) * t165 + qJD(3);
t139 = pkin(3) * t162 + t144;
t161 = t158 * t127 - t156 * t130;
t145 = -qJD(4) + t164;
t159 = qJD(1) ^ 2;
t157 = cos(qJ(5));
t155 = sin(qJ(5));
t150 = t154 ^ 2;
t149 = t152 ^ 2;
t148 = -qJD(1) * pkin(1) + qJD(2);
t142 = -qJD(5) + t145;
t135 = (-t151 * t156 + t153 * t158) * t165;
t134 = (t151 * t158 + t153 * t156) * t165;
t131 = -t151 * t163 + t137;
t129 = t134 * pkin(4) + t139;
t124 = -t155 * t134 + t157 * t135;
t123 = t157 * t134 + t155 * t135;
t122 = -t134 * pkin(7) + t166;
t121 = -t145 * pkin(4) - t135 * pkin(7) + t161;
t1 = [t159 / 0.2e1, 0, 0, -t148 * t164, t148 * t165, (t149 + t150) * t159 * qJ(2), t148 ^ 2 / 0.2e1 + (t150 / 0.2e1 + t149 / 0.2e1) * qJ(2) ^ 2 * t159, (t144 * t151 * t152 - t131 * t154) * qJD(1), (t132 * t154 + t144 * t167) * qJD(1), (-t131 * t153 - t132 * t151) * t165, t132 ^ 2 / 0.2e1 + t131 ^ 2 / 0.2e1 + t144 ^ 2 / 0.2e1, t135 ^ 2 / 0.2e1, -t135 * t134, -t135 * t145, t134 * t145, t145 ^ 2 / 0.2e1, t139 * t134 - t161 * t145, t139 * t135 + t166 * t145, t124 ^ 2 / 0.2e1, -t124 * t123, -t124 * t142, t123 * t142, t142 ^ 2 / 0.2e1, -(t157 * t121 - t155 * t122) * t142 + t129 * t123, (t155 * t121 + t157 * t122) * t142 + t129 * t124;];
T_reg = t1;
