% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR9_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR9_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:21:57
% EndTime: 2019-12-31 20:21:57
% DurationCPUTime: 0.13s
% Computational Cost: add. (221->39), mult. (558->86), div. (0->0), fcn. (396->8), ass. (0->39)
t162 = qJD(1) ^ 2;
t173 = t162 / 0.2e1;
t172 = cos(qJ(4));
t171 = pkin(6) + qJ(3);
t161 = cos(qJ(2));
t170 = t161 * t162;
t155 = sin(pkin(9));
t156 = cos(pkin(9));
t167 = qJD(1) * t161;
t159 = sin(qJ(2));
t168 = qJD(1) * t159;
t146 = -t155 * t168 + t156 * t167;
t147 = (t155 * t161 + t156 * t159) * qJD(1);
t152 = qJD(3) + (-pkin(2) * t161 - pkin(1)) * qJD(1);
t134 = -t146 * pkin(3) - t147 * pkin(7) + t152;
t150 = qJD(2) * pkin(2) - t171 * t168;
t151 = t171 * t167;
t139 = t155 * t150 + t156 * t151;
t137 = qJD(2) * pkin(7) + t139;
t158 = sin(qJ(4));
t169 = t158 * t134 + t172 * t137;
t166 = qJD(1) * qJD(2);
t165 = t159 * t166;
t164 = t161 * t166;
t163 = t172 * t134 - t137 * t158;
t138 = t150 * t156 - t155 * t151;
t145 = qJD(4) - t146;
t136 = -qJD(2) * pkin(3) - t138;
t160 = cos(qJ(5));
t157 = sin(qJ(5));
t143 = qJD(5) + t145;
t142 = t158 * qJD(2) + t172 * t147;
t141 = -t172 * qJD(2) + t147 * t158;
t131 = -t141 * t157 + t142 * t160;
t130 = t160 * t141 + t142 * t157;
t129 = pkin(4) * t141 + t136;
t128 = -pkin(8) * t141 + t169;
t127 = pkin(4) * t145 - pkin(8) * t142 + t163;
t1 = [t173, 0, 0, t159 ^ 2 * t173, t159 * t170, t165, t164, qJD(2) ^ 2 / 0.2e1, pkin(1) * t170 - pkin(6) * t165, -pkin(1) * t159 * t162 - pkin(6) * t164, -t138 * t147 + t139 * t146, t139 ^ 2 / 0.2e1 + t138 ^ 2 / 0.2e1 + t152 ^ 2 / 0.2e1, t142 ^ 2 / 0.2e1, -t142 * t141, t142 * t145, -t141 * t145, t145 ^ 2 / 0.2e1, t136 * t141 + t163 * t145, t136 * t142 - t169 * t145, t131 ^ 2 / 0.2e1, -t131 * t130, t131 * t143, -t130 * t143, t143 ^ 2 / 0.2e1, (t127 * t160 - t128 * t157) * t143 + t129 * t130, -(t127 * t157 + t128 * t160) * t143 + t129 * t131;];
T_reg = t1;
