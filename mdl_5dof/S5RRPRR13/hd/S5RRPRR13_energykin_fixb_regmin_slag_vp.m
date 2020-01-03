% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRR13
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
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR13_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR13_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:33:51
% EndTime: 2019-12-31 20:33:51
% DurationCPUTime: 0.14s
% Computational Cost: add. (260->42), mult. (618->91), div. (0->0), fcn. (432->8), ass. (0->38)
t163 = qJD(1) ^ 2;
t173 = t163 / 0.2e1;
t172 = cos(qJ(4));
t162 = cos(qJ(2));
t171 = t162 * t163;
t160 = sin(qJ(2));
t144 = (-pkin(2) * t162 - qJ(3) * t160 - pkin(1)) * qJD(1);
t168 = t162 * qJD(1);
t149 = pkin(6) * t168 + qJD(2) * qJ(3);
t156 = sin(pkin(9));
t157 = cos(pkin(9));
t138 = t157 * t144 - t156 * t149;
t169 = qJD(1) * t160;
t146 = t156 * qJD(2) + t157 * t169;
t132 = -pkin(3) * t168 - t146 * pkin(7) + t138;
t139 = t156 * t144 + t157 * t149;
t145 = -t157 * qJD(2) + t156 * t169;
t134 = -t145 * pkin(7) + t139;
t159 = sin(qJ(4));
t170 = t159 * t132 + t172 * t134;
t167 = qJD(1) * qJD(2);
t166 = t160 * t167;
t165 = t162 * t167;
t164 = t172 * t132 - t159 * t134;
t152 = -qJD(4) + t168;
t148 = -qJD(2) * pkin(2) + pkin(6) * t169 + qJD(3);
t140 = t145 * pkin(3) + t148;
t161 = cos(qJ(5));
t158 = sin(qJ(5));
t150 = -qJD(5) + t152;
t137 = -t159 * t145 + t172 * t146;
t136 = t172 * t145 + t159 * t146;
t129 = t136 * pkin(4) + t140;
t128 = -t158 * t136 + t161 * t137;
t127 = t161 * t136 + t158 * t137;
t126 = -t136 * pkin(8) + t170;
t125 = -t152 * pkin(4) - t137 * pkin(8) + t164;
t1 = [t173, 0, 0, t160 ^ 2 * t173, t160 * t171, t166, t165, qJD(2) ^ 2 / 0.2e1, pkin(1) * t171 - pkin(6) * t166, -t163 * pkin(1) * t160 - pkin(6) * t165, -t138 * t168 + t148 * t145, t139 * t168 + t148 * t146, -t138 * t146 - t139 * t145, t139 ^ 2 / 0.2e1 + t138 ^ 2 / 0.2e1 + t148 ^ 2 / 0.2e1, t137 ^ 2 / 0.2e1, -t137 * t136, -t137 * t152, t136 * t152, t152 ^ 2 / 0.2e1, t140 * t136 - t164 * t152, t140 * t137 + t170 * t152, t128 ^ 2 / 0.2e1, -t128 * t127, -t128 * t150, t127 * t150, t150 ^ 2 / 0.2e1, -(t161 * t125 - t158 * t126) * t150 + t129 * t127, (t158 * t125 + t161 * t126) * t150 + t129 * t128;];
T_reg = t1;
