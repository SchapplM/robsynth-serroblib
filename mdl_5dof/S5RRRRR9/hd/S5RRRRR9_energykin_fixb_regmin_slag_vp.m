% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x31]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRR9_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR9_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:29:51
% EndTime: 2019-12-31 22:29:51
% DurationCPUTime: 0.14s
% Computational Cost: add. (259->41), mult. (586->90), div. (0->0), fcn. (423->8), ass. (0->39)
t165 = qJD(1) ^ 2;
t177 = t165 / 0.2e1;
t176 = cos(qJ(4));
t164 = cos(qJ(2));
t175 = t164 * t165;
t160 = sin(qJ(3));
t163 = cos(qJ(3));
t161 = sin(qJ(2));
t172 = qJD(1) * t161;
t146 = t160 * qJD(2) + t163 * t172;
t171 = t164 * qJD(1);
t154 = -qJD(3) + t171;
t144 = (-pkin(2) * t164 - pkin(7) * t161 - pkin(1)) * qJD(1);
t151 = pkin(6) * t171 + qJD(2) * pkin(7);
t166 = t163 * t144 - t160 * t151;
t134 = -t154 * pkin(3) - t146 * pkin(8) + t166;
t145 = -t163 * qJD(2) + t160 * t172;
t173 = t160 * t144 + t163 * t151;
t136 = -t145 * pkin(8) + t173;
t159 = sin(qJ(4));
t174 = t159 * t134 + t176 * t136;
t170 = qJD(1) * qJD(2);
t169 = t161 * t170;
t168 = t164 * t170;
t150 = -qJD(2) * pkin(2) + pkin(6) * t172;
t167 = t176 * t134 - t159 * t136;
t152 = -qJD(4) + t154;
t140 = t145 * pkin(3) + t150;
t162 = cos(qJ(5));
t158 = sin(qJ(5));
t148 = -qJD(5) + t152;
t139 = -t159 * t145 + t176 * t146;
t138 = t176 * t145 + t159 * t146;
t131 = t138 * pkin(4) + t140;
t130 = -t158 * t138 + t162 * t139;
t129 = t162 * t138 + t158 * t139;
t128 = -t138 * pkin(9) + t174;
t127 = -t152 * pkin(4) - t139 * pkin(9) + t167;
t1 = [t177, 0, 0, t161 ^ 2 * t177, t161 * t175, t169, t168, qJD(2) ^ 2 / 0.2e1, pkin(1) * t175 - pkin(6) * t169, -t165 * pkin(1) * t161 - pkin(6) * t168, t146 ^ 2 / 0.2e1, -t146 * t145, -t146 * t154, t145 * t154, t154 ^ 2 / 0.2e1, t150 * t145 - t166 * t154, t150 * t146 + t173 * t154, t139 ^ 2 / 0.2e1, -t139 * t138, -t152 * t139, t152 * t138, t152 ^ 2 / 0.2e1, t140 * t138 - t167 * t152, t140 * t139 + t174 * t152, t130 ^ 2 / 0.2e1, -t130 * t129, -t130 * t148, t129 * t148, t148 ^ 2 / 0.2e1, -(t162 * t127 - t158 * t128) * t148 + t131 * t129, (t158 * t127 + t162 * t128) * t148 + t131 * t130;];
T_reg = t1;
