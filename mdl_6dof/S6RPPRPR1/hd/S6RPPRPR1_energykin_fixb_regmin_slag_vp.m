% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRPR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:40:01
% EndTime: 2019-03-09 01:40:02
% DurationCPUTime: 0.10s
% Computational Cost: add. (340->51), mult. (794->104), div. (0->0), fcn. (541->10), ass. (0->40)
t180 = cos(pkin(11));
t166 = sin(pkin(9));
t159 = (pkin(1) * t166 + qJ(3)) * qJD(1);
t167 = cos(pkin(10));
t163 = t167 * qJD(2);
t165 = sin(pkin(10));
t149 = t163 + (-pkin(7) * qJD(1) - t159) * t165;
t152 = t165 * qJD(2) + t167 * t159;
t177 = qJD(1) * t167;
t150 = pkin(7) * t177 + t152;
t170 = sin(qJ(4));
t172 = cos(qJ(4));
t179 = t170 * t149 + t172 * t150;
t137 = qJD(4) * qJ(5) + t179;
t168 = cos(pkin(9));
t176 = -pkin(1) * t168 - pkin(2);
t154 = qJD(3) + (-pkin(3) * t167 + t176) * qJD(1);
t178 = qJD(1) * t165;
t155 = t170 * t178 - t172 * t177;
t156 = (t165 * t172 + t167 * t170) * qJD(1);
t142 = t155 * pkin(4) - t156 * qJ(5) + t154;
t164 = sin(pkin(11));
t133 = t180 * t137 + t164 * t142;
t132 = -t164 * t137 + t180 * t142;
t175 = t172 * t149 - t170 * t150;
t136 = -qJD(4) * pkin(4) + qJD(5) - t175;
t173 = qJD(1) ^ 2;
t171 = cos(qJ(6));
t169 = sin(qJ(6));
t158 = t176 * qJD(1) + qJD(3);
t153 = qJD(6) + t155;
t151 = -t165 * t159 + t163;
t148 = t164 * qJD(4) + t180 * t156;
t147 = -t180 * qJD(4) + t164 * t156;
t139 = -t169 * t147 + t171 * t148;
t138 = t171 * t147 + t169 * t148;
t134 = t147 * pkin(5) + t136;
t131 = -t147 * pkin(8) + t133;
t130 = t155 * pkin(5) - t148 * pkin(8) + t132;
t1 = [t173 / 0.2e1, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t166 ^ 2 / 0.2e1 + t168 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t173, -t158 * t177, t158 * t178 (-t151 * t165 + t152 * t167) * qJD(1), t152 ^ 2 / 0.2e1 + t151 ^ 2 / 0.2e1 + t158 ^ 2 / 0.2e1, t156 ^ 2 / 0.2e1, -t156 * t155, t156 * qJD(4), -t155 * qJD(4), qJD(4) ^ 2 / 0.2e1, t175 * qJD(4) + t154 * t155, -t179 * qJD(4) + t154 * t156, t132 * t155 + t136 * t147, -t133 * t155 + t136 * t148, -t132 * t148 - t133 * t147, t133 ^ 2 / 0.2e1 + t132 ^ 2 / 0.2e1 + t136 ^ 2 / 0.2e1, t139 ^ 2 / 0.2e1, -t139 * t138, t139 * t153, -t138 * t153, t153 ^ 2 / 0.2e1 (t171 * t130 - t169 * t131) * t153 + t134 * t138 -(t169 * t130 + t171 * t131) * t153 + t134 * t139;];
T_reg  = t1;
