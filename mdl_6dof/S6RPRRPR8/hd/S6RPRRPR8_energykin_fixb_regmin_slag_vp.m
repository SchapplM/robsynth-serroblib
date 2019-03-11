% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPR8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR8_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:24:50
% EndTime: 2019-03-09 05:24:50
% DurationCPUTime: 0.12s
% Computational Cost: add. (313->46), mult. (631->95), div. (0->0), fcn. (412->8), ass. (0->38)
t170 = qJD(1) ^ 2;
t178 = t170 / 0.2e1;
t177 = cos(qJ(4));
t176 = t170 * qJ(2);
t166 = sin(qJ(4));
t169 = cos(qJ(3));
t174 = qJD(1) * t169;
t157 = t166 * qJD(3) + t177 * t174;
t167 = sin(qJ(3));
t160 = t167 * qJD(1) + qJD(4);
t153 = (pkin(3) * t167 - pkin(8) * t169 + qJ(2)) * qJD(1);
t159 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t154 = qJD(3) * pkin(8) + t167 * t159;
t171 = t177 * t153 - t166 * t154;
t143 = t160 * pkin(4) - t157 * qJ(5) + t171;
t156 = -t177 * qJD(3) + t166 * t174;
t175 = t166 * t153 + t177 * t154;
t145 = -t156 * qJ(5) + t175;
t163 = sin(pkin(10));
t164 = cos(pkin(10));
t137 = t163 * t143 + t164 * t145;
t173 = qJD(3) * t159;
t172 = qJD(1) * qJD(3);
t136 = t164 * t143 - t163 * t145;
t155 = -qJD(3) * pkin(3) - t169 * t159;
t149 = t156 * pkin(4) + qJD(5) + t155;
t168 = cos(qJ(6));
t165 = sin(qJ(6));
t161 = -qJD(1) * pkin(1) + qJD(2);
t158 = qJD(6) + t160;
t148 = -t163 * t156 + t164 * t157;
t147 = -t164 * t156 - t163 * t157;
t140 = -t147 * pkin(5) + t149;
t139 = t165 * t147 + t168 * t148;
t138 = -t168 * t147 + t165 * t148;
t135 = t147 * pkin(9) + t137;
t134 = t160 * pkin(5) - t148 * pkin(9) + t136;
t1 = [t178, 0, 0, t161 * qJD(1), t176, qJ(2) ^ 2 * t178 + t161 ^ 2 / 0.2e1, t169 ^ 2 * t178, -t169 * t170 * t167, t169 * t172, -t167 * t172, qJD(3) ^ 2 / 0.2e1, t167 * t176 + t169 * t173, -t167 * t173 + t169 * t176, t157 ^ 2 / 0.2e1, -t157 * t156, t157 * t160, -t156 * t160, t160 ^ 2 / 0.2e1, t155 * t156 + t171 * t160, t155 * t157 - t175 * t160, -t136 * t148 + t137 * t147, t137 ^ 2 / 0.2e1 + t136 ^ 2 / 0.2e1 + t149 ^ 2 / 0.2e1, t139 ^ 2 / 0.2e1, -t139 * t138, t139 * t158, -t138 * t158, t158 ^ 2 / 0.2e1 (t168 * t134 - t165 * t135) * t158 + t140 * t138 -(t165 * t134 + t168 * t135) * t158 + t140 * t139;];
T_reg  = t1;
