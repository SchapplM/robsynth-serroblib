% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPPR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR6_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:54:37
% EndTime: 2019-03-09 02:54:37
% DurationCPUTime: 0.13s
% Computational Cost: add. (346->46), mult. (700->95), div. (0->0), fcn. (440->8), ass. (0->37)
t173 = qJD(1) ^ 2;
t179 = t173 / 0.2e1;
t178 = cos(pkin(10));
t177 = t173 * qJ(2);
t172 = cos(qJ(3));
t161 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t174 = -qJ(4) * qJD(1) + t161;
t155 = qJD(3) * pkin(3) + t174 * t172;
t170 = sin(qJ(3));
t156 = t174 * t170;
t167 = sin(pkin(9));
t168 = cos(pkin(9));
t147 = t167 * t155 + t168 * t156;
t143 = qJD(3) * qJ(5) + t147;
t158 = (t167 * t172 + t168 * t170) * qJD(1);
t159 = (-t167 * t170 + t168 * t172) * qJD(1);
t160 = qJD(4) + (pkin(3) * t170 + qJ(2)) * qJD(1);
t148 = t158 * pkin(4) - t159 * qJ(5) + t160;
t166 = sin(pkin(10));
t137 = t178 * t143 + t166 * t148;
t176 = qJD(3) * t161;
t175 = qJD(1) * qJD(3);
t136 = -t166 * t143 + t178 * t148;
t146 = t168 * t155 - t167 * t156;
t142 = -qJD(3) * pkin(4) + qJD(5) - t146;
t171 = cos(qJ(6));
t169 = sin(qJ(6));
t163 = -qJD(1) * pkin(1) + qJD(2);
t157 = qJD(6) + t158;
t151 = t166 * qJD(3) + t178 * t159;
t150 = -t178 * qJD(3) + t166 * t159;
t140 = -t169 * t150 + t171 * t151;
t139 = t171 * t150 + t169 * t151;
t138 = t150 * pkin(5) + t142;
t135 = -t150 * pkin(8) + t137;
t134 = t158 * pkin(5) - t151 * pkin(8) + t136;
t1 = [t179, 0, 0, t163 * qJD(1), t177, qJ(2) ^ 2 * t179 + t163 ^ 2 / 0.2e1, t172 ^ 2 * t179, -t172 * t173 * t170, t172 * t175, -t170 * t175, qJD(3) ^ 2 / 0.2e1, t170 * t177 + t172 * t176, -t170 * t176 + t172 * t177, -t146 * t159 - t147 * t158, t147 ^ 2 / 0.2e1 + t146 ^ 2 / 0.2e1 + t160 ^ 2 / 0.2e1, t136 * t158 + t142 * t150, -t137 * t158 + t142 * t151, -t136 * t151 - t137 * t150, t137 ^ 2 / 0.2e1 + t136 ^ 2 / 0.2e1 + t142 ^ 2 / 0.2e1, t140 ^ 2 / 0.2e1, -t140 * t139, t140 * t157, -t139 * t157, t157 ^ 2 / 0.2e1 (t171 * t134 - t169 * t135) * t157 + t138 * t139 -(t169 * t134 + t171 * t135) * t157 + t138 * t140;];
T_reg  = t1;
