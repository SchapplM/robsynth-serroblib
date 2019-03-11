% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPPR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:42:38
% EndTime: 2019-03-09 02:42:38
% DurationCPUTime: 0.09s
% Computational Cost: add. (231->45), mult. (532->96), div. (0->0), fcn. (323->8), ass. (0->38)
t177 = pkin(4) + pkin(8);
t167 = qJD(1) ^ 2;
t176 = t167 / 0.2e1;
t160 = sin(pkin(9));
t154 = (pkin(1) * t160 + pkin(7)) * qJD(1);
t166 = cos(qJ(3));
t158 = t166 * qJD(2);
t164 = sin(qJ(3));
t144 = qJD(3) * pkin(3) + t158 + (-qJ(4) * qJD(1) - t154) * t164;
t173 = qJD(1) * t166;
t175 = t164 * qJD(2) + t166 * t154;
t147 = qJ(4) * t173 + t175;
t159 = sin(pkin(10));
t161 = cos(pkin(10));
t139 = t159 * t144 + t161 * t147;
t174 = qJD(1) * t164;
t172 = qJD(1) * qJD(3);
t162 = cos(pkin(9));
t171 = -pkin(1) * t162 - pkin(2);
t138 = t161 * t144 - t159 * t147;
t137 = -qJD(3) * qJ(5) - t139;
t170 = qJD(5) - t138;
t152 = (t159 * t166 + t161 * t164) * qJD(1);
t150 = qJD(4) + (-pkin(3) * t166 + t171) * qJD(1);
t169 = -t152 * qJ(5) + t150;
t165 = cos(qJ(6));
t163 = sin(qJ(6));
t155 = t171 * qJD(1);
t151 = t159 * t174 - t161 * t173;
t149 = qJD(6) + t152;
t146 = t165 * qJD(3) + t163 * t151;
t145 = t163 * qJD(3) - t165 * t151;
t140 = t151 * pkin(4) + t169;
t136 = -qJD(3) * pkin(4) + t170;
t135 = t177 * t151 + t169;
t134 = -t151 * pkin(5) - t137;
t133 = t152 * pkin(5) - t177 * qJD(3) + t170;
t1 = [t176, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t160 ^ 2 / 0.2e1 + t162 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t167, t164 ^ 2 * t176, t164 * t167 * t166, t164 * t172, t166 * t172, qJD(3) ^ 2 / 0.2e1, -t155 * t173 + (-t164 * t154 + t158) * qJD(3), -t175 * qJD(3) + t155 * t174, -t138 * t152 - t139 * t151, t139 ^ 2 / 0.2e1 + t138 ^ 2 / 0.2e1 + t150 ^ 2 / 0.2e1, t136 * t152 + t137 * t151, t136 * qJD(3) - t140 * t151, -t137 * qJD(3) - t140 * t152, t140 ^ 2 / 0.2e1 + t137 ^ 2 / 0.2e1 + t136 ^ 2 / 0.2e1, t146 ^ 2 / 0.2e1, -t146 * t145, t146 * t149, -t145 * t149, t149 ^ 2 / 0.2e1 (t165 * t133 - t163 * t135) * t149 + t134 * t145 -(t163 * t133 + t165 * t135) * t149 + t134 * t146;];
T_reg  = t1;
