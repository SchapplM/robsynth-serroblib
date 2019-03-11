% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPP5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:44:41
% EndTime: 2019-03-09 04:44:41
% DurationCPUTime: 0.11s
% Computational Cost: add. (389->50), mult. (931->97), div. (0->0), fcn. (646->6), ass. (0->38)
t180 = -pkin(4) - pkin(5);
t179 = cos(qJ(4));
t178 = pkin(7) + qJ(2);
t165 = sin(qJ(3));
t166 = cos(qJ(3));
t163 = cos(pkin(9));
t173 = qJD(1) * t163;
t162 = sin(pkin(9));
t174 = qJD(1) * t162;
t152 = t165 * t174 - t166 * t173;
t153 = (t162 * t166 + t163 * t165) * qJD(1);
t156 = qJD(2) + (-pkin(2) * t163 - pkin(1)) * qJD(1);
t139 = t152 * pkin(3) - t153 * pkin(8) + t156;
t154 = t178 * t174;
t155 = t178 * t173;
t175 = -t165 * t154 + t166 * t155;
t143 = qJD(3) * pkin(8) + t175;
t164 = sin(qJ(4));
t177 = t164 * t139 + t179 * t143;
t176 = -t166 * t154 - t165 * t155;
t147 = qJD(4) + t152;
t136 = t147 * qJ(5) + t177;
t172 = qJD(3) * pkin(3) + t176;
t171 = t179 * t139 - t164 * t143;
t170 = qJD(5) - t171;
t145 = t164 * qJD(3) + t179 * t153;
t169 = t145 * qJ(5) + t172;
t167 = qJD(1) ^ 2;
t161 = t163 ^ 2;
t160 = t162 ^ 2;
t158 = -qJD(1) * pkin(1) + qJD(2);
t144 = -t179 * qJD(3) + t164 * t153;
t137 = t144 * pkin(4) - t169;
t135 = -t147 * pkin(4) + t170;
t134 = t180 * t144 + qJD(6) + t169;
t133 = t144 * qJ(6) + t136;
t132 = -t145 * qJ(6) + t180 * t147 + t170;
t1 = [t167 / 0.2e1, 0, 0, -t158 * t173, t158 * t174 (t160 + t161) * t167 * qJ(2), t158 ^ 2 / 0.2e1 + (t161 / 0.2e1 + t160 / 0.2e1) * qJ(2) ^ 2 * t167, t153 ^ 2 / 0.2e1, -t153 * t152, t153 * qJD(3), -t152 * qJD(3), qJD(3) ^ 2 / 0.2e1, t176 * qJD(3) + t156 * t152, -t175 * qJD(3) + t156 * t153, t145 ^ 2 / 0.2e1, -t145 * t144, t145 * t147, -t144 * t147, t147 ^ 2 / 0.2e1, -t144 * t172 + t171 * t147, -t145 * t172 - t177 * t147, -t135 * t147 + t137 * t144, t135 * t145 - t136 * t144, t136 * t147 - t137 * t145, t136 ^ 2 / 0.2e1 + t137 ^ 2 / 0.2e1 + t135 ^ 2 / 0.2e1, -t132 * t147 - t134 * t144, t133 * t147 + t134 * t145, -t132 * t145 + t133 * t144, t133 ^ 2 / 0.2e1 + t132 ^ 2 / 0.2e1 + t134 ^ 2 / 0.2e1;];
T_reg  = t1;
