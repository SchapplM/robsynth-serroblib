% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x33]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP1_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:58:13
% EndTime: 2019-03-10 00:58:13
% DurationCPUTime: 0.14s
% Computational Cost: add. (454->48), mult. (1041->101), div. (0->0), fcn. (773->8), ass. (0->45)
t190 = -pkin(8) - pkin(7);
t174 = qJD(1) ^ 2;
t189 = t174 / 0.2e1;
t188 = cos(qJ(3));
t187 = cos(qJ(5));
t173 = cos(qJ(2));
t186 = t173 * t174;
t167 = qJD(2) + qJD(3);
t166 = qJD(4) + t167;
t170 = sin(qJ(3));
t171 = sin(qJ(2));
t159 = (t170 * t173 + t188 * t171) * qJD(1);
t182 = qJD(1) * t171;
t161 = qJD(2) * pkin(2) + t190 * t182;
t181 = qJD(1) * t173;
t162 = t190 * t181;
t175 = t188 * t161 + t170 * t162;
t145 = t167 * pkin(3) - t159 * pkin(9) + t175;
t158 = t170 * t182 - t188 * t181;
t183 = t170 * t161 - t188 * t162;
t148 = -t158 * pkin(9) + t183;
t169 = sin(qJ(4));
t172 = cos(qJ(4));
t184 = t169 * t145 + t172 * t148;
t140 = t166 * pkin(10) + t184;
t152 = t172 * t158 + t169 * t159;
t153 = -t169 * t158 + t172 * t159;
t163 = (-pkin(2) * t173 - pkin(1)) * qJD(1);
t154 = t158 * pkin(3) + t163;
t143 = t152 * pkin(4) - t153 * pkin(10) + t154;
t168 = sin(qJ(5));
t185 = t187 * t140 + t168 * t143;
t180 = qJD(1) * qJD(2);
t179 = t171 * t180;
t178 = t173 * t180;
t177 = -t168 * t140 + t187 * t143;
t176 = t172 * t145 - t169 * t148;
t139 = -t166 * pkin(4) - t176;
t151 = qJD(5) + t152;
t150 = t187 * t153 + t168 * t166;
t149 = t168 * t153 - t187 * t166;
t137 = t149 * pkin(5) + qJD(6) + t139;
t136 = -t149 * qJ(6) + t185;
t135 = t151 * pkin(5) - t150 * qJ(6) + t177;
t1 = [t189, 0, 0, t171 ^ 2 * t189, t171 * t186, t179, t178, qJD(2) ^ 2 / 0.2e1, pkin(1) * t186 - pkin(7) * t179, -t174 * pkin(1) * t171 - pkin(7) * t178, t159 ^ 2 / 0.2e1, -t159 * t158, t159 * t167, -t158 * t167, t167 ^ 2 / 0.2e1, t163 * t158 + t175 * t167, t163 * t159 - t183 * t167, t153 ^ 2 / 0.2e1, -t153 * t152, t153 * t166, -t152 * t166, t166 ^ 2 / 0.2e1, t154 * t152 + t176 * t166, t154 * t153 - t184 * t166, t150 ^ 2 / 0.2e1, -t150 * t149, t150 * t151, -t149 * t151, t151 ^ 2 / 0.2e1, t139 * t149 + t177 * t151, t139 * t150 - t185 * t151, -t135 * t150 - t136 * t149, t136 ^ 2 / 0.2e1 + t135 ^ 2 / 0.2e1 + t137 ^ 2 / 0.2e1;];
T_reg  = t1;
