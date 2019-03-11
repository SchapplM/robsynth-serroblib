% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:08:44
% EndTime: 2019-03-09 06:08:44
% DurationCPUTime: 0.13s
% Computational Cost: add. (401->50), mult. (1021->100), div. (0->0), fcn. (770->8), ass. (0->42)
t187 = cos(qJ(3));
t186 = cos(qJ(5));
t185 = pkin(7) + qJ(2);
t168 = qJD(3) + qJD(4);
t169 = sin(pkin(10));
t170 = cos(pkin(10));
t173 = sin(qJ(3));
t159 = (t187 * t169 + t170 * t173) * qJD(1);
t181 = qJD(1) * t169;
t160 = t185 * t181;
t180 = qJD(1) * t170;
t161 = t185 * t180;
t177 = -t187 * t160 - t173 * t161;
t146 = qJD(3) * pkin(3) - t159 * pkin(8) + t177;
t158 = t173 * t181 - t187 * t180;
t182 = -t173 * t160 + t187 * t161;
t147 = -t158 * pkin(8) + t182;
t172 = sin(qJ(4));
t174 = cos(qJ(4));
t183 = t172 * t146 + t174 * t147;
t139 = t168 * pkin(9) + t183;
t151 = t174 * t158 + t172 * t159;
t152 = -t172 * t158 + t174 * t159;
t162 = qJD(2) + (-pkin(2) * t170 - pkin(1)) * qJD(1);
t153 = t158 * pkin(3) + t162;
t142 = t151 * pkin(4) - t152 * pkin(9) + t153;
t171 = sin(qJ(5));
t184 = t186 * t139 + t171 * t142;
t179 = -t171 * t139 + t186 * t142;
t178 = t174 * t146 - t172 * t147;
t138 = -t168 * pkin(4) - t178;
t175 = qJD(1) ^ 2;
t167 = t170 ^ 2;
t166 = t169 ^ 2;
t165 = -qJD(1) * pkin(1) + qJD(2);
t150 = qJD(5) + t151;
t149 = t186 * t152 + t171 * t168;
t148 = t171 * t152 - t186 * t168;
t136 = t148 * pkin(5) + qJD(6) + t138;
t135 = -t148 * qJ(6) + t184;
t134 = t150 * pkin(5) - t149 * qJ(6) + t179;
t1 = [t175 / 0.2e1, 0, 0, -t165 * t180, t165 * t181 (t166 + t167) * t175 * qJ(2), t165 ^ 2 / 0.2e1 + (t167 / 0.2e1 + t166 / 0.2e1) * qJ(2) ^ 2 * t175, t159 ^ 2 / 0.2e1, -t158 * t159, qJD(3) * t159, -t158 * qJD(3), qJD(3) ^ 2 / 0.2e1, t177 * qJD(3) + t162 * t158, -t182 * qJD(3) + t162 * t159, t152 ^ 2 / 0.2e1, -t152 * t151, t152 * t168, -t151 * t168, t168 ^ 2 / 0.2e1, t153 * t151 + t178 * t168, t153 * t152 - t183 * t168, t149 ^ 2 / 0.2e1, -t149 * t148, t149 * t150, -t148 * t150, t150 ^ 2 / 0.2e1, t138 * t148 + t179 * t150, t138 * t149 - t184 * t150, -t134 * t149 - t135 * t148, t135 ^ 2 / 0.2e1 + t134 ^ 2 / 0.2e1 + t136 ^ 2 / 0.2e1;];
T_reg  = t1;
