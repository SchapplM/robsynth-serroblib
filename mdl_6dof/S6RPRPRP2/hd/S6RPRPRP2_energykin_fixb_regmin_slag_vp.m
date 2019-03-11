% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:06:12
% EndTime: 2019-03-09 03:06:12
% DurationCPUTime: 0.09s
% Computational Cost: add. (313->44), mult. (699->96), div. (0->0), fcn. (447->8), ass. (0->37)
t174 = qJD(1) ^ 2;
t183 = t174 / 0.2e1;
t167 = sin(pkin(9));
t160 = (pkin(1) * t167 + pkin(7)) * qJD(1);
t173 = cos(qJ(3));
t165 = t173 * qJD(2);
t171 = sin(qJ(3));
t151 = qJD(3) * pkin(3) + t165 + (-qJ(4) * qJD(1) - t160) * t171;
t179 = qJD(1) * t173;
t181 = t171 * qJD(2) + t173 * t160;
t154 = qJ(4) * t179 + t181;
t166 = sin(pkin(10));
t168 = cos(pkin(10));
t145 = t166 * t151 + t168 * t154;
t143 = qJD(3) * pkin(8) + t145;
t169 = cos(pkin(9));
t177 = -pkin(1) * t169 - pkin(2);
t156 = qJD(4) + (-pkin(3) * t173 + t177) * qJD(1);
t180 = qJD(1) * t171;
t157 = -t166 * t180 + t168 * t179;
t158 = (t166 * t173 + t168 * t171) * qJD(1);
t147 = -t157 * pkin(4) - t158 * pkin(8) + t156;
t170 = sin(qJ(5));
t172 = cos(qJ(5));
t182 = t172 * t143 + t170 * t147;
t178 = qJD(1) * qJD(3);
t144 = t168 * t151 - t166 * t154;
t176 = -t170 * t143 + t172 * t147;
t142 = -qJD(3) * pkin(4) - t144;
t161 = t177 * qJD(1);
t155 = qJD(5) - t157;
t153 = t170 * qJD(3) + t172 * t158;
t152 = -t172 * qJD(3) + t170 * t158;
t140 = t152 * pkin(5) - t153 * qJ(6) + t142;
t139 = t155 * qJ(6) + t182;
t138 = -t155 * pkin(5) + qJD(6) - t176;
t1 = [t183, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t167 ^ 2 / 0.2e1 + t169 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t174, t171 ^ 2 * t183, t171 * t174 * t173, t171 * t178, t173 * t178, qJD(3) ^ 2 / 0.2e1, -t161 * t179 + (-t171 * t160 + t165) * qJD(3), -t181 * qJD(3) + t161 * t180, -t144 * t158 + t145 * t157, t145 ^ 2 / 0.2e1 + t144 ^ 2 / 0.2e1 + t156 ^ 2 / 0.2e1, t153 ^ 2 / 0.2e1, -t153 * t152, t153 * t155, -t152 * t155, t155 ^ 2 / 0.2e1, t142 * t152 + t176 * t155, t142 * t153 - t182 * t155, -t138 * t155 + t140 * t152, t138 * t153 - t139 * t152, t139 * t155 - t140 * t153, t139 ^ 2 / 0.2e1 + t140 ^ 2 / 0.2e1 + t138 ^ 2 / 0.2e1;];
T_reg  = t1;
