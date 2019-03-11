% Calculate minimal parameter regressor of fixed base kinetic energy for
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% qJD [7x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% 
% Output:
% T_reg [1x45]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 08:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S7RRRRRRR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(7,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_energykin_fixb_regmin_slag_vp: qJ has to be [7x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [7 1]), ...
  'S7RRRRRRR1_energykin_fixb_regmin_slag_vp: qJD has to be [7x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_energykin_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 07:42:52
% EndTime: 2019-03-10 07:42:52
% DurationCPUTime: 0.17s
% Computational Cost: add. (598->43), mult. (1207->106), div. (0->0), fcn. (1040->12), ass. (0->44)
t177 = sin(qJ(3));
t182 = cos(qJ(3));
t178 = sin(qJ(2));
t188 = qJD(1) * t178;
t170 = -qJD(2) * t177 + t182 * t188;
t184 = qJD(1) ^ 2;
t193 = t184 / 0.2e1;
t192 = cos(qJ(6));
t169 = -qJD(2) * t182 - t177 * t188;
t166 = t169 * pkin(2);
t176 = sin(qJ(4));
t191 = t166 * t176;
t181 = cos(qJ(4));
t190 = t166 * t181;
t183 = cos(qJ(2));
t171 = t183 * qJD(1) + qJD(3);
t162 = t170 * t181 - t171 * t176;
t165 = t170 * pkin(2);
t157 = -pkin(3) * t162 - t165;
t167 = -qJD(4) - t169;
t158 = -pkin(3) * t167 + t190;
t175 = sin(qJ(5));
t180 = cos(qJ(5));
t151 = t157 * t175 + t158 * t180;
t174 = sin(qJ(6));
t189 = t192 * t151 + t174 * t191;
t186 = qJD(1) * qJD(2);
t156 = t162 * t180 - t167 * t175;
t161 = t170 * t176 + t171 * t181;
t160 = qJD(5) + t161;
t152 = t156 * t174 - t192 * t160;
t155 = t162 * t175 + t180 * t167;
t150 = -t157 * t180 + t175 * t158;
t179 = cos(qJ(7));
t173 = sin(qJ(7));
t154 = qJD(6) + t155;
t153 = t192 * t156 + t174 * t160;
t148 = -qJD(7) + t152;
t146 = -t174 * t151 + t192 * t191;
t145 = t153 * t179 - t154 * t173;
t144 = t153 * t173 + t154 * t179;
t143 = -pkin(4) * t154 + t189;
t142 = pkin(4) * t153 + t150;
t1 = [t193, 0, 0, t178 ^ 2 * t193, t178 * t184 * t183, t178 * t186, t183 * t186, qJD(2) ^ 2 / 0.2e1, 0, 0, t170 ^ 2 / 0.2e1, t170 * t169, t170 * t171, t169 * t171, t171 ^ 2 / 0.2e1, -t165 * t171, -t166 * t171, t162 ^ 2 / 0.2e1, -t162 * t161, -t162 * t167, t161 * t167, t167 ^ 2 / 0.2e1, -t161 * t165 + t167 * t191, -t162 * t165 + t167 * t190, t156 ^ 2 / 0.2e1, -t156 * t155, t156 * t160, -t155 * t160, t160 ^ 2 / 0.2e1, -t150 * t160 + t155 * t191, -t151 * t160 + t156 * t191, t153 ^ 2 / 0.2e1, -t153 * t152, t153 * t154, -t152 * t154, t154 ^ 2 / 0.2e1, t146 * t154 + t150 * t152, t150 * t153 - t189 * t154, t145 ^ 2 / 0.2e1, -t145 * t144, -t145 * t148, t144 * t148, t148 ^ 2 / 0.2e1 -(-t142 * t179 - t143 * t173) * t148 + t146 * t144 (-t173 * t142 + t143 * t179) * t148 + t146 * t145;];
T_reg  = t1;
