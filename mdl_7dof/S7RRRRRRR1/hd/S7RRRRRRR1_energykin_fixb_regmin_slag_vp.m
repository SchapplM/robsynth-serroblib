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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-26 21:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
% StartTime: 2018-11-26 20:34:13
% EndTime: 2018-11-26 20:34:13
% DurationCPUTime: 0.18s
% Computational Cost: add. (598->43), mult. (1207->106), div. (0->0), fcn. (1040->12), ass. (0->44)
t177 = sin(qJ(3));
t182 = cos(qJ(3));
t178 = sin(qJ(2));
t188 = qJD(1) * t178;
t170 = -t177 * qJD(2) + t182 * t188;
t184 = qJD(1) ^ 2;
t193 = t184 / 0.2e1;
t192 = cos(qJ(6));
t169 = -t182 * qJD(2) - t177 * t188;
t166 = t169 * pkin(2);
t176 = sin(qJ(4));
t191 = t166 * t176;
t181 = cos(qJ(4));
t190 = t181 * t166;
t183 = cos(qJ(2));
t171 = t183 * qJD(1) + qJD(3);
t162 = t181 * t170 - t176 * t171;
t165 = t170 * pkin(2);
t157 = -t162 * pkin(3) - t165;
t167 = -qJD(4) - t169;
t158 = -t167 * pkin(3) + t190;
t175 = sin(qJ(5));
t180 = cos(qJ(5));
t151 = t175 * t157 + t180 * t158;
t174 = sin(qJ(6));
t189 = t192 * t151 + t174 * t191;
t186 = qJD(1) * qJD(2);
t156 = t180 * t162 - t175 * t167;
t161 = t176 * t170 + t181 * t171;
t160 = qJD(5) + t161;
t152 = t174 * t156 - t192 * t160;
t155 = t175 * t162 + t180 * t167;
t150 = -t180 * t157 + t175 * t158;
t179 = cos(qJ(7));
t173 = sin(qJ(7));
t154 = qJD(6) + t155;
t153 = t192 * t156 + t174 * t160;
t148 = -qJD(7) + t152;
t146 = -t174 * t151 + t192 * t191;
t145 = t179 * t153 - t173 * t154;
t144 = t173 * t153 + t179 * t154;
t143 = -t154 * pkin(4) + t189;
t142 = t153 * pkin(4) + t150;
t1 = [t193, 0, 0, t178 ^ 2 * t193, t178 * t184 * t183, t178 * t186, t183 * t186, qJD(2) ^ 2 / 0.2e1, 0, 0, t170 ^ 2 / 0.2e1, t170 * t169, t170 * t171, t169 * t171, t171 ^ 2 / 0.2e1, -t165 * t171, -t166 * t171, t162 ^ 2 / 0.2e1, -t162 * t161, -t162 * t167, t161 * t167, t167 ^ 2 / 0.2e1, -t165 * t161 + t167 * t191, -t165 * t162 + t167 * t190, t156 ^ 2 / 0.2e1, -t156 * t155, t156 * t160, -t155 * t160, t160 ^ 2 / 0.2e1, -t150 * t160 + t155 * t191, -t151 * t160 + t156 * t191, t153 ^ 2 / 0.2e1, -t153 * t152, t153 * t154, -t152 * t154, t154 ^ 2 / 0.2e1, t146 * t154 + t150 * t152, t150 * t153 - t189 * t154, t145 ^ 2 / 0.2e1, -t145 * t144, -t145 * t148, t144 * t148, t148 ^ 2 / 0.2e1 -(-t179 * t142 - t173 * t143) * t148 + t146 * t144 (-t173 * t142 + t179 * t143) * t148 + t146 * t145;];
T_reg  = t1;
