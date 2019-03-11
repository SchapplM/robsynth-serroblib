% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
% 
% Output:
% T_reg [1x27]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPRP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:27:46
% EndTime: 2019-03-09 08:27:47
% DurationCPUTime: 0.13s
% Computational Cost: add. (572->50), mult. (1352->103), div. (0->0), fcn. (960->8), ass. (0->43)
t192 = qJD(1) ^ 2;
t202 = t192 / 0.2e1;
t201 = pkin(7) + qJ(3);
t191 = cos(qJ(2));
t200 = t191 * t192;
t185 = sin(pkin(9));
t187 = cos(pkin(9));
t197 = qJD(1) * t191;
t189 = sin(qJ(2));
t198 = qJD(1) * t189;
t175 = t185 * t198 - t187 * t197;
t176 = (t185 * t191 + t187 * t189) * qJD(1);
t181 = qJD(3) + (-pkin(2) * t191 - pkin(1)) * qJD(1);
t164 = t175 * pkin(3) - t176 * qJ(4) + t181;
t179 = qJD(2) * pkin(2) - t201 * t198;
t180 = t201 * t197;
t169 = t185 * t179 + t187 * t180;
t167 = qJD(2) * qJ(4) + t169;
t184 = sin(pkin(10));
t186 = cos(pkin(10));
t157 = t186 * t164 - t184 * t167;
t172 = t184 * qJD(2) + t186 * t176;
t154 = t175 * pkin(4) - t172 * pkin(8) + t157;
t158 = t184 * t164 + t186 * t167;
t171 = -t186 * qJD(2) + t184 * t176;
t156 = -t171 * pkin(8) + t158;
t188 = sin(qJ(5));
t190 = cos(qJ(5));
t199 = t188 * t154 + t190 * t156;
t196 = qJD(1) * qJD(2);
t195 = t189 * t196;
t194 = t191 * t196;
t168 = t187 * t179 - t185 * t180;
t193 = t190 * t154 - t188 * t156;
t166 = -qJD(2) * pkin(3) + qJD(4) - t168;
t159 = t171 * pkin(4) + t166;
t174 = qJD(5) + t175;
t161 = -t188 * t171 + t190 * t172;
t160 = t190 * t171 + t188 * t172;
t152 = t160 * pkin(5) - t161 * qJ(6) + t159;
t151 = t174 * qJ(6) + t199;
t150 = -t174 * pkin(5) + qJD(6) - t193;
t1 = [t202, 0, 0, t189 ^ 2 * t202, t189 * t200, t195, t194, qJD(2) ^ 2 / 0.2e1, pkin(1) * t200 - pkin(7) * t195, -t192 * pkin(1) * t189 - pkin(7) * t194, -t168 * t176 - t169 * t175, t169 ^ 2 / 0.2e1 + t168 ^ 2 / 0.2e1 + t181 ^ 2 / 0.2e1, t157 * t175 + t166 * t171, -t158 * t175 + t166 * t172, -t157 * t172 - t158 * t171, t158 ^ 2 / 0.2e1 + t157 ^ 2 / 0.2e1 + t166 ^ 2 / 0.2e1, t161 ^ 2 / 0.2e1, -t161 * t160, t161 * t174, -t160 * t174, t174 ^ 2 / 0.2e1, t159 * t160 + t193 * t174, t159 * t161 - t199 * t174, -t150 * t174 + t152 * t160, t150 * t161 - t151 * t160, t151 * t174 - t152 * t161, t151 ^ 2 / 0.2e1 + t152 ^ 2 / 0.2e1 + t150 ^ 2 / 0.2e1;];
T_reg  = t1;
