% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRRP7
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
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRP7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP7_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:20:59
% EndTime: 2019-03-09 06:20:59
% DurationCPUTime: 0.13s
% Computational Cost: add. (519->52), mult. (1242->104), div. (0->0), fcn. (922->8), ass. (0->42)
t196 = cos(qJ(4));
t195 = pkin(7) + qJ(2);
t178 = sin(pkin(10));
t179 = cos(pkin(10));
t182 = sin(qJ(3));
t184 = cos(qJ(3));
t169 = (t178 * t184 + t179 * t182) * qJD(1);
t181 = sin(qJ(4));
t162 = t181 * qJD(3) + t196 * t169;
t190 = qJD(1) * t179;
t191 = qJD(1) * t178;
t168 = t182 * t191 - t184 * t190;
t164 = qJD(4) + t168;
t172 = qJD(2) + (-pkin(2) * t179 - pkin(1)) * qJD(1);
t156 = t168 * pkin(3) - t169 * pkin(8) + t172;
t170 = t195 * t191;
t171 = t195 * t190;
t192 = -t182 * t170 + t184 * t171;
t159 = qJD(3) * pkin(8) + t192;
t189 = t196 * t156 - t181 * t159;
t148 = t164 * pkin(4) - t162 * pkin(9) + t189;
t161 = -t196 * qJD(3) + t181 * t169;
t193 = t181 * t156 + t196 * t159;
t150 = -t161 * pkin(9) + t193;
t180 = sin(qJ(5));
t183 = cos(qJ(5));
t194 = t180 * t148 + t183 * t150;
t188 = -t184 * t170 - t182 * t171;
t187 = t183 * t148 - t180 * t150;
t158 = -qJD(3) * pkin(3) - t188;
t151 = t161 * pkin(4) + t158;
t185 = qJD(1) ^ 2;
t177 = t179 ^ 2;
t176 = t178 ^ 2;
t174 = -qJD(1) * pkin(1) + qJD(2);
t163 = qJD(5) + t164;
t153 = -t180 * t161 + t183 * t162;
t152 = t183 * t161 + t180 * t162;
t146 = t152 * pkin(5) - t153 * qJ(6) + t151;
t145 = t163 * qJ(6) + t194;
t144 = -t163 * pkin(5) + qJD(6) - t187;
t1 = [t185 / 0.2e1, 0, 0, -t174 * t190, t174 * t191 (t176 + t177) * t185 * qJ(2), t174 ^ 2 / 0.2e1 + (t177 / 0.2e1 + t176 / 0.2e1) * qJ(2) ^ 2 * t185, t169 ^ 2 / 0.2e1, -t169 * t168, t169 * qJD(3), -t168 * qJD(3), qJD(3) ^ 2 / 0.2e1, t188 * qJD(3) + t172 * t168, -t192 * qJD(3) + t172 * t169, t162 ^ 2 / 0.2e1, -t162 * t161, t162 * t164, -t161 * t164, t164 ^ 2 / 0.2e1, t158 * t161 + t189 * t164, t158 * t162 - t193 * t164, t153 ^ 2 / 0.2e1, -t153 * t152, t153 * t163, -t152 * t163, t163 ^ 2 / 0.2e1, t151 * t152 + t187 * t163, t151 * t153 - t194 * t163, -t144 * t163 + t146 * t152, t144 * t153 - t145 * t152, t145 * t163 - t146 * t153, t145 ^ 2 / 0.2e1 + t146 ^ 2 / 0.2e1 + t144 ^ 2 / 0.2e1;];
T_reg  = t1;
