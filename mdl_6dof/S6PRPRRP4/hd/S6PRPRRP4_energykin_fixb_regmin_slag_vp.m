% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRRP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:12:18
% EndTime: 2019-03-08 20:12:18
% DurationCPUTime: 0.11s
% Computational Cost: add. (298->46), mult. (714->94), div. (0->0), fcn. (536->10), ass. (0->38)
t183 = sin(qJ(2));
t194 = qJD(1) * sin(pkin(6));
t172 = qJD(2) * qJ(3) + t183 * t194;
t179 = cos(pkin(11));
t193 = qJD(1) * cos(pkin(6));
t174 = t179 * t193;
t177 = sin(pkin(11));
t160 = t174 + (-pkin(8) * qJD(2) - t172) * t177;
t165 = t179 * t172 + t177 * t193;
t191 = qJD(2) * t179;
t161 = pkin(8) * t191 + t165;
t182 = sin(qJ(4));
t185 = cos(qJ(4));
t195 = t182 * t160 + t185 * t161;
t154 = qJD(4) * pkin(9) + t195;
t186 = cos(qJ(2));
t188 = -t186 * t194 + qJD(3);
t166 = (-pkin(3) * t179 - pkin(2)) * qJD(2) + t188;
t192 = qJD(2) * t177;
t168 = t182 * t192 - t185 * t191;
t169 = (t177 * t185 + t179 * t182) * qJD(2);
t156 = t168 * pkin(4) - t169 * pkin(9) + t166;
t181 = sin(qJ(5));
t184 = cos(qJ(5));
t196 = t184 * t154 + t181 * t156;
t190 = qJD(2) * t194;
t189 = t185 * t160 - t182 * t161;
t187 = -t181 * t154 + t184 * t156;
t153 = -qJD(4) * pkin(4) - t189;
t171 = -qJD(2) * pkin(2) + t188;
t167 = qJD(5) + t168;
t164 = -t177 * t172 + t174;
t163 = t181 * qJD(4) + t184 * t169;
t162 = -t184 * qJD(4) + t181 * t169;
t151 = t162 * pkin(5) - t163 * qJ(6) + t153;
t150 = t167 * qJ(6) + t196;
t149 = -t167 * pkin(5) + qJD(6) - t187;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, t186 * t190, -t183 * t190, -t171 * t191, t171 * t192 (-t164 * t177 + t165 * t179) * qJD(2), t165 ^ 2 / 0.2e1 + t164 ^ 2 / 0.2e1 + t171 ^ 2 / 0.2e1, t169 ^ 2 / 0.2e1, -t169 * t168, t169 * qJD(4), -t168 * qJD(4), qJD(4) ^ 2 / 0.2e1, t189 * qJD(4) + t166 * t168, -t195 * qJD(4) + t166 * t169, t163 ^ 2 / 0.2e1, -t163 * t162, t163 * t167, -t162 * t167, t167 ^ 2 / 0.2e1, t153 * t162 + t187 * t167, t153 * t163 - t196 * t167, -t149 * t167 + t151 * t162, t149 * t163 - t150 * t162, t150 * t167 - t151 * t163, t150 ^ 2 / 0.2e1 + t151 ^ 2 / 0.2e1 + t149 ^ 2 / 0.2e1;];
T_reg  = t1;
