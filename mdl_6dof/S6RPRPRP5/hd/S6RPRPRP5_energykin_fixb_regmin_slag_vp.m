% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRPRP5
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
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRP5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:16:18
% EndTime: 2019-03-09 03:16:18
% DurationCPUTime: 0.13s
% Computational Cost: add. (543->53), mult. (1328->105), div. (0->0), fcn. (966->8), ass. (0->41)
t194 = pkin(7) + qJ(2);
t183 = sin(qJ(3));
t185 = cos(qJ(3));
t181 = cos(pkin(9));
t190 = qJD(1) * t181;
t179 = sin(pkin(9));
t191 = qJD(1) * t179;
t168 = t183 * t191 - t185 * t190;
t169 = (t179 * t185 + t181 * t183) * qJD(1);
t172 = qJD(2) + (-pkin(2) * t181 - pkin(1)) * qJD(1);
t157 = t168 * pkin(3) - t169 * qJ(4) + t172;
t170 = t194 * t191;
t171 = t194 * t190;
t192 = -t183 * t170 + t185 * t171;
t160 = qJD(3) * qJ(4) + t192;
t178 = sin(pkin(10));
t180 = cos(pkin(10));
t150 = t180 * t157 - t178 * t160;
t163 = t178 * qJD(3) + t180 * t169;
t147 = t168 * pkin(4) - t163 * pkin(8) + t150;
t151 = t178 * t157 + t180 * t160;
t162 = -t180 * qJD(3) + t178 * t169;
t149 = -t162 * pkin(8) + t151;
t182 = sin(qJ(5));
t184 = cos(qJ(5));
t193 = t182 * t147 + t184 * t149;
t189 = -t185 * t170 - t183 * t171;
t188 = t184 * t147 - t182 * t149;
t159 = -qJD(3) * pkin(3) + qJD(4) - t189;
t152 = t162 * pkin(4) + t159;
t186 = qJD(1) ^ 2;
t177 = t181 ^ 2;
t176 = t179 ^ 2;
t174 = -qJD(1) * pkin(1) + qJD(2);
t164 = qJD(5) + t168;
t154 = -t182 * t162 + t184 * t163;
t153 = t184 * t162 + t182 * t163;
t145 = t153 * pkin(5) - t154 * qJ(6) + t152;
t144 = t164 * qJ(6) + t193;
t143 = -t164 * pkin(5) + qJD(6) - t188;
t1 = [t186 / 0.2e1, 0, 0, -t174 * t190, t174 * t191 (t176 + t177) * t186 * qJ(2), t174 ^ 2 / 0.2e1 + (t177 / 0.2e1 + t176 / 0.2e1) * qJ(2) ^ 2 * t186, t169 ^ 2 / 0.2e1, -t169 * t168, t169 * qJD(3), -t168 * qJD(3), qJD(3) ^ 2 / 0.2e1, t189 * qJD(3) + t172 * t168, -t192 * qJD(3) + t172 * t169, t150 * t168 + t159 * t162, -t151 * t168 + t159 * t163, -t150 * t163 - t151 * t162, t151 ^ 2 / 0.2e1 + t150 ^ 2 / 0.2e1 + t159 ^ 2 / 0.2e1, t154 ^ 2 / 0.2e1, -t154 * t153, t154 * t164, -t153 * t164, t164 ^ 2 / 0.2e1, t152 * t153 + t188 * t164, t152 * t154 - t193 * t164, -t143 * t164 + t145 * t153, t143 * t154 - t144 * t153, t144 * t164 - t145 * t154, t144 ^ 2 / 0.2e1 + t145 ^ 2 / 0.2e1 + t143 ^ 2 / 0.2e1;];
T_reg  = t1;
