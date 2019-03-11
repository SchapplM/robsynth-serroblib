% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRRP3
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
% Datum: 2019-03-10 01:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:09:29
% EndTime: 2019-03-10 01:09:29
% DurationCPUTime: 0.14s
% Computational Cost: add. (473->48), mult. (1002->101), div. (0->0), fcn. (731->8), ass. (0->45)
t195 = -pkin(8) - pkin(7);
t180 = qJD(1) ^ 2;
t194 = t180 / 0.2e1;
t193 = cos(qJ(5));
t179 = cos(qJ(2));
t192 = t179 * t180;
t175 = sin(qJ(3));
t176 = sin(qJ(2));
t178 = cos(qJ(3));
t164 = (t175 * t179 + t176 * t178) * qJD(1);
t172 = qJD(2) + qJD(3);
t174 = sin(qJ(4));
t177 = cos(qJ(4));
t159 = t177 * t164 + t174 * t172;
t187 = qJD(1) * t179;
t188 = qJD(1) * t176;
t163 = t175 * t188 - t178 * t187;
t162 = qJD(4) + t163;
t169 = (-pkin(2) * t179 - pkin(1)) * qJD(1);
t153 = t163 * pkin(3) - t164 * pkin(9) + t169;
t167 = qJD(2) * pkin(2) + t195 * t188;
t168 = t195 * t187;
t189 = t175 * t167 - t178 * t168;
t156 = t172 * pkin(9) + t189;
t182 = t177 * t153 - t174 * t156;
t144 = t162 * pkin(4) - t159 * pkin(10) + t182;
t158 = t174 * t164 - t177 * t172;
t190 = t174 * t153 + t177 * t156;
t147 = -t158 * pkin(10) + t190;
t173 = sin(qJ(5));
t191 = t173 * t144 + t193 * t147;
t186 = qJD(1) * qJD(2);
t185 = t176 * t186;
t184 = t179 * t186;
t183 = t193 * t144 - t173 * t147;
t181 = t178 * t167 + t175 * t168;
t155 = -t172 * pkin(3) - t181;
t148 = t158 * pkin(4) + t155;
t160 = qJD(5) + t162;
t150 = -t173 * t158 + t193 * t159;
t149 = t193 * t158 + t173 * t159;
t145 = t149 * pkin(5) + qJD(6) + t148;
t141 = -t149 * qJ(6) + t191;
t140 = t160 * pkin(5) - t150 * qJ(6) + t183;
t1 = [t194, 0, 0, t176 ^ 2 * t194, t176 * t192, t185, t184, qJD(2) ^ 2 / 0.2e1, pkin(1) * t192 - pkin(7) * t185, -t180 * pkin(1) * t176 - pkin(7) * t184, t164 ^ 2 / 0.2e1, -t164 * t163, t164 * t172, -t163 * t172, t172 ^ 2 / 0.2e1, t169 * t163 + t181 * t172, t169 * t164 - t189 * t172, t159 ^ 2 / 0.2e1, -t159 * t158, t159 * t162, -t158 * t162, t162 ^ 2 / 0.2e1, t155 * t158 + t182 * t162, t155 * t159 - t190 * t162, t150 ^ 2 / 0.2e1, -t150 * t149, t150 * t160, -t149 * t160, t160 ^ 2 / 0.2e1, t148 * t149 + t183 * t160, t148 * t150 - t191 * t160, -t140 * t150 - t141 * t149, t141 ^ 2 / 0.2e1 + t140 ^ 2 / 0.2e1 + t145 ^ 2 / 0.2e1;];
T_reg  = t1;
