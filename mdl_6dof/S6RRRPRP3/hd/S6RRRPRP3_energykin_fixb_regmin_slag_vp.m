% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:41:30
% EndTime: 2019-03-09 16:41:30
% DurationCPUTime: 0.14s
% Computational Cost: add. (640->51), mult. (1356->106), div. (0->0), fcn. (969->8), ass. (0->44)
t197 = -pkin(8) - pkin(7);
t185 = qJD(1) ^ 2;
t196 = t185 / 0.2e1;
t184 = cos(qJ(2));
t195 = t184 * t185;
t180 = sin(qJ(3));
t183 = cos(qJ(3));
t191 = qJD(1) * t184;
t181 = sin(qJ(2));
t192 = qJD(1) * t181;
t167 = t180 * t192 - t183 * t191;
t168 = (t180 * t184 + t181 * t183) * qJD(1);
t173 = (-pkin(2) * t184 - pkin(1)) * qJD(1);
t158 = t167 * pkin(3) - t168 * qJ(4) + t173;
t176 = qJD(2) + qJD(3);
t171 = qJD(2) * pkin(2) + t197 * t192;
t172 = t197 * t191;
t193 = t180 * t171 - t183 * t172;
t161 = t176 * qJ(4) + t193;
t177 = sin(pkin(10));
t178 = cos(pkin(10));
t151 = t178 * t158 - t177 * t161;
t164 = t178 * t168 + t177 * t176;
t148 = t167 * pkin(4) - t164 * pkin(9) + t151;
t152 = t177 * t158 + t178 * t161;
t163 = t177 * t168 - t178 * t176;
t150 = -t163 * pkin(9) + t152;
t179 = sin(qJ(5));
t182 = cos(qJ(5));
t194 = t179 * t148 + t182 * t150;
t190 = qJD(1) * qJD(2);
t189 = t181 * t190;
t188 = t184 * t190;
t187 = t183 * t171 + t180 * t172;
t186 = t182 * t148 - t179 * t150;
t160 = -t176 * pkin(3) + qJD(4) - t187;
t153 = t163 * pkin(4) + t160;
t166 = qJD(5) + t167;
t155 = -t179 * t163 + t182 * t164;
t154 = t182 * t163 + t179 * t164;
t146 = t154 * pkin(5) - t155 * qJ(6) + t153;
t145 = t166 * qJ(6) + t194;
t144 = -t166 * pkin(5) + qJD(6) - t186;
t1 = [t196, 0, 0, t181 ^ 2 * t196, t181 * t195, t189, t188, qJD(2) ^ 2 / 0.2e1, pkin(1) * t195 - pkin(7) * t189, -t185 * pkin(1) * t181 - pkin(7) * t188, t168 ^ 2 / 0.2e1, -t168 * t167, t168 * t176, -t167 * t176, t176 ^ 2 / 0.2e1, t173 * t167 + t187 * t176, t173 * t168 - t193 * t176, t151 * t167 + t160 * t163, -t152 * t167 + t160 * t164, -t151 * t164 - t152 * t163, t152 ^ 2 / 0.2e1 + t151 ^ 2 / 0.2e1 + t160 ^ 2 / 0.2e1, t155 ^ 2 / 0.2e1, -t155 * t154, t155 * t166, -t154 * t166, t166 ^ 2 / 0.2e1, t153 * t154 + t186 * t166, t153 * t155 - t194 * t166, -t144 * t166 + t146 * t154, t144 * t155 - t145 * t154, t145 * t166 - t146 * t155, t145 ^ 2 / 0.2e1 + t146 ^ 2 / 0.2e1 + t144 ^ 2 / 0.2e1;];
T_reg  = t1;
