% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRRP6
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
% Datum: 2019-03-09 06:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRP6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP6_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:16:43
% EndTime: 2019-03-09 06:16:43
% DurationCPUTime: 0.10s
% Computational Cost: add. (399->50), mult. (985->100), div. (0->0), fcn. (728->8), ass. (0->42)
t190 = cos(qJ(5));
t189 = pkin(7) + qJ(2);
t172 = sin(pkin(10));
t173 = cos(pkin(10));
t176 = sin(qJ(3));
t178 = cos(qJ(3));
t163 = (t172 * t178 + t173 * t176) * qJD(1);
t175 = sin(qJ(4));
t177 = cos(qJ(4));
t156 = t175 * qJD(3) + t177 * t163;
t184 = qJD(1) * t173;
t185 = qJD(1) * t172;
t162 = t176 * t185 - t178 * t184;
t158 = qJD(4) + t162;
t166 = qJD(2) + (-pkin(2) * t173 - pkin(1)) * qJD(1);
t150 = t162 * pkin(3) - t163 * pkin(8) + t166;
t164 = t189 * t185;
t165 = t189 * t184;
t186 = -t176 * t164 + t178 * t165;
t153 = qJD(3) * pkin(8) + t186;
t182 = t177 * t150 - t175 * t153;
t141 = t158 * pkin(4) - t156 * pkin(9) + t182;
t155 = -t177 * qJD(3) + t175 * t163;
t187 = t175 * t150 + t177 * t153;
t144 = -t155 * pkin(9) + t187;
t174 = sin(qJ(5));
t188 = t174 * t141 + t190 * t144;
t183 = t190 * t141 - t174 * t144;
t181 = -t178 * t164 - t176 * t165;
t152 = -qJD(3) * pkin(3) - t181;
t145 = t155 * pkin(4) + t152;
t179 = qJD(1) ^ 2;
t171 = t173 ^ 2;
t170 = t172 ^ 2;
t168 = -qJD(1) * pkin(1) + qJD(2);
t157 = qJD(5) + t158;
t147 = -t174 * t155 + t190 * t156;
t146 = t190 * t155 + t174 * t156;
t142 = t146 * pkin(5) + qJD(6) + t145;
t138 = -t146 * qJ(6) + t188;
t137 = t157 * pkin(5) - t147 * qJ(6) + t183;
t1 = [t179 / 0.2e1, 0, 0, -t168 * t184, t168 * t185 (t170 + t171) * t179 * qJ(2), t168 ^ 2 / 0.2e1 + (t171 / 0.2e1 + t170 / 0.2e1) * qJ(2) ^ 2 * t179, t163 ^ 2 / 0.2e1, -t163 * t162, t163 * qJD(3), -t162 * qJD(3), qJD(3) ^ 2 / 0.2e1, t181 * qJD(3) + t166 * t162, -t186 * qJD(3) + t166 * t163, t156 ^ 2 / 0.2e1, -t156 * t155, t156 * t158, -t155 * t158, t158 ^ 2 / 0.2e1, t152 * t155 + t182 * t158, t152 * t156 - t187 * t158, t147 ^ 2 / 0.2e1, -t147 * t146, t147 * t157, -t146 * t157, t157 ^ 2 / 0.2e1, t145 * t146 + t183 * t157, t145 * t147 - t188 * t157, -t137 * t147 - t138 * t146, t138 ^ 2 / 0.2e1 + t137 ^ 2 / 0.2e1 + t142 ^ 2 / 0.2e1;];
T_reg  = t1;
