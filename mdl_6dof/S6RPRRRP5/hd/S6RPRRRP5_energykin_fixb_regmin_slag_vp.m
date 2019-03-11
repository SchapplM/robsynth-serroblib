% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRRP5
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
% Datum: 2019-03-09 06:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRP5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:12:24
% EndTime: 2019-03-09 06:12:24
% DurationCPUTime: 0.14s
% Computational Cost: add. (531->52), mult. (1326->104), div. (0->0), fcn. (1008->8), ass. (0->42)
t195 = cos(qJ(3));
t194 = pkin(7) + qJ(2);
t176 = qJD(3) + qJD(4);
t177 = sin(pkin(10));
t178 = cos(pkin(10));
t181 = sin(qJ(3));
t167 = (t195 * t177 + t178 * t181) * qJD(1);
t190 = qJD(1) * t177;
t168 = t194 * t190;
t189 = qJD(1) * t178;
t169 = t194 * t189;
t187 = -t195 * t168 - t181 * t169;
t154 = qJD(3) * pkin(3) - t167 * pkin(8) + t187;
t166 = t181 * t190 - t195 * t189;
t191 = -t181 * t168 + t195 * t169;
t155 = -t166 * pkin(8) + t191;
t180 = sin(qJ(4));
t183 = cos(qJ(4));
t192 = t180 * t154 + t183 * t155;
t148 = t176 * pkin(9) + t192;
t159 = t183 * t166 + t180 * t167;
t160 = -t180 * t166 + t183 * t167;
t170 = qJD(2) + (-pkin(2) * t178 - pkin(1)) * qJD(1);
t161 = t166 * pkin(3) + t170;
t150 = t159 * pkin(4) - t160 * pkin(9) + t161;
t179 = sin(qJ(5));
t182 = cos(qJ(5));
t193 = t182 * t148 + t179 * t150;
t188 = t183 * t154 - t180 * t155;
t186 = -t179 * t148 + t182 * t150;
t147 = -t176 * pkin(4) - t188;
t184 = qJD(1) ^ 2;
t175 = t178 ^ 2;
t174 = t177 ^ 2;
t173 = -qJD(1) * pkin(1) + qJD(2);
t158 = qJD(5) + t159;
t157 = t182 * t160 + t179 * t176;
t156 = t179 * t160 - t182 * t176;
t145 = t156 * pkin(5) - t157 * qJ(6) + t147;
t144 = t158 * qJ(6) + t193;
t143 = -t158 * pkin(5) + qJD(6) - t186;
t1 = [t184 / 0.2e1, 0, 0, -t173 * t189, t173 * t190 (t174 + t175) * t184 * qJ(2), t173 ^ 2 / 0.2e1 + (t175 / 0.2e1 + t174 / 0.2e1) * qJ(2) ^ 2 * t184, t167 ^ 2 / 0.2e1, -t167 * t166, t167 * qJD(3), -t166 * qJD(3), qJD(3) ^ 2 / 0.2e1, t187 * qJD(3) + t170 * t166, -t191 * qJD(3) + t170 * t167, t160 ^ 2 / 0.2e1, -t160 * t159, t160 * t176, -t159 * t176, t176 ^ 2 / 0.2e1, t161 * t159 + t188 * t176, t161 * t160 - t192 * t176, t157 ^ 2 / 0.2e1, -t157 * t156, t157 * t158, -t156 * t158, t158 ^ 2 / 0.2e1, t147 * t156 + t186 * t158, t147 * t157 - t193 * t158, -t143 * t158 + t145 * t156, t143 * t157 - t144 * t156, t144 * t158 - t145 * t157, t144 ^ 2 / 0.2e1 + t145 ^ 2 / 0.2e1 + t143 ^ 2 / 0.2e1;];
T_reg  = t1;
