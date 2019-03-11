% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPRP1
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
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:33:00
% EndTime: 2019-03-09 16:33:00
% DurationCPUTime: 0.14s
% Computational Cost: add. (441->47), mult. (1050->98), div. (0->0), fcn. (762->8), ass. (0->44)
t194 = -pkin(8) - pkin(7);
t180 = qJD(1) ^ 2;
t193 = t180 / 0.2e1;
t192 = cos(qJ(3));
t191 = cos(qJ(5));
t179 = cos(qJ(2));
t190 = t179 * t180;
t177 = sin(qJ(3));
t178 = sin(qJ(2));
t166 = (t177 * t179 + t192 * t178) * qJD(1);
t173 = qJD(2) + qJD(3);
t187 = qJD(1) * t178;
t168 = qJD(2) * pkin(2) + t194 * t187;
t186 = qJD(1) * t179;
t169 = t194 * t186;
t181 = t192 * t168 + t177 * t169;
t152 = t173 * pkin(3) - t166 * qJ(4) + t181;
t165 = t177 * t187 - t192 * t186;
t188 = t177 * t168 - t192 * t169;
t155 = -t165 * qJ(4) + t188;
t174 = sin(pkin(10));
t175 = cos(pkin(10));
t147 = t174 * t152 + t175 * t155;
t145 = t173 * pkin(9) + t147;
t159 = -t175 * t165 - t174 * t166;
t160 = -t174 * t165 + t175 * t166;
t170 = (-pkin(2) * t179 - pkin(1)) * qJD(1);
t161 = t165 * pkin(3) + qJD(4) + t170;
t150 = -t159 * pkin(4) - t160 * pkin(9) + t161;
t176 = sin(qJ(5));
t189 = t191 * t145 + t176 * t150;
t185 = qJD(1) * qJD(2);
t184 = t178 * t185;
t183 = t179 * t185;
t182 = -t176 * t145 + t191 * t150;
t146 = t175 * t152 - t174 * t155;
t144 = -t173 * pkin(4) - t146;
t158 = qJD(5) - t159;
t157 = t191 * t160 + t176 * t173;
t156 = t176 * t160 - t191 * t173;
t142 = t156 * pkin(5) + qJD(6) + t144;
t141 = -t156 * qJ(6) + t189;
t140 = t158 * pkin(5) - t157 * qJ(6) + t182;
t1 = [t193, 0, 0, t178 ^ 2 * t193, t178 * t190, t184, t183, qJD(2) ^ 2 / 0.2e1, pkin(1) * t190 - pkin(7) * t184, -t180 * pkin(1) * t178 - pkin(7) * t183, t166 ^ 2 / 0.2e1, -t166 * t165, t166 * t173, -t165 * t173, t173 ^ 2 / 0.2e1, t170 * t165 + t181 * t173, t170 * t166 - t188 * t173, -t146 * t160 + t147 * t159, t147 ^ 2 / 0.2e1 + t146 ^ 2 / 0.2e1 + t161 ^ 2 / 0.2e1, t157 ^ 2 / 0.2e1, -t157 * t156, t157 * t158, -t156 * t158, t158 ^ 2 / 0.2e1, t144 * t156 + t182 * t158, t144 * t157 - t189 * t158, -t140 * t157 - t141 * t156, t141 ^ 2 / 0.2e1 + t140 ^ 2 / 0.2e1 + t142 ^ 2 / 0.2e1;];
T_reg  = t1;
