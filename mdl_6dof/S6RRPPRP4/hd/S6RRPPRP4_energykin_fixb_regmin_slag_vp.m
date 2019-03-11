% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPRP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:39:40
% EndTime: 2019-03-09 08:39:41
% DurationCPUTime: 0.12s
% Computational Cost: add. (363->49), mult. (800->100), div. (0->0), fcn. (494->6), ass. (0->36)
t185 = qJD(1) ^ 2;
t194 = t185 / 0.2e1;
t184 = cos(qJ(2));
t193 = t184 * t185;
t182 = sin(qJ(2));
t166 = (-pkin(2) * t184 - qJ(3) * t182 - pkin(1)) * qJD(1);
t190 = t184 * qJD(1);
t173 = pkin(7) * t190 + qJD(2) * qJ(3);
t178 = sin(pkin(9));
t179 = cos(pkin(9));
t161 = t179 * t166 - t178 * t173;
t157 = pkin(3) * t190 + qJD(4) - t161;
t191 = qJD(1) * t182;
t169 = t178 * qJD(2) + t179 * t191;
t152 = pkin(4) * t190 - t169 * pkin(8) + t157;
t162 = t178 * t166 + t179 * t173;
t158 = -qJ(4) * t190 + t162;
t168 = -t179 * qJD(2) + t178 * t191;
t155 = t168 * pkin(8) + t158;
t181 = sin(qJ(5));
t183 = cos(qJ(5));
t192 = t181 * t152 + t183 * t155;
t189 = qJD(1) * qJD(2);
t172 = -qJD(2) * pkin(2) + pkin(7) * t191 + qJD(3);
t188 = t182 * t189;
t187 = t184 * t189;
t156 = t168 * pkin(3) - t169 * qJ(4) + t172;
t186 = t183 * t152 - t181 * t155;
t154 = -t168 * pkin(4) - t156;
t174 = qJD(5) + t190;
t160 = t181 * t168 + t183 * t169;
t159 = -t183 * t168 + t181 * t169;
t150 = t159 * pkin(5) - t160 * qJ(6) + t154;
t149 = t174 * qJ(6) + t192;
t148 = -t174 * pkin(5) + qJD(6) - t186;
t1 = [t194, 0, 0, t182 ^ 2 * t194, t182 * t193, t188, t187, qJD(2) ^ 2 / 0.2e1, pkin(1) * t193 - pkin(7) * t188, -t185 * pkin(1) * t182 - pkin(7) * t187, -t161 * t190 + t172 * t168, t162 * t190 + t172 * t169, -t161 * t169 - t162 * t168, t162 ^ 2 / 0.2e1 + t161 ^ 2 / 0.2e1 + t172 ^ 2 / 0.2e1, t156 * t168 + t157 * t190, t157 * t169 - t158 * t168, -t156 * t169 - t158 * t190, t158 ^ 2 / 0.2e1 + t156 ^ 2 / 0.2e1 + t157 ^ 2 / 0.2e1, t160 ^ 2 / 0.2e1, -t160 * t159, t160 * t174, -t159 * t174, t174 ^ 2 / 0.2e1, t154 * t159 + t186 * t174, t154 * t160 - t192 * t174, -t148 * t174 + t150 * t159, t148 * t160 - t149 * t159, t149 * t174 - t150 * t160, t149 ^ 2 / 0.2e1 + t150 ^ 2 / 0.2e1 + t148 ^ 2 / 0.2e1;];
T_reg  = t1;
