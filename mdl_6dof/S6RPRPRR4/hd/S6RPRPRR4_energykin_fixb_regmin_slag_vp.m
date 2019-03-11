% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:46:26
% EndTime: 2019-03-09 03:46:26
% DurationCPUTime: 0.12s
% Computational Cost: add. (192->45), mult. (417->99), div. (0->0), fcn. (235->8), ass. (0->39)
t182 = -pkin(3) - pkin(8);
t169 = qJD(1) ^ 2;
t181 = t169 / 0.2e1;
t161 = sin(pkin(10));
t154 = (pkin(1) * t161 + pkin(7)) * qJD(1);
t165 = sin(qJ(3));
t168 = cos(qJ(3));
t173 = t168 * qJD(2) - t165 * t154;
t172 = qJD(4) - t173;
t177 = t165 * qJD(1);
t142 = pkin(4) * t177 + t182 * qJD(3) + t172;
t162 = cos(pkin(10));
t175 = -pkin(1) * t162 - pkin(2);
t171 = -qJ(4) * t165 + t175;
t145 = (t182 * t168 + t171) * qJD(1);
t164 = sin(qJ(5));
t167 = cos(qJ(5));
t180 = t164 * t142 + t167 * t145;
t179 = t165 * qJD(2) + t168 * t154;
t178 = qJD(1) * t168;
t176 = qJD(1) * qJD(3);
t147 = -qJD(3) * qJ(4) - t179;
t174 = t167 * t142 - t164 * t145;
t157 = qJD(5) + t177;
t144 = pkin(4) * t178 - t147;
t166 = cos(qJ(6));
t163 = sin(qJ(6));
t156 = qJD(6) + t157;
t155 = t175 * qJD(1);
t153 = t167 * qJD(3) - t164 * t178;
t152 = t164 * qJD(3) + t167 * t178;
t148 = (-pkin(3) * t168 + t171) * qJD(1);
t146 = -qJD(3) * pkin(3) + t172;
t141 = -t163 * t152 + t166 * t153;
t140 = t166 * t152 + t163 * t153;
t137 = t152 * pkin(5) + t144;
t136 = -t152 * pkin(9) + t180;
t135 = t157 * pkin(5) - t153 * pkin(9) + t174;
t1 = [t181, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t161 ^ 2 / 0.2e1 + t162 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t169, t165 ^ 2 * t181, t168 * t169 * t165, t165 * t176, t168 * t176, qJD(3) ^ 2 / 0.2e1, t173 * qJD(3) - t155 * t178, -t179 * qJD(3) + t155 * t177 (t146 * t165 - t147 * t168) * qJD(1), t146 * qJD(3) + t148 * t178, -t147 * qJD(3) - t148 * t177, t148 ^ 2 / 0.2e1 + t147 ^ 2 / 0.2e1 + t146 ^ 2 / 0.2e1, t153 ^ 2 / 0.2e1, -t153 * t152, t153 * t157, -t152 * t157, t157 ^ 2 / 0.2e1, t144 * t152 + t174 * t157, t144 * t153 - t180 * t157, t141 ^ 2 / 0.2e1, -t141 * t140, t141 * t156, -t140 * t156, t156 ^ 2 / 0.2e1 (t166 * t135 - t163 * t136) * t156 + t137 * t140 -(t163 * t135 + t166 * t136) * t156 + t137 * t141;];
T_reg  = t1;
