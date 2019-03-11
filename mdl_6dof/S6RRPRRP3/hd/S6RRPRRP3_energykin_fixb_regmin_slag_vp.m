% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:50:33
% EndTime: 2019-03-09 11:50:33
% DurationCPUTime: 0.13s
% Computational Cost: add. (417->47), mult. (998->98), div. (0->0), fcn. (722->8), ass. (0->44)
t188 = qJD(1) ^ 2;
t201 = t188 / 0.2e1;
t200 = cos(qJ(5));
t199 = pkin(7) + qJ(3);
t187 = cos(qJ(2));
t198 = t187 * t188;
t181 = sin(pkin(10));
t182 = cos(pkin(10));
t185 = sin(qJ(2));
t173 = (t181 * t187 + t182 * t185) * qJD(1);
t184 = sin(qJ(4));
t186 = cos(qJ(4));
t168 = t184 * qJD(2) + t186 * t173;
t194 = qJD(1) * t187;
t195 = qJD(1) * t185;
t172 = -t181 * t195 + t182 * t194;
t171 = qJD(4) - t172;
t178 = qJD(3) + (-pkin(2) * t187 - pkin(1)) * qJD(1);
t160 = -t172 * pkin(3) - t173 * pkin(8) + t178;
t176 = qJD(2) * pkin(2) - t199 * t195;
t177 = t199 * t194;
t165 = t181 * t176 + t182 * t177;
t163 = qJD(2) * pkin(8) + t165;
t189 = t186 * t160 - t184 * t163;
t151 = t171 * pkin(4) - t168 * pkin(9) + t189;
t167 = -t186 * qJD(2) + t184 * t173;
t196 = t184 * t160 + t186 * t163;
t154 = -t167 * pkin(9) + t196;
t183 = sin(qJ(5));
t197 = t183 * t151 + t200 * t154;
t193 = qJD(1) * qJD(2);
t192 = t185 * t193;
t191 = t187 * t193;
t190 = t200 * t151 - t183 * t154;
t164 = t182 * t176 - t181 * t177;
t162 = -qJD(2) * pkin(3) - t164;
t155 = t167 * pkin(4) + t162;
t169 = qJD(5) + t171;
t157 = -t183 * t167 + t200 * t168;
t156 = t200 * t167 + t183 * t168;
t152 = t156 * pkin(5) + qJD(6) + t155;
t148 = -t156 * qJ(6) + t197;
t147 = t169 * pkin(5) - t157 * qJ(6) + t190;
t1 = [t201, 0, 0, t185 ^ 2 * t201, t185 * t198, t192, t191, qJD(2) ^ 2 / 0.2e1, pkin(1) * t198 - pkin(7) * t192, -t188 * pkin(1) * t185 - pkin(7) * t191, -t164 * t173 + t165 * t172, t165 ^ 2 / 0.2e1 + t164 ^ 2 / 0.2e1 + t178 ^ 2 / 0.2e1, t168 ^ 2 / 0.2e1, -t168 * t167, t168 * t171, -t167 * t171, t171 ^ 2 / 0.2e1, t162 * t167 + t189 * t171, t162 * t168 - t196 * t171, t157 ^ 2 / 0.2e1, -t157 * t156, t157 * t169, -t156 * t169, t169 ^ 2 / 0.2e1, t155 * t156 + t190 * t169, t155 * t157 - t197 * t169, -t147 * t157 - t148 * t156, t148 ^ 2 / 0.2e1 + t147 ^ 2 / 0.2e1 + t152 ^ 2 / 0.2e1;];
T_reg  = t1;
