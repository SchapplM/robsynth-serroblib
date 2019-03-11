% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRRP8
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
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRP8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP8_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:25:15
% EndTime: 2019-03-09 12:25:15
% DurationCPUTime: 0.13s
% Computational Cost: add. (624->52), mult. (1394->107), div. (0->0), fcn. (998->8), ass. (0->43)
t189 = qJD(1) ^ 2;
t202 = t189 / 0.2e1;
t201 = cos(qJ(4));
t200 = cos(pkin(10));
t188 = cos(qJ(2));
t199 = t188 * t189;
t183 = sin(pkin(10));
t186 = sin(qJ(2));
t196 = qJD(1) * t186;
t172 = -t200 * qJD(2) + t183 * t196;
t173 = t183 * qJD(2) + t200 * t196;
t185 = sin(qJ(4));
t164 = -t185 * t172 + t201 * t173;
t195 = t188 * qJD(1);
t179 = -qJD(4) + t195;
t171 = (-pkin(2) * t188 - qJ(3) * t186 - pkin(1)) * qJD(1);
t176 = pkin(7) * t195 + qJD(2) * qJ(3);
t165 = t200 * t171 - t183 * t176;
t159 = -pkin(3) * t195 - t173 * pkin(8) + t165;
t166 = t183 * t171 + t200 * t176;
t161 = -t172 * pkin(8) + t166;
t191 = t201 * t159 - t185 * t161;
t151 = -t179 * pkin(4) - t164 * pkin(9) + t191;
t163 = t201 * t172 + t185 * t173;
t197 = t185 * t159 + t201 * t161;
t153 = -t163 * pkin(9) + t197;
t184 = sin(qJ(5));
t187 = cos(qJ(5));
t198 = t184 * t151 + t187 * t153;
t194 = qJD(1) * qJD(2);
t193 = t186 * t194;
t192 = t188 * t194;
t175 = -qJD(2) * pkin(2) + pkin(7) * t196 + qJD(3);
t190 = t187 * t151 - t184 * t153;
t167 = t172 * pkin(3) + t175;
t156 = t163 * pkin(4) + t167;
t177 = -qJD(5) + t179;
t155 = -t184 * t163 + t187 * t164;
t154 = t187 * t163 + t184 * t164;
t149 = t154 * pkin(5) - t155 * qJ(6) + t156;
t148 = -t177 * qJ(6) + t198;
t147 = t177 * pkin(5) + qJD(6) - t190;
t1 = [t202, 0, 0, t186 ^ 2 * t202, t186 * t199, t193, t192, qJD(2) ^ 2 / 0.2e1, pkin(1) * t199 - pkin(7) * t193, -t189 * pkin(1) * t186 - pkin(7) * t192, -t165 * t195 + t175 * t172, t166 * t195 + t175 * t173, -t165 * t173 - t166 * t172, t166 ^ 2 / 0.2e1 + t165 ^ 2 / 0.2e1 + t175 ^ 2 / 0.2e1, t164 ^ 2 / 0.2e1, -t164 * t163, -t164 * t179, t163 * t179, t179 ^ 2 / 0.2e1, t167 * t163 - t191 * t179, t167 * t164 + t197 * t179, t155 ^ 2 / 0.2e1, -t155 * t154, -t155 * t177, t154 * t177, t177 ^ 2 / 0.2e1, t156 * t154 - t190 * t177, t156 * t155 + t198 * t177, t147 * t177 + t149 * t154, t147 * t155 - t148 * t154, -t148 * t177 - t149 * t155, t148 ^ 2 / 0.2e1 + t149 ^ 2 / 0.2e1 + t147 ^ 2 / 0.2e1;];
T_reg  = t1;
