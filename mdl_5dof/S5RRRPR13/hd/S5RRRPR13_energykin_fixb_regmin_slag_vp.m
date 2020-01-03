% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPR13_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR13_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:46:32
% EndTime: 2019-12-31 21:46:32
% DurationCPUTime: 0.12s
% Computational Cost: add. (262->42), mult. (649->89), div. (0->0), fcn. (470->8), ass. (0->39)
t188 = pkin(3) + pkin(9);
t166 = sin(pkin(5));
t174 = qJD(1) ^ 2;
t187 = t166 ^ 2 * t174;
t183 = cos(pkin(5)) * qJD(1);
t164 = qJD(2) + t183;
t170 = sin(qJ(2));
t173 = cos(qJ(2));
t184 = qJD(1) * t166;
t179 = t173 * t184;
t182 = pkin(1) * t183;
t185 = pkin(7) * t179 + t170 * t182;
t152 = t164 * pkin(8) + t185;
t154 = (-pkin(2) * t173 - pkin(8) * t170 - pkin(1)) * t184;
t169 = sin(qJ(3));
t172 = cos(qJ(3));
t186 = t172 * t152 + t169 * t154;
t181 = t173 * t187;
t180 = t170 * t184;
t178 = -t169 * t152 + t172 * t154;
t159 = -qJD(3) + t179;
t145 = t159 * qJ(4) - t186;
t177 = qJD(4) - t178;
t176 = -pkin(7) * t180 + t173 * t182;
t157 = t169 * t164 + t172 * t180;
t151 = -t164 * pkin(2) - t176;
t175 = -t157 * qJ(4) + t151;
t171 = cos(qJ(5));
t168 = sin(qJ(5));
t156 = -t172 * t164 + t169 * t180;
t155 = qJD(5) + t157;
t147 = t168 * t156 - t171 * t159;
t146 = -t171 * t156 - t168 * t159;
t144 = t159 * pkin(3) + t177;
t143 = t156 * pkin(3) + t175;
t142 = -t156 * pkin(4) - t145;
t141 = t188 * t156 + t175;
t140 = t157 * pkin(4) + t188 * t159 + t177;
t1 = [t174 / 0.2e1, 0, 0, t170 ^ 2 * t187 / 0.2e1, t170 * t181, t164 * t180, t164 * t179, t164 ^ 2 / 0.2e1, pkin(1) * t181 + t176 * t164, -pkin(1) * t170 * t187 - t185 * t164, t157 ^ 2 / 0.2e1, -t157 * t156, -t157 * t159, t156 * t159, t159 ^ 2 / 0.2e1, t151 * t156 - t178 * t159, t151 * t157 + t186 * t159, t144 * t157 + t145 * t156, -t143 * t156 - t144 * t159, -t143 * t157 + t145 * t159, t143 ^ 2 / 0.2e1 + t145 ^ 2 / 0.2e1 + t144 ^ 2 / 0.2e1, t147 ^ 2 / 0.2e1, -t147 * t146, t147 * t155, -t146 * t155, t155 ^ 2 / 0.2e1, (t171 * t140 - t168 * t141) * t155 + t142 * t146, -(t168 * t140 + t171 * t141) * t155 + t142 * t147;];
T_reg = t1;
