% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRP11_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP11_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:18:31
% EndTime: 2019-12-31 22:18:31
% DurationCPUTime: 0.10s
% Computational Cost: add. (361->41), mult. (860->89), div. (0->0), fcn. (641->8), ass. (0->38)
t171 = sin(pkin(5));
t179 = qJD(1) ^ 2;
t192 = t171 ^ 2 * t179;
t187 = cos(pkin(5)) * qJD(1);
t169 = qJD(2) + t187;
t178 = cos(qJ(2));
t175 = sin(qJ(2));
t188 = qJD(1) * t171;
t184 = t175 * t188;
t186 = pkin(1) * t187;
t180 = -pkin(7) * t184 + t178 * t186;
t156 = -t169 * pkin(2) - t180;
t174 = sin(qJ(3));
t177 = cos(qJ(3));
t160 = -t177 * t169 + t174 * t184;
t161 = t174 * t169 + t177 * t184;
t147 = t160 * pkin(3) - t161 * pkin(9) + t156;
t183 = t178 * t188;
t164 = -qJD(3) + t183;
t189 = pkin(7) * t183 + t175 * t186;
t157 = t169 * pkin(8) + t189;
t158 = (-pkin(2) * t178 - pkin(8) * t175 - pkin(1)) * t188;
t190 = t177 * t157 + t174 * t158;
t150 = -t164 * pkin(9) + t190;
t173 = sin(qJ(4));
t176 = cos(qJ(4));
t191 = t173 * t147 + t176 * t150;
t185 = t178 * t192;
t182 = -t174 * t157 + t177 * t158;
t181 = t176 * t147 - t173 * t150;
t149 = t164 * pkin(3) - t182;
t159 = qJD(4) + t160;
t152 = t176 * t161 - t173 * t164;
t151 = t173 * t161 + t176 * t164;
t145 = t151 * pkin(4) - t152 * qJ(5) + t149;
t144 = t159 * qJ(5) + t191;
t143 = -t159 * pkin(4) + qJD(5) - t181;
t1 = [t179 / 0.2e1, 0, 0, t175 ^ 2 * t192 / 0.2e1, t175 * t185, t169 * t184, t169 * t183, t169 ^ 2 / 0.2e1, pkin(1) * t185 + t180 * t169, -pkin(1) * t175 * t192 - t189 * t169, t161 ^ 2 / 0.2e1, -t161 * t160, -t161 * t164, t160 * t164, t164 ^ 2 / 0.2e1, t156 * t160 - t182 * t164, t156 * t161 + t190 * t164, t152 ^ 2 / 0.2e1, -t152 * t151, t152 * t159, -t151 * t159, t159 ^ 2 / 0.2e1, t149 * t151 + t181 * t159, t149 * t152 - t191 * t159, -t143 * t159 + t145 * t151, t143 * t152 - t144 * t151, t144 * t159 - t145 * t152, t144 ^ 2 / 0.2e1 + t145 ^ 2 / 0.2e1 + t143 ^ 2 / 0.2e1;];
T_reg = t1;
