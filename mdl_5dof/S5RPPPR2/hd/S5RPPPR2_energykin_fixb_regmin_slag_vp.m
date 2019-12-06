% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPPR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:31:30
% EndTime: 2019-12-05 17:31:31
% DurationCPUTime: 0.10s
% Computational Cost: add. (208->43), mult. (553->95), div. (0->0), fcn. (353->8), ass. (0->35)
t161 = sin(pkin(7));
t163 = cos(pkin(8));
t173 = t161 * t163;
t164 = cos(pkin(7));
t149 = qJD(2) + (-pkin(2) * t164 - qJ(3) * t161 - pkin(1)) * qJD(1);
t160 = sin(pkin(8));
t171 = qJD(1) * t164;
t170 = qJ(2) * t171;
t142 = t160 * t149 + t163 * t170;
t138 = -qJ(4) * t171 + t142;
t172 = qJD(1) * t161;
t153 = qJ(2) * t172 + qJD(3);
t144 = (pkin(3) * t160 - qJ(4) * t163) * t172 + t153;
t159 = sin(pkin(9));
t162 = cos(pkin(9));
t135 = t162 * t138 + t159 * t144;
t169 = t160 * t172;
t141 = t163 * t149 - t160 * t170;
t134 = -t159 * t138 + t162 * t144;
t137 = pkin(3) * t171 + qJD(4) - t141;
t147 = (t159 * t173 + t162 * t164) * qJD(1);
t167 = qJD(1) ^ 2;
t166 = cos(qJ(5));
t165 = sin(qJ(5));
t158 = t164 ^ 2;
t157 = t161 ^ 2;
t156 = -qJD(1) * pkin(1) + qJD(2);
t148 = (-t159 * t164 + t162 * t173) * qJD(1);
t145 = qJD(5) + t147;
t140 = t166 * t148 + t165 * t169;
t139 = t165 * t148 - t166 * t169;
t133 = pkin(6) * t169 + t135;
t132 = -pkin(4) * t169 - t134;
t131 = t147 * pkin(4) - t148 * pkin(6) + t137;
t1 = [t167 / 0.2e1, 0, 0, -t156 * t171, t156 * t172, (t157 + t158) * t167 * qJ(2), t156 ^ 2 / 0.2e1 + (t158 / 0.2e1 + t157 / 0.2e1) * qJ(2) ^ 2 * t167, (t153 * t160 * t161 - t141 * t164) * qJD(1), (t142 * t164 + t153 * t173) * qJD(1), (-t141 * t163 - t142 * t160) * t172, t142 ^ 2 / 0.2e1 + t141 ^ 2 / 0.2e1 + t153 ^ 2 / 0.2e1, t134 * t169 + t137 * t147, -t135 * t169 + t137 * t148, -t134 * t148 - t135 * t147, t135 ^ 2 / 0.2e1 + t134 ^ 2 / 0.2e1 + t137 ^ 2 / 0.2e1, t140 ^ 2 / 0.2e1, -t140 * t139, t140 * t145, -t139 * t145, t145 ^ 2 / 0.2e1, (t166 * t131 - t165 * t133) * t145 + t132 * t139, -(t165 * t131 + t166 * t133) * t145 + t132 * t140;];
T_reg = t1;
