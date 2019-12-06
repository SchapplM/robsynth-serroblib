% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% T_reg [1x19]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PPRRR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:19:42
% EndTime: 2019-12-05 15:19:43
% DurationCPUTime: 0.10s
% Computational Cost: add. (94->28), mult. (265->72), div. (0->0), fcn. (208->12), ass. (0->36)
t152 = cos(pkin(5)) * qJD(1) + qJD(2);
t158 = sin(pkin(6));
t161 = cos(pkin(6));
t160 = cos(pkin(11));
t159 = sin(pkin(5));
t180 = qJD(1) * t159;
t174 = t160 * t180;
t185 = t152 * t158 + t161 * t174;
t164 = sin(qJ(3));
t167 = cos(qJ(3));
t157 = sin(pkin(11));
t175 = t157 * t180;
t184 = -t164 * t175 + t185 * t167;
t168 = qJD(3) ^ 2;
t183 = t168 / 0.2e1;
t176 = t185 * t164 + t167 * t175;
t143 = qJD(3) * pkin(8) + t176;
t145 = t161 * t152 - t158 * t174;
t163 = sin(qJ(4));
t166 = cos(qJ(4));
t181 = t166 * t143 + t163 * t145;
t179 = qJD(3) * t163;
t178 = t166 * qJD(3);
t177 = qJD(3) * qJD(4);
t172 = -t163 * t143 + t166 * t145;
t169 = qJD(1) ^ 2;
t165 = cos(qJ(5));
t162 = sin(qJ(5));
t153 = -qJD(5) + t178;
t149 = t162 * qJD(4) + t165 * t179;
t148 = -t165 * qJD(4) + t162 * t179;
t142 = -qJD(3) * pkin(3) - t184;
t140 = (-pkin(4) * t166 - pkin(9) * t163 - pkin(3)) * qJD(3) - t184;
t139 = qJD(4) * pkin(9) + t181;
t138 = -qJD(4) * pkin(4) - t172;
t1 = [t169 / 0.2e1, t152 ^ 2 / 0.2e1 + (t157 ^ 2 / 0.2e1 + t160 ^ 2 / 0.2e1) * t169 * t159 ^ 2, t183, t184 * qJD(3), -t176 * qJD(3), t163 ^ 2 * t183, t163 * t168 * t166, t163 * t177, t166 * t177, qJD(4) ^ 2 / 0.2e1, t172 * qJD(4) - t142 * t178, -t181 * qJD(4) + t142 * t179, t149 ^ 2 / 0.2e1, -t149 * t148, -t149 * t153, t148 * t153, t153 ^ 2 / 0.2e1, -(-t162 * t139 + t165 * t140) * t153 + t138 * t148, (t165 * t139 + t162 * t140) * t153 + t138 * t149;];
T_reg = t1;
