% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPR9_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR9_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:24:37
% EndTime: 2019-12-31 21:24:37
% DurationCPUTime: 0.10s
% Computational Cost: add. (258->40), mult. (595->87), div. (0->0), fcn. (412->8), ass. (0->38)
t162 = qJD(1) ^ 2;
t172 = t162 / 0.2e1;
t171 = cos(qJ(3));
t161 = cos(qJ(2));
t170 = t161 * t162;
t158 = sin(qJ(3));
t159 = sin(qJ(2));
t168 = qJD(1) * t159;
t145 = t158 * qJD(2) + t171 * t168;
t167 = t161 * qJD(1);
t151 = -qJD(3) + t167;
t143 = (-pkin(2) * t161 - pkin(7) * t159 - pkin(1)) * qJD(1);
t148 = pkin(6) * t167 + qJD(2) * pkin(7);
t163 = t171 * t143 - t158 * t148;
t134 = -t151 * pkin(3) - t145 * qJ(4) + t163;
t144 = -t171 * qJD(2) + t158 * t168;
t169 = t158 * t143 + t171 * t148;
t136 = -t144 * qJ(4) + t169;
t155 = sin(pkin(9));
t156 = cos(pkin(9));
t128 = t155 * t134 + t156 * t136;
t166 = qJD(1) * qJD(2);
t165 = t159 * t166;
t164 = t161 * t166;
t147 = -qJD(2) * pkin(2) + pkin(6) * t168;
t127 = t156 * t134 - t155 * t136;
t140 = t144 * pkin(3) + qJD(4) + t147;
t160 = cos(qJ(5));
t157 = sin(qJ(5));
t149 = -qJD(5) + t151;
t139 = -t155 * t144 + t156 * t145;
t138 = -t156 * t144 - t155 * t145;
t131 = -t138 * pkin(4) + t140;
t130 = t157 * t138 + t160 * t139;
t129 = -t160 * t138 + t157 * t139;
t126 = t138 * pkin(8) + t128;
t125 = -t151 * pkin(4) - t139 * pkin(8) + t127;
t1 = [t172, 0, 0, t159 ^ 2 * t172, t159 * t170, t165, t164, qJD(2) ^ 2 / 0.2e1, pkin(1) * t170 - pkin(6) * t165, -t162 * pkin(1) * t159 - pkin(6) * t164, t145 ^ 2 / 0.2e1, -t145 * t144, -t145 * t151, t144 * t151, t151 ^ 2 / 0.2e1, t147 * t144 - t163 * t151, t147 * t145 + t169 * t151, -t127 * t139 + t128 * t138, t128 ^ 2 / 0.2e1 + t127 ^ 2 / 0.2e1 + t140 ^ 2 / 0.2e1, t130 ^ 2 / 0.2e1, -t130 * t129, -t130 * t149, t129 * t149, t149 ^ 2 / 0.2e1, -(t160 * t125 - t157 * t126) * t149 + t131 * t129, (t157 * t125 + t160 * t126) * t149 + t131 * t130;];
T_reg = t1;
