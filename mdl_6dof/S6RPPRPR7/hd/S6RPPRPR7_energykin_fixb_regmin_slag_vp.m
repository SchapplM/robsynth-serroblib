% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRPR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR7_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:53:43
% EndTime: 2019-03-09 01:53:43
% DurationCPUTime: 0.12s
% Computational Cost: add. (325->45), mult. (687->97), div. (0->0), fcn. (446->8), ass. (0->38)
t164 = qJD(1) ^ 2;
t170 = t164 / 0.2e1;
t169 = cos(pkin(10));
t158 = sin(pkin(9));
t149 = qJD(2) + (-pkin(1) - qJ(3)) * qJD(1);
t166 = -pkin(7) * qJD(1) + t149;
t142 = t166 * t158;
t159 = cos(pkin(9));
t143 = t166 * t159;
t161 = sin(qJ(4));
t163 = cos(qJ(4));
t168 = t163 * t142 + t161 * t143;
t134 = qJD(4) * qJ(5) + t168;
t145 = (t158 * t163 + t159 * t161) * qJD(1);
t146 = (-t158 * t161 + t159 * t163) * qJD(1);
t151 = qJD(1) * qJ(2) + qJD(3);
t167 = qJD(1) * t158;
t147 = pkin(3) * t167 + t151;
t135 = t145 * pkin(4) - t146 * qJ(5) + t147;
t157 = sin(pkin(10));
t126 = t169 * t134 + t157 * t135;
t125 = -t157 * t134 + t169 * t135;
t165 = -t161 * t142 + t163 * t143;
t133 = -qJD(4) * pkin(4) + qJD(5) - t165;
t162 = cos(qJ(6));
t160 = sin(qJ(6));
t155 = t159 ^ 2;
t154 = t158 ^ 2;
t152 = -qJD(1) * pkin(1) + qJD(2);
t144 = qJD(6) + t145;
t138 = t157 * qJD(4) + t169 * t146;
t137 = -t169 * qJD(4) + t157 * t146;
t129 = -t160 * t137 + t162 * t138;
t128 = t162 * t137 + t160 * t138;
t127 = t137 * pkin(5) + t133;
t124 = -t137 * pkin(8) + t126;
t123 = t145 * pkin(5) - t138 * pkin(8) + t125;
t1 = [t170, 0, 0, t152 * qJD(1), t164 * qJ(2), qJ(2) ^ 2 * t170 + t152 ^ 2 / 0.2e1, t151 * t167, t151 * t159 * qJD(1) (-t154 - t155) * t149 * qJD(1), t151 ^ 2 / 0.2e1 + (t154 / 0.2e1 + t155 / 0.2e1) * t149 ^ 2, t146 ^ 2 / 0.2e1, -t146 * t145, t146 * qJD(4), -t145 * qJD(4), qJD(4) ^ 2 / 0.2e1, t165 * qJD(4) + t147 * t145, -t168 * qJD(4) + t147 * t146, t125 * t145 + t133 * t137, -t126 * t145 + t133 * t138, -t125 * t138 - t126 * t137, t126 ^ 2 / 0.2e1 + t125 ^ 2 / 0.2e1 + t133 ^ 2 / 0.2e1, t129 ^ 2 / 0.2e1, -t129 * t128, t129 * t144, -t128 * t144, t144 ^ 2 / 0.2e1 (t162 * t123 - t160 * t124) * t144 + t127 * t128 -(t160 * t123 + t162 * t124) * t144 + t127 * t129;];
T_reg  = t1;
