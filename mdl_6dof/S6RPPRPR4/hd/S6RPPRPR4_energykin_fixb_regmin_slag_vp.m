% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPPRPR4
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
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRPR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:47:08
% EndTime: 2019-03-09 01:47:08
% DurationCPUTime: 0.10s
% Computational Cost: add. (216->44), mult. (413->89), div. (0->0), fcn. (222->8), ass. (0->36)
t161 = qJD(1) ^ 2;
t167 = t161 / 0.2e1;
t143 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t153 = sin(pkin(9));
t155 = cos(pkin(9));
t165 = qJ(2) * qJD(1);
t137 = t153 * t143 + t155 * t165;
t135 = -qJD(1) * pkin(7) + t137;
t160 = cos(qJ(4));
t151 = t160 * qJD(3);
t158 = sin(qJ(4));
t128 = qJD(4) * pkin(4) + t151 + (qJ(5) * qJD(1) - t135) * t158;
t163 = qJD(1) * t160;
t166 = t158 * qJD(3) + t160 * t135;
t129 = -qJ(5) * t163 + t166;
t152 = sin(pkin(10));
t154 = cos(pkin(10));
t124 = t152 * t128 + t154 * t129;
t164 = qJD(1) * t158;
t162 = qJD(1) * qJD(4);
t136 = t155 * t143 - t153 * t165;
t134 = qJD(1) * pkin(3) - t136;
t140 = t152 * t164 - t154 * t163;
t123 = t154 * t128 - t152 * t129;
t130 = pkin(4) * t163 + qJD(5) + t134;
t159 = cos(qJ(6));
t157 = sin(qJ(6));
t148 = -qJD(1) * pkin(1) + qJD(2);
t141 = (t152 * t160 + t154 * t158) * qJD(1);
t138 = -qJD(6) + t140;
t133 = t157 * qJD(4) - t159 * t141;
t132 = -t159 * qJD(4) - t157 * t141;
t125 = -t140 * pkin(5) + t141 * pkin(8) + t130;
t122 = qJD(4) * pkin(8) + t124;
t121 = -qJD(4) * pkin(5) - t123;
t1 = [t167, 0, 0, -t148 * qJD(1), t161 * qJ(2), qJ(2) ^ 2 * t167 + t148 ^ 2 / 0.2e1, -t136 * qJD(1), t137 * qJD(1), t137 ^ 2 / 0.2e1 + t136 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t158 ^ 2 * t167, t158 * t161 * t160, -t158 * t162, -t160 * t162, qJD(4) ^ 2 / 0.2e1, t134 * t163 + (-t158 * t135 + t151) * qJD(4), -t166 * qJD(4) - t134 * t164, t123 * t141 + t124 * t140, t124 ^ 2 / 0.2e1 + t123 ^ 2 / 0.2e1 + t130 ^ 2 / 0.2e1, t133 ^ 2 / 0.2e1, -t133 * t132, -t133 * t138, t132 * t138, t138 ^ 2 / 0.2e1 -(-t157 * t122 + t159 * t125) * t138 + t121 * t132 (t159 * t122 + t157 * t125) * t138 + t121 * t133;];
T_reg  = t1;
