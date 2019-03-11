% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR7_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:21:00
% EndTime: 2019-03-09 05:21:00
% DurationCPUTime: 0.10s
% Computational Cost: add. (331->45), mult. (672->94), div. (0->0), fcn. (440->8), ass. (0->38)
t163 = qJD(1) ^ 2;
t170 = t163 / 0.2e1;
t169 = t163 * qJ(2);
t158 = sin(qJ(4));
t159 = sin(qJ(3));
t161 = cos(qJ(4));
t162 = cos(qJ(3));
t147 = (-t158 * t159 + t161 * t162) * qJD(1);
t153 = qJD(3) + qJD(4);
t149 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t165 = -pkin(8) * qJD(1) + t149;
t143 = qJD(3) * pkin(3) + t165 * t162;
t145 = t165 * t159;
t164 = t161 * t143 - t158 * t145;
t131 = t153 * pkin(4) - t147 * qJ(5) + t164;
t146 = (t158 * t162 + t159 * t161) * qJD(1);
t168 = t158 * t143 + t161 * t145;
t133 = -t146 * qJ(5) + t168;
t155 = sin(pkin(10));
t156 = cos(pkin(10));
t128 = t155 * t131 + t156 * t133;
t148 = (pkin(3) * t159 + qJ(2)) * qJD(1);
t167 = qJD(3) * t149;
t166 = qJD(1) * qJD(3);
t137 = -t156 * t146 - t155 * t147;
t127 = t156 * t131 - t155 * t133;
t139 = t146 * pkin(4) + qJD(5) + t148;
t160 = cos(qJ(6));
t157 = sin(qJ(6));
t152 = -qJD(1) * pkin(1) + qJD(2);
t138 = -t155 * t146 + t156 * t147;
t136 = qJD(6) - t137;
t135 = t160 * t138 + t157 * t153;
t134 = t157 * t138 - t160 * t153;
t129 = -t137 * pkin(5) - t138 * pkin(9) + t139;
t126 = t153 * pkin(9) + t128;
t125 = -t153 * pkin(5) - t127;
t1 = [t170, 0, 0, t152 * qJD(1), t169, qJ(2) ^ 2 * t170 + t152 ^ 2 / 0.2e1, t162 ^ 2 * t170, -t162 * t163 * t159, t162 * t166, -t159 * t166, qJD(3) ^ 2 / 0.2e1, t159 * t169 + t162 * t167, -t159 * t167 + t162 * t169, t147 ^ 2 / 0.2e1, -t147 * t146, t147 * t153, -t146 * t153, t153 ^ 2 / 0.2e1, t148 * t146 + t164 * t153, t148 * t147 - t168 * t153, -t127 * t138 + t128 * t137, t128 ^ 2 / 0.2e1 + t127 ^ 2 / 0.2e1 + t139 ^ 2 / 0.2e1, t135 ^ 2 / 0.2e1, -t135 * t134, t135 * t136, -t134 * t136, t136 ^ 2 / 0.2e1 (-t157 * t126 + t160 * t129) * t136 + t125 * t134 -(t160 * t126 + t157 * t129) * t136 + t125 * t135;];
T_reg  = t1;
