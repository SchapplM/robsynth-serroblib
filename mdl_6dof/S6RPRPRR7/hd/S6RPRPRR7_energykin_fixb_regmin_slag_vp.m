% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR7_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:56:03
% EndTime: 2019-03-09 03:56:03
% DurationCPUTime: 0.10s
% Computational Cost: add. (309->45), mult. (655->94), div. (0->0), fcn. (442->8), ass. (0->38)
t164 = qJD(1) ^ 2;
t171 = t164 / 0.2e1;
t170 = t164 * qJ(2);
t163 = cos(qJ(3));
t150 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t166 = -qJ(4) * qJD(1) + t150;
t144 = qJD(3) * pkin(3) + t166 * t163;
t160 = sin(qJ(3));
t146 = t166 * t160;
t156 = sin(pkin(10));
t157 = cos(pkin(10));
t135 = t157 * t144 - t156 * t146;
t148 = (-t156 * t160 + t157 * t163) * qJD(1);
t131 = qJD(3) * pkin(4) - t148 * pkin(8) + t135;
t136 = t156 * t144 + t157 * t146;
t147 = (-t156 * t163 - t157 * t160) * qJD(1);
t132 = t147 * pkin(8) + t136;
t159 = sin(qJ(5));
t162 = cos(qJ(5));
t169 = t159 * t131 + t162 * t132;
t168 = qJD(3) * t150;
t167 = qJD(1) * qJD(3);
t149 = qJD(4) + (pkin(3) * t160 + qJ(2)) * qJD(1);
t138 = -t162 * t147 + t159 * t148;
t165 = t162 * t131 - t159 * t132;
t140 = -t147 * pkin(4) + t149;
t161 = cos(qJ(6));
t158 = sin(qJ(6));
t154 = qJD(3) + qJD(5);
t153 = -qJD(1) * pkin(1) + qJD(2);
t139 = t159 * t147 + t162 * t148;
t137 = qJD(6) + t138;
t134 = t161 * t139 + t158 * t154;
t133 = t158 * t139 - t161 * t154;
t128 = t138 * pkin(5) - t139 * pkin(9) + t140;
t127 = t154 * pkin(9) + t169;
t126 = -t154 * pkin(5) - t165;
t1 = [t171, 0, 0, t153 * qJD(1), t170, qJ(2) ^ 2 * t171 + t153 ^ 2 / 0.2e1, t163 ^ 2 * t171, -t163 * t164 * t160, t163 * t167, -t160 * t167, qJD(3) ^ 2 / 0.2e1, t160 * t170 + t163 * t168, -t160 * t168 + t163 * t170, -t135 * t148 + t136 * t147, t136 ^ 2 / 0.2e1 + t135 ^ 2 / 0.2e1 + t149 ^ 2 / 0.2e1, t139 ^ 2 / 0.2e1, -t139 * t138, t139 * t154, -t138 * t154, t154 ^ 2 / 0.2e1, t140 * t138 + t165 * t154, t140 * t139 - t169 * t154, t134 ^ 2 / 0.2e1, -t134 * t133, t134 * t137, -t133 * t137, t137 ^ 2 / 0.2e1 (-t158 * t127 + t161 * t128) * t137 + t126 * t133 -(t161 * t127 + t158 * t128) * t137 + t126 * t134;];
T_reg  = t1;
