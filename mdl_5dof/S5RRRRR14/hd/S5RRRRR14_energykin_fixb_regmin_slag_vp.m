% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x27]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRR14_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR14_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:43:31
% EndTime: 2024-09-27 18:43:31
% DurationCPUTime: 0.03s
% Computational Cost: add. (303->36), mult. (500->86), div. (0->0), fcn. (356->10), ass. (0->41)
t175 = cos(qJ(4));
t174 = pkin(1) * qJD(1);
t152 = qJD(1) + qJD(2);
t166 = cos(qJ(2)) * t174;
t147 = t152 * pkin(2) + t166;
t154 = cos(pkin(5));
t173 = t147 * t154;
t150 = t152 ^ 2;
t153 = sin(pkin(5));
t151 = t153 ^ 2;
t172 = t150 * t151;
t171 = t152 * t153;
t160 = cos(qJ(3));
t170 = t152 * t160;
t142 = t160 * t173;
t167 = sin(qJ(2)) * t174;
t144 = pkin(8) * t171 + t167;
t149 = t154 * t152 + qJD(3);
t157 = sin(qJ(3));
t133 = t149 * pkin(3) + t142 + (-pkin(9) * t171 - t144) * t157;
t163 = t153 * t170;
t168 = t160 * t144 + t157 * t173;
t136 = pkin(9) * t163 + t168;
t156 = sin(qJ(4));
t169 = t156 * t133 + t175 * t136;
t165 = t147 * t151 * t152;
t164 = t157 * t171;
t162 = t175 * t133 - t156 * t136;
t148 = qJD(4) + t149;
t138 = (-pkin(3) * t170 - t147) * t153;
t159 = cos(qJ(5));
t155 = sin(qJ(5));
t146 = qJD(5) + t148;
t140 = (t156 * t160 + t175 * t157) * t171;
t139 = t156 * t164 - t175 * t163;
t134 = t139 * pkin(4) + t138;
t132 = -t155 * t139 + t159 * t140;
t131 = t159 * t139 + t155 * t140;
t128 = -t139 * pkin(10) + t169;
t127 = t148 * pkin(4) - t140 * pkin(10) + t162;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t150 / 0.2e1, t152 * t166, -t152 * t167, t157 ^ 2 * t172 / 0.2e1, t157 * t160 * t172, t149 * t164, t149 * t163, t149 ^ 2 / 0.2e1, t160 * t165 + (-t157 * t144 + t142) * t149, -t168 * t149 - t157 * t165, t140 ^ 2 / 0.2e1, -t140 * t139, t140 * t148, -t139 * t148, t148 ^ 2 / 0.2e1, t138 * t139 + t162 * t148, t138 * t140 - t169 * t148, t132 ^ 2 / 0.2e1, -t132 * t131, t132 * t146, -t131 * t146, t146 ^ 2 / 0.2e1, (t159 * t127 - t155 * t128) * t146 + t134 * t131, -(t155 * t127 + t159 * t128) * t146 + t134 * t132;];
T_reg = t1;
