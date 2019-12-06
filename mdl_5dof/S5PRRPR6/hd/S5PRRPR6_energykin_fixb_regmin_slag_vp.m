% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRPR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR6_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:32:46
% EndTime: 2019-12-05 16:32:46
% DurationCPUTime: 0.09s
% Computational Cost: add. (152->35), mult. (372->80), div. (0->0), fcn. (259->10), ass. (0->35)
t157 = qJD(2) ^ 2;
t168 = t157 / 0.2e1;
t167 = cos(pkin(10));
t153 = sin(qJ(2));
t165 = qJD(1) * sin(pkin(5));
t142 = qJD(2) * pkin(7) + t153 * t165;
t152 = sin(qJ(3));
t155 = cos(qJ(3));
t164 = qJD(1) * cos(pkin(5));
t166 = t155 * t142 + t152 * t164;
t135 = qJD(3) * qJ(4) + t166;
t156 = cos(qJ(2));
t160 = t156 * t165;
t136 = -t160 + (-pkin(3) * t155 - qJ(4) * t152 - pkin(2)) * qJD(2);
t148 = sin(pkin(10));
t127 = t167 * t135 + t148 * t136;
t163 = qJD(2) * t152;
t162 = t155 * qJD(2);
t161 = qJD(2) * qJD(3);
t159 = qJD(2) * t165;
t126 = -t148 * t135 + t167 * t136;
t158 = -t152 * t142 + t155 * t164;
t132 = -qJD(3) * pkin(3) + qJD(4) - t158;
t154 = cos(qJ(5));
t151 = sin(qJ(5));
t145 = -qJD(5) + t162;
t143 = -qJD(2) * pkin(2) - t160;
t141 = t148 * qJD(3) + t167 * t163;
t140 = -t167 * qJD(3) + t148 * t163;
t130 = -t151 * t140 + t154 * t141;
t129 = t154 * t140 + t151 * t141;
t128 = t140 * pkin(4) + t132;
t125 = -t140 * pkin(8) + t127;
t124 = -pkin(4) * t162 - t141 * pkin(8) + t126;
t1 = [qJD(1) ^ 2 / 0.2e1, t168, t156 * t159, -t153 * t159, t152 ^ 2 * t168, t152 * t157 * t155, t152 * t161, t155 * t161, qJD(3) ^ 2 / 0.2e1, t158 * qJD(3) - t143 * t162, -t166 * qJD(3) + t143 * t163, -t126 * t162 + t132 * t140, t127 * t162 + t132 * t141, -t126 * t141 - t127 * t140, t127 ^ 2 / 0.2e1 + t126 ^ 2 / 0.2e1 + t132 ^ 2 / 0.2e1, t130 ^ 2 / 0.2e1, -t130 * t129, -t130 * t145, t129 * t145, t145 ^ 2 / 0.2e1, -(t154 * t124 - t151 * t125) * t145 + t128 * t129, t128 * t130 + (t151 * t124 + t154 * t125) * t145;];
T_reg = t1;
