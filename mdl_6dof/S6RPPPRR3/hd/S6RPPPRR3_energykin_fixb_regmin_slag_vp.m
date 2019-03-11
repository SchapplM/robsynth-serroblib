% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
% 
% Output:
% T_reg [1x27]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPPRR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:33:58
% EndTime: 2019-03-09 01:33:58
% DurationCPUTime: 0.08s
% Computational Cost: add. (209->42), mult. (409->89), div. (0->0), fcn. (228->8), ass. (0->36)
t152 = qJD(1) ^ 2;
t158 = t152 / 0.2e1;
t134 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t144 = sin(pkin(9));
t146 = cos(pkin(9));
t156 = qJ(2) * qJD(1);
t128 = t144 * t134 + t146 * t156;
t126 = -qJD(1) * qJ(4) + t128;
t145 = cos(pkin(10));
t141 = t145 * qJD(3);
t143 = sin(pkin(10));
t117 = t141 + (pkin(7) * qJD(1) - t126) * t143;
t120 = t143 * qJD(3) + t145 * t126;
t154 = qJD(1) * t145;
t118 = -pkin(7) * t154 + t120;
t149 = sin(qJ(5));
t151 = cos(qJ(5));
t157 = t149 * t117 + t151 * t118;
t155 = qJD(1) * t143;
t127 = t146 * t134 - t144 * t156;
t131 = t149 * t155 - t151 * t154;
t153 = t151 * t117 - t149 * t118;
t125 = qJD(1) * pkin(3) + qJD(4) - t127;
t121 = pkin(4) * t154 + t125;
t150 = cos(qJ(6));
t148 = sin(qJ(6));
t139 = -qJD(1) * pkin(1) + qJD(2);
t132 = (-t143 * t151 - t145 * t149) * qJD(1);
t129 = -qJD(6) + t131;
t124 = t148 * qJD(5) + t150 * t132;
t123 = -t150 * qJD(5) + t148 * t132;
t119 = -t143 * t126 + t141;
t114 = -t131 * pkin(5) - t132 * pkin(8) + t121;
t113 = qJD(5) * pkin(8) + t157;
t112 = -qJD(5) * pkin(5) - t153;
t1 = [t158, 0, 0, -t139 * qJD(1), t152 * qJ(2), qJ(2) ^ 2 * t158 + t139 ^ 2 / 0.2e1, -t127 * qJD(1), t128 * qJD(1), t128 ^ 2 / 0.2e1 + t127 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t125 * t154, -t125 * t155 (t119 * t143 - t120 * t145) * qJD(1), t120 ^ 2 / 0.2e1 + t119 ^ 2 / 0.2e1 + t125 ^ 2 / 0.2e1, t132 ^ 2 / 0.2e1, t132 * t131, t132 * qJD(5), t131 * qJD(5), qJD(5) ^ 2 / 0.2e1, t153 * qJD(5) - t121 * t131, -t157 * qJD(5) + t121 * t132, t124 ^ 2 / 0.2e1, -t124 * t123, -t124 * t129, t123 * t129, t129 ^ 2 / 0.2e1 -(-t148 * t113 + t150 * t114) * t129 + t112 * t123 (t150 * t113 + t148 * t114) * t129 + t112 * t124;];
T_reg  = t1;
