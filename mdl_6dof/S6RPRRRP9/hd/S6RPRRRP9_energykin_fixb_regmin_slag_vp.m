% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRP9_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP9_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:28:30
% EndTime: 2019-03-09 06:28:30
% DurationCPUTime: 0.10s
% Computational Cost: add. (241->42), mult. (480->87), div. (0->0), fcn. (292->6), ass. (0->34)
t147 = qJD(1) ^ 2;
t157 = t147 / 0.2e1;
t156 = cos(qJ(5));
t155 = t147 * qJ(2);
t143 = sin(qJ(4));
t145 = cos(qJ(4));
t146 = cos(qJ(3));
t152 = qJD(1) * t146;
t136 = t143 * qJD(3) + t145 * t152;
t144 = sin(qJ(3));
t139 = t144 * qJD(1) + qJD(4);
t132 = (pkin(3) * t144 - pkin(8) * t146 + qJ(2)) * qJD(1);
t138 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t133 = qJD(3) * pkin(8) + t144 * t138;
t148 = t145 * t132 - t143 * t133;
t122 = t139 * pkin(4) - t136 * pkin(9) + t148;
t135 = -t145 * qJD(3) + t143 * t152;
t153 = t143 * t132 + t145 * t133;
t124 = -t135 * pkin(9) + t153;
t142 = sin(qJ(5));
t154 = t142 * t122 + t156 * t124;
t151 = qJD(3) * t138;
t150 = qJD(1) * qJD(3);
t149 = t156 * t122 - t142 * t124;
t134 = -qJD(3) * pkin(3) - t146 * t138;
t127 = t135 * pkin(4) + t134;
t140 = -qJD(1) * pkin(1) + qJD(2);
t137 = qJD(5) + t139;
t126 = -t142 * t135 + t156 * t136;
t125 = t156 * t135 + t142 * t136;
t119 = t125 * pkin(5) + qJD(6) + t127;
t118 = -t125 * qJ(6) + t154;
t117 = t137 * pkin(5) - t126 * qJ(6) + t149;
t1 = [t157, 0, 0, t140 * qJD(1), t155, qJ(2) ^ 2 * t157 + t140 ^ 2 / 0.2e1, t146 ^ 2 * t157, -t146 * t147 * t144, t146 * t150, -t144 * t150, qJD(3) ^ 2 / 0.2e1, t144 * t155 + t146 * t151, -t144 * t151 + t146 * t155, t136 ^ 2 / 0.2e1, -t136 * t135, t136 * t139, -t135 * t139, t139 ^ 2 / 0.2e1, t134 * t135 + t148 * t139, t134 * t136 - t153 * t139, t126 ^ 2 / 0.2e1, -t126 * t125, t126 * t137, -t125 * t137, t137 ^ 2 / 0.2e1, t127 * t125 + t149 * t137, t127 * t126 - t154 * t137, -t117 * t126 - t118 * t125, t118 ^ 2 / 0.2e1 + t117 ^ 2 / 0.2e1 + t119 ^ 2 / 0.2e1;];
T_reg  = t1;
