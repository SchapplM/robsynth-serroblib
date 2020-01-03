% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRP9_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP9_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:06:03
% EndTime: 2019-12-31 22:06:03
% DurationCPUTime: 0.12s
% Computational Cost: add. (260->38), mult. (565->83), div. (0->0), fcn. (368->6), ass. (0->34)
t146 = qJD(1) ^ 2;
t158 = t146 / 0.2e1;
t157 = cos(qJ(3));
t145 = cos(qJ(2));
t156 = t145 * t146;
t142 = sin(qJ(3));
t143 = sin(qJ(2));
t153 = qJD(1) * t143;
t131 = t142 * qJD(2) + t157 * t153;
t152 = t145 * qJD(1);
t137 = -qJD(3) + t152;
t129 = (-pkin(2) * t145 - pkin(7) * t143 - pkin(1)) * qJD(1);
t134 = pkin(6) * t152 + qJD(2) * pkin(7);
t148 = t157 * t129 - t142 * t134;
t120 = -t137 * pkin(3) - t131 * pkin(8) + t148;
t130 = -t157 * qJD(2) + t142 * t153;
t154 = t142 * t129 + t157 * t134;
t122 = -t130 * pkin(8) + t154;
t141 = sin(qJ(4));
t144 = cos(qJ(4));
t155 = t141 * t120 + t144 * t122;
t151 = qJD(1) * qJD(2);
t150 = t143 * t151;
t149 = t145 * t151;
t133 = -qJD(2) * pkin(2) + pkin(6) * t153;
t147 = t144 * t120 - t141 * t122;
t125 = t130 * pkin(3) + t133;
t135 = -qJD(4) + t137;
t124 = -t141 * t130 + t144 * t131;
t123 = t144 * t130 + t141 * t131;
t118 = t123 * pkin(4) - t124 * qJ(5) + t125;
t117 = -t135 * qJ(5) + t155;
t116 = t135 * pkin(4) + qJD(5) - t147;
t1 = [t158, 0, 0, t143 ^ 2 * t158, t143 * t156, t150, t149, qJD(2) ^ 2 / 0.2e1, pkin(1) * t156 - pkin(6) * t150, -t146 * pkin(1) * t143 - pkin(6) * t149, t131 ^ 2 / 0.2e1, -t131 * t130, -t131 * t137, t130 * t137, t137 ^ 2 / 0.2e1, t133 * t130 - t148 * t137, t133 * t131 + t154 * t137, t124 ^ 2 / 0.2e1, -t124 * t123, -t124 * t135, t123 * t135, t135 ^ 2 / 0.2e1, t125 * t123 - t147 * t135, t125 * t124 + t155 * t135, t116 * t135 + t118 * t123, t116 * t124 - t117 * t123, -t117 * t135 - t118 * t124, t117 ^ 2 / 0.2e1 + t118 ^ 2 / 0.2e1 + t116 ^ 2 / 0.2e1;];
T_reg = t1;
