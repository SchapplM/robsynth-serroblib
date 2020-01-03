% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP11_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP11_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:54:31
% EndTime: 2019-12-31 18:54:31
% DurationCPUTime: 0.12s
% Computational Cost: add. (216->39), mult. (548->81), div. (0->0), fcn. (368->6), ass. (0->32)
t148 = pkin(6) + qJ(2);
t137 = sin(qJ(3));
t139 = cos(qJ(3));
t135 = cos(pkin(8));
t144 = qJD(1) * t135;
t134 = sin(pkin(8));
t145 = qJD(1) * t134;
t124 = t137 * t145 - t139 * t144;
t125 = (t134 * t139 + t135 * t137) * qJD(1);
t128 = qJD(2) + (-pkin(2) * t135 - pkin(1)) * qJD(1);
t114 = t124 * pkin(3) - t125 * pkin(7) + t128;
t126 = t148 * t145;
t127 = t148 * t144;
t146 = -t137 * t126 + t139 * t127;
t117 = qJD(3) * pkin(7) + t146;
t136 = sin(qJ(4));
t138 = cos(qJ(4));
t147 = t136 * t114 + t138 * t117;
t143 = -t139 * t126 - t137 * t127;
t142 = t138 * t114 - t136 * t117;
t116 = -qJD(3) * pkin(3) - t143;
t140 = qJD(1) ^ 2;
t133 = t135 ^ 2;
t132 = t134 ^ 2;
t130 = -qJD(1) * pkin(1) + qJD(2);
t120 = qJD(4) + t124;
t119 = t136 * qJD(3) + t138 * t125;
t118 = -t138 * qJD(3) + t136 * t125;
t112 = t118 * pkin(4) - t119 * qJ(5) + t116;
t111 = t120 * qJ(5) + t147;
t110 = -t120 * pkin(4) + qJD(5) - t142;
t1 = [t140 / 0.2e1, 0, 0, -t130 * t144, t130 * t145, (t132 + t133) * t140 * qJ(2), t130 ^ 2 / 0.2e1 + (t133 / 0.2e1 + t132 / 0.2e1) * qJ(2) ^ 2 * t140, t125 ^ 2 / 0.2e1, -t125 * t124, t125 * qJD(3), -t124 * qJD(3), qJD(3) ^ 2 / 0.2e1, t143 * qJD(3) + t128 * t124, -t146 * qJD(3) + t128 * t125, t119 ^ 2 / 0.2e1, -t119 * t118, t119 * t120, -t118 * t120, t120 ^ 2 / 0.2e1, t116 * t118 + t142 * t120, t116 * t119 - t147 * t120, -t110 * t120 + t112 * t118, t110 * t119 - t111 * t118, t111 * t120 - t112 * t119, t111 ^ 2 / 0.2e1 + t112 ^ 2 / 0.2e1 + t110 ^ 2 / 0.2e1;];
T_reg = t1;
