% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPR8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR8_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:20:25
% EndTime: 2019-12-31 21:20:25
% DurationCPUTime: 0.08s
% Computational Cost: add. (178->38), mult. (416->82), div. (0->0), fcn. (262->6), ass. (0->36)
t148 = pkin(3) + pkin(8);
t147 = -pkin(7) - pkin(6);
t135 = qJD(1) ^ 2;
t146 = t135 / 0.2e1;
t134 = cos(qJ(2));
t145 = t134 * t135;
t131 = sin(qJ(2));
t143 = qJD(1) * t131;
t124 = qJD(2) * pkin(2) + t147 * t143;
t142 = qJD(1) * t134;
t125 = t147 * t142;
t130 = sin(qJ(3));
t133 = cos(qJ(3));
t144 = t130 * t124 - t133 * t125;
t141 = qJD(1) * qJD(2);
t140 = t131 * t141;
t139 = t134 * t141;
t138 = t133 * t124 + t130 * t125;
t128 = qJD(2) + qJD(3);
t114 = -t128 * qJ(4) - t144;
t126 = (-pkin(2) * t134 - pkin(1)) * qJD(1);
t137 = qJD(4) - t138;
t121 = (t130 * t134 + t131 * t133) * qJD(1);
t136 = -t121 * qJ(4) + t126;
t132 = cos(qJ(5));
t129 = sin(qJ(5));
t120 = t130 * t143 - t133 * t142;
t119 = qJD(5) + t121;
t116 = t129 * t120 + t132 * t128;
t115 = -t132 * t120 + t129 * t128;
t113 = -t128 * pkin(3) + t137;
t112 = t120 * pkin(3) + t136;
t111 = -t120 * pkin(4) - t114;
t110 = t148 * t120 + t136;
t109 = t121 * pkin(4) - t148 * t128 + t137;
t1 = [t146, 0, 0, t131 ^ 2 * t146, t131 * t145, t140, t139, qJD(2) ^ 2 / 0.2e1, pkin(1) * t145 - pkin(6) * t140, -t135 * pkin(1) * t131 - pkin(6) * t139, t121 ^ 2 / 0.2e1, -t121 * t120, t121 * t128, -t120 * t128, t128 ^ 2 / 0.2e1, t126 * t120 + t138 * t128, t126 * t121 - t144 * t128, t113 * t121 + t114 * t120, -t112 * t120 + t113 * t128, -t112 * t121 - t114 * t128, t112 ^ 2 / 0.2e1 + t114 ^ 2 / 0.2e1 + t113 ^ 2 / 0.2e1, t116 ^ 2 / 0.2e1, -t116 * t115, t116 * t119, -t115 * t119, t119 ^ 2 / 0.2e1, (t132 * t109 - t129 * t110) * t119 + t111 * t115, -(t129 * t109 + t132 * t110) * t119 + t111 * t116;];
T_reg = t1;
