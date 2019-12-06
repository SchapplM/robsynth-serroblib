% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:38:50
% EndTime: 2019-12-05 18:38:50
% DurationCPUTime: 0.13s
% Computational Cost: add. (258->39), mult. (663->86), div. (0->0), fcn. (480->8), ass. (0->39)
t155 = -pkin(7) - pkin(6);
t144 = qJD(1) ^ 2;
t154 = t144 / 0.2e1;
t153 = cos(qJ(3));
t143 = cos(qJ(2));
t152 = t143 * t144;
t140 = sin(qJ(3));
t141 = sin(qJ(2));
t129 = (t140 * t143 + t153 * t141) * qJD(1);
t136 = qJD(2) + qJD(3);
t150 = qJD(1) * t141;
t131 = qJD(2) * pkin(2) + t155 * t150;
t149 = qJD(1) * t143;
t132 = t155 * t149;
t145 = t153 * t131 + t140 * t132;
t119 = t136 * pkin(3) - t129 * qJ(4) + t145;
t128 = t140 * t150 - t153 * t149;
t151 = t140 * t131 - t153 * t132;
t121 = -t128 * qJ(4) + t151;
t137 = sin(pkin(9));
t138 = cos(pkin(9));
t113 = t137 * t119 + t138 * t121;
t148 = qJD(1) * qJD(2);
t147 = t141 * t148;
t146 = t143 * t148;
t112 = t138 * t119 - t137 * t121;
t133 = (-pkin(2) * t143 - pkin(1)) * qJD(1);
t125 = t128 * pkin(3) + qJD(4) + t133;
t142 = cos(qJ(5));
t139 = sin(qJ(5));
t135 = qJD(5) + t136;
t124 = -t137 * t128 + t138 * t129;
t123 = -t138 * t128 - t137 * t129;
t116 = -t123 * pkin(4) + t125;
t115 = t139 * t123 + t142 * t124;
t114 = -t142 * t123 + t139 * t124;
t111 = t123 * pkin(8) + t113;
t110 = t136 * pkin(4) - t124 * pkin(8) + t112;
t1 = [t154, 0, 0, t141 ^ 2 * t154, t141 * t152, t147, t146, qJD(2) ^ 2 / 0.2e1, pkin(1) * t152 - pkin(6) * t147, -t144 * pkin(1) * t141 - pkin(6) * t146, t129 ^ 2 / 0.2e1, -t129 * t128, t129 * t136, -t128 * t136, t136 ^ 2 / 0.2e1, t133 * t128 + t136 * t145, t133 * t129 - t151 * t136, -t112 * t124 + t113 * t123, t113 ^ 2 / 0.2e1 + t112 ^ 2 / 0.2e1 + t125 ^ 2 / 0.2e1, t115 ^ 2 / 0.2e1, -t115 * t114, t115 * t135, -t114 * t135, t135 ^ 2 / 0.2e1, t116 * t114 + (t142 * t110 - t139 * t111) * t135, t116 * t115 - (t139 * t110 + t142 * t111) * t135;];
T_reg = t1;
