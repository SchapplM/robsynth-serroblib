% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% T_reg [1x20]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRP7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:56:25
% EndTime: 2019-12-05 16:56:25
% DurationCPUTime: 0.07s
% Computational Cost: add. (94->29), mult. (239->68), div. (0->0), fcn. (155->8), ass. (0->31)
t138 = qJD(2) ^ 2;
t151 = t138 / 0.2e1;
t150 = cos(qJ(4));
t135 = sin(qJ(2));
t147 = qJD(1) * sin(pkin(5));
t125 = qJD(2) * pkin(7) + t135 * t147;
t134 = sin(qJ(3));
t136 = cos(qJ(3));
t146 = qJD(1) * cos(pkin(5));
t148 = t136 * t125 + t134 * t146;
t117 = qJD(3) * pkin(8) + t148;
t137 = cos(qJ(2));
t142 = t137 * t147;
t120 = -t142 + (-pkin(3) * t136 - pkin(8) * t134 - pkin(2)) * qJD(2);
t133 = sin(qJ(4));
t149 = t150 * t117 + t133 * t120;
t145 = qJD(2) * t134;
t144 = t136 * qJD(2);
t143 = qJD(2) * qJD(3);
t141 = qJD(2) * t147;
t140 = -t133 * t117 + t150 * t120;
t139 = -t134 * t125 + t136 * t146;
t116 = -qJD(3) * pkin(3) - t139;
t128 = -qJD(4) + t144;
t126 = -qJD(2) * pkin(2) - t142;
t124 = t133 * qJD(3) + t150 * t145;
t123 = -t150 * qJD(3) + t133 * t145;
t114 = t123 * pkin(4) + qJD(5) + t116;
t113 = -t123 * qJ(5) + t149;
t112 = -t128 * pkin(4) - t124 * qJ(5) + t140;
t1 = [qJD(1) ^ 2 / 0.2e1, t151, t137 * t141, -t135 * t141, t134 ^ 2 * t151, t134 * t138 * t136, t134 * t143, t136 * t143, qJD(3) ^ 2 / 0.2e1, t139 * qJD(3) - t126 * t144, -t148 * qJD(3) + t126 * t145, t124 ^ 2 / 0.2e1, -t124 * t123, -t124 * t128, t123 * t128, t128 ^ 2 / 0.2e1, t116 * t123 - t140 * t128, t116 * t124 + t149 * t128, -t112 * t124 - t113 * t123, t113 ^ 2 / 0.2e1 + t112 ^ 2 / 0.2e1 + t114 ^ 2 / 0.2e1;];
T_reg = t1;
