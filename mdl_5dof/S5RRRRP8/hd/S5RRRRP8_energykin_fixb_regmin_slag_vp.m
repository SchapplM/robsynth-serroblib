% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRRP8
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
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRP8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP8_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:02:07
% EndTime: 2019-12-31 22:02:08
% DurationCPUTime: 0.12s
% Computational Cost: add. (194->36), mult. (448->79), div. (0->0), fcn. (292->6), ass. (0->34)
t141 = qJD(1) ^ 2;
t153 = t141 / 0.2e1;
t152 = cos(qJ(4));
t140 = cos(qJ(2));
t151 = t140 * t141;
t137 = sin(qJ(3));
t139 = cos(qJ(3));
t138 = sin(qJ(2));
t148 = qJD(1) * t138;
t126 = t137 * qJD(2) + t139 * t148;
t147 = t140 * qJD(1);
t132 = -qJD(3) + t147;
t124 = (-pkin(2) * t140 - pkin(7) * t138 - pkin(1)) * qJD(1);
t129 = pkin(6) * t147 + qJD(2) * pkin(7);
t142 = t139 * t124 - t137 * t129;
t115 = -t132 * pkin(3) - t126 * pkin(8) + t142;
t125 = -t139 * qJD(2) + t137 * t148;
t149 = t137 * t124 + t139 * t129;
t117 = -t125 * pkin(8) + t149;
t136 = sin(qJ(4));
t150 = t136 * t115 + t152 * t117;
t146 = qJD(1) * qJD(2);
t145 = t138 * t146;
t144 = t140 * t146;
t128 = -qJD(2) * pkin(2) + pkin(6) * t148;
t143 = t152 * t115 - t136 * t117;
t120 = t125 * pkin(3) + t128;
t130 = -qJD(4) + t132;
t119 = -t136 * t125 + t152 * t126;
t118 = t152 * t125 + t136 * t126;
t112 = t118 * pkin(4) + qJD(5) + t120;
t111 = -t118 * qJ(5) + t150;
t110 = -t130 * pkin(4) - t119 * qJ(5) + t143;
t1 = [t153, 0, 0, t138 ^ 2 * t153, t138 * t151, t145, t144, qJD(2) ^ 2 / 0.2e1, pkin(1) * t151 - pkin(6) * t145, -t141 * pkin(1) * t138 - pkin(6) * t144, t126 ^ 2 / 0.2e1, -t126 * t125, -t126 * t132, t125 * t132, t132 ^ 2 / 0.2e1, t128 * t125 - t142 * t132, t128 * t126 + t149 * t132, t119 ^ 2 / 0.2e1, -t119 * t118, -t119 * t130, t118 * t130, t130 ^ 2 / 0.2e1, t120 * t118 - t143 * t130, t120 * t119 + t150 * t130, -t110 * t119 - t111 * t118, t111 ^ 2 / 0.2e1 + t110 ^ 2 / 0.2e1 + t112 ^ 2 / 0.2e1;];
T_reg = t1;
