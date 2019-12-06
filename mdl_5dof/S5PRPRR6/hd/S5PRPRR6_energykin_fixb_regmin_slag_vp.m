% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPRR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR6_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:57:52
% EndTime: 2019-12-05 15:57:53
% DurationCPUTime: 0.08s
% Computational Cost: add. (122->36), mult. (325->78), div. (0->0), fcn. (240->10), ass. (0->33)
t147 = sin(qJ(2));
t157 = qJD(1) * sin(pkin(5));
t136 = qJD(2) * qJ(3) + t147 * t157;
t143 = cos(pkin(10));
t156 = qJD(1) * cos(pkin(5));
t138 = t143 * t156;
t141 = sin(pkin(10));
t124 = t138 + (-pkin(7) * qJD(2) - t136) * t141;
t129 = t143 * t136 + t141 * t156;
t154 = qJD(2) * t143;
t125 = pkin(7) * t154 + t129;
t146 = sin(qJ(4));
t149 = cos(qJ(4));
t158 = t146 * t124 + t149 * t125;
t155 = qJD(2) * t141;
t153 = qJD(2) * t157;
t132 = t146 * t155 - t149 * t154;
t150 = cos(qJ(2));
t152 = -t150 * t157 + qJD(3);
t151 = t149 * t124 - t146 * t125;
t130 = (-pkin(3) * t143 - pkin(2)) * qJD(2) + t152;
t148 = cos(qJ(5));
t145 = sin(qJ(5));
t135 = -qJD(2) * pkin(2) + t152;
t133 = (t141 * t149 + t143 * t146) * qJD(2);
t131 = qJD(5) + t132;
t128 = -t141 * t136 + t138;
t127 = t145 * qJD(4) + t148 * t133;
t126 = -t148 * qJD(4) + t145 * t133;
t121 = t132 * pkin(4) - t133 * pkin(8) + t130;
t120 = qJD(4) * pkin(8) + t158;
t119 = -qJD(4) * pkin(4) - t151;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, t150 * t153, -t147 * t153, -t135 * t154, t135 * t155, (-t128 * t141 + t129 * t143) * qJD(2), t129 ^ 2 / 0.2e1 + t128 ^ 2 / 0.2e1 + t135 ^ 2 / 0.2e1, t133 ^ 2 / 0.2e1, -t133 * t132, t133 * qJD(4), -t132 * qJD(4), qJD(4) ^ 2 / 0.2e1, t151 * qJD(4) + t130 * t132, -t158 * qJD(4) + t130 * t133, t127 ^ 2 / 0.2e1, -t127 * t126, t127 * t131, -t126 * t131, t131 ^ 2 / 0.2e1, (-t145 * t120 + t148 * t121) * t131 + t119 * t126, -(t148 * t120 + t145 * t121) * t131 + t119 * t127;];
T_reg = t1;
