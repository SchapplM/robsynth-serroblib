% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRPRR4
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
% T_reg [1x19]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPRR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:51:34
% EndTime: 2019-12-05 15:51:34
% DurationCPUTime: 0.06s
% Computational Cost: add. (76->27), mult. (201->67), div. (0->0), fcn. (136->10), ass. (0->31)
t147 = qJD(2) ^ 2;
t156 = t147 / 0.2e1;
t146 = cos(qJ(2));
t154 = qJD(1) * sin(pkin(5));
t130 = qJD(2) * pkin(2) + t146 * t154;
t138 = sin(pkin(10));
t140 = cos(pkin(10));
t143 = sin(qJ(2));
t150 = t143 * t154;
t126 = t138 * t130 + t140 * t150;
t124 = qJD(2) * pkin(7) + t126;
t134 = cos(pkin(5)) * qJD(1) + qJD(3);
t142 = sin(qJ(4));
t145 = cos(qJ(4));
t155 = t145 * t124 + t142 * t134;
t153 = qJD(2) * t142;
t152 = t145 * qJD(2);
t151 = qJD(2) * qJD(4);
t149 = qJD(2) * t154;
t125 = t140 * t130 - t138 * t150;
t148 = -t142 * t124 + t145 * t134;
t144 = cos(qJ(5));
t141 = sin(qJ(5));
t135 = -qJD(5) + t152;
t129 = t141 * qJD(4) + t144 * t153;
t128 = -t144 * qJD(4) + t141 * t153;
t123 = -qJD(2) * pkin(3) - t125;
t121 = (-pkin(4) * t145 - pkin(8) * t142 - pkin(3)) * qJD(2) - t125;
t120 = qJD(4) * pkin(8) + t155;
t119 = -qJD(4) * pkin(4) - t148;
t1 = [qJD(1) ^ 2 / 0.2e1, t156, t146 * t149, -t143 * t149, t126 ^ 2 / 0.2e1 + t125 ^ 2 / 0.2e1 + t134 ^ 2 / 0.2e1, t142 ^ 2 * t156, t142 * t147 * t145, t142 * t151, t145 * t151, qJD(4) ^ 2 / 0.2e1, t148 * qJD(4) - t123 * t152, -t155 * qJD(4) + t123 * t153, t129 ^ 2 / 0.2e1, -t129 * t128, -t129 * t135, t128 * t135, t135 ^ 2 / 0.2e1, -(-t141 * t120 + t144 * t121) * t135 + t119 * t128, (t144 * t120 + t141 * t121) * t135 + t119 * t129;];
T_reg = t1;
