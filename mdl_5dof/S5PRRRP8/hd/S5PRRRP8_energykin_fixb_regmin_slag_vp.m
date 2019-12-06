% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRRP8
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
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRP8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP8_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:00:35
% EndTime: 2019-12-05 17:00:35
% DurationCPUTime: 0.12s
% Computational Cost: add. (132->31), mult. (309->72), div. (0->0), fcn. (202->8), ass. (0->31)
t145 = qJD(2) ^ 2;
t157 = t145 / 0.2e1;
t141 = sin(qJ(2));
t154 = qJD(1) * sin(pkin(5));
t131 = qJD(2) * pkin(7) + t141 * t154;
t140 = sin(qJ(3));
t143 = cos(qJ(3));
t153 = qJD(1) * cos(pkin(5));
t155 = t143 * t131 + t140 * t153;
t124 = qJD(3) * pkin(8) + t155;
t144 = cos(qJ(2));
t149 = t144 * t154;
t126 = -t149 + (-pkin(3) * t143 - pkin(8) * t140 - pkin(2)) * qJD(2);
t139 = sin(qJ(4));
t142 = cos(qJ(4));
t156 = t142 * t124 + t139 * t126;
t152 = qJD(2) * t140;
t151 = t143 * qJD(2);
t150 = qJD(2) * qJD(3);
t148 = qJD(2) * t154;
t147 = -t140 * t131 + t143 * t153;
t146 = -t139 * t124 + t142 * t126;
t123 = -qJD(3) * pkin(3) - t147;
t134 = -qJD(4) + t151;
t132 = -qJD(2) * pkin(2) - t149;
t130 = t139 * qJD(3) + t142 * t152;
t129 = -t142 * qJD(3) + t139 * t152;
t121 = t129 * pkin(4) - t130 * qJ(5) + t123;
t120 = -t134 * qJ(5) + t156;
t119 = t134 * pkin(4) + qJD(5) - t146;
t1 = [qJD(1) ^ 2 / 0.2e1, t157, t144 * t148, -t141 * t148, t140 ^ 2 * t157, t140 * t145 * t143, t140 * t150, t143 * t150, qJD(3) ^ 2 / 0.2e1, t147 * qJD(3) - t132 * t151, -t155 * qJD(3) + t132 * t152, t130 ^ 2 / 0.2e1, -t130 * t129, -t130 * t134, t129 * t134, t134 ^ 2 / 0.2e1, t123 * t129 - t146 * t134, t123 * t130 + t156 * t134, t119 * t134 + t121 * t129, t119 * t130 - t120 * t129, -t120 * t134 - t121 * t130, t120 ^ 2 / 0.2e1 + t121 ^ 2 / 0.2e1 + t119 ^ 2 / 0.2e1;];
T_reg = t1;
