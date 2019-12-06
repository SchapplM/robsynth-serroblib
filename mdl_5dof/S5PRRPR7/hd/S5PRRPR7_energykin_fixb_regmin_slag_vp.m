% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRPR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR7_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:37:33
% EndTime: 2019-12-05 16:37:33
% DurationCPUTime: 0.09s
% Computational Cost: add. (148->35), mult. (364->81), div. (0->0), fcn. (251->10), ass. (0->35)
t160 = qJD(2) ^ 2;
t171 = t160 / 0.2e1;
t156 = sin(qJ(2));
t169 = qJD(1) * sin(pkin(5));
t146 = qJD(2) * pkin(7) + t156 * t169;
t155 = sin(qJ(3));
t158 = cos(qJ(3));
t168 = qJD(1) * cos(pkin(5));
t170 = t158 * t146 + t155 * t168;
t137 = qJD(3) * qJ(4) + t170;
t159 = cos(qJ(2));
t163 = t159 * t169;
t140 = -t163 + (-pkin(3) * t158 - qJ(4) * t155 - pkin(2)) * qJD(2);
t150 = sin(pkin(10));
t152 = cos(pkin(10));
t133 = t152 * t137 + t150 * t140;
t167 = qJD(2) * t155;
t166 = qJD(2) * t158;
t165 = (-qJD(2) * pkin(2) - t163) * qJD(2);
t164 = qJD(2) * qJD(3);
t162 = qJD(2) * t169;
t161 = -t155 * t146 + t158 * t168;
t144 = -t152 * qJD(3) + t150 * t167;
t132 = -t150 * t137 + t152 * t140;
t135 = -qJD(3) * pkin(3) + qJD(4) - t161;
t157 = cos(qJ(5));
t154 = sin(qJ(5));
t145 = t150 * qJD(3) + t152 * t167;
t141 = qJD(5) + t144;
t139 = t157 * t145 - t154 * t166;
t138 = t154 * t145 + t157 * t166;
t131 = t144 * pkin(4) - t145 * pkin(8) + t135;
t130 = -pkin(8) * t166 + t133;
t129 = pkin(4) * t166 - t132;
t1 = [qJD(1) ^ 2 / 0.2e1, t171, t159 * t162, -t156 * t162, t155 ^ 2 * t171, t155 * t160 * t158, t155 * t164, t158 * t164, qJD(3) ^ 2 / 0.2e1, t161 * qJD(3) - t158 * t165, -t170 * qJD(3) + t155 * t165, -t132 * t166 + t135 * t144, t133 * t166 + t135 * t145, -t132 * t145 - t133 * t144, t133 ^ 2 / 0.2e1 + t132 ^ 2 / 0.2e1 + t135 ^ 2 / 0.2e1, t139 ^ 2 / 0.2e1, -t139 * t138, t139 * t141, -t138 * t141, t141 ^ 2 / 0.2e1, (-t154 * t130 + t157 * t131) * t141 + t129 * t138, -(t157 * t130 + t154 * t131) * t141 + t129 * t139;];
T_reg = t1;
