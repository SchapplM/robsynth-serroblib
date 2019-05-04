% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRRR10V2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
% 
% Output:
% T_reg [1x38]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRR10V2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10V2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:51:47
% EndTime: 2019-04-11 14:51:47
% DurationCPUTime: 0.14s
% Computational Cost: add. (372->42), mult. (775->99), div. (0->0), fcn. (619->10), ass. (0->40)
t168 = qJD(1) ^ 2;
t175 = t168 / 0.2e1;
t174 = pkin(2) * qJD(2);
t167 = cos(qJ(2));
t173 = t167 * t168;
t161 = sin(qJ(3));
t162 = sin(qJ(2));
t166 = cos(qJ(3));
t149 = (t161 * t162 - t166 * t167) * qJD(1);
t150 = (t161 * t167 + t162 * t166) * qJD(1);
t154 = (-pkin(2) * t167 - pkin(1)) * qJD(1);
t142 = t149 * pkin(3) - t150 * pkin(5) + t154;
t157 = qJD(2) + qJD(3);
t170 = t161 * t174;
t152 = t157 * pkin(5) + t170;
t160 = sin(qJ(4));
t165 = cos(qJ(4));
t138 = t160 * t142 + t165 * t152;
t169 = t166 * t174;
t153 = -t157 * pkin(3) - t169;
t159 = sin(qJ(5));
t164 = cos(qJ(5));
t172 = t164 * t138 + t159 * t153;
t171 = qJD(1) * qJD(2);
t146 = t165 * t150 + t160 * t157;
t148 = qJD(4) + t149;
t140 = t159 * t146 - t164 * t148;
t145 = t160 * t150 - t165 * t157;
t137 = -t165 * t142 + t160 * t152;
t163 = cos(qJ(6));
t158 = sin(qJ(6));
t144 = qJD(5) + t145;
t141 = t164 * t146 + t159 * t148;
t139 = qJD(6) + t140;
t135 = t159 * t138 - t164 * t153;
t134 = t163 * t141 + t158 * t144;
t133 = t158 * t141 - t163 * t144;
t132 = t144 * pkin(6) + t172;
t131 = -t141 * pkin(6) + t137;
t1 = [t175, 0, 0, t162 ^ 2 * t175, t162 * t173, t162 * t171, t167 * t171, qJD(2) ^ 2 / 0.2e1, pkin(1) * t173, -t168 * pkin(1) * t162, t150 ^ 2 / 0.2e1, -t150 * t149, t150 * t157, -t149 * t157, t157 ^ 2 / 0.2e1, t154 * t149 + t157 * t169, t154 * t150 - t157 * t170, t146 ^ 2 / 0.2e1, -t146 * t145, t146 * t148, -t145 * t148, t148 ^ 2 / 0.2e1, -t137 * t148 + t153 * t145, -t138 * t148 + t153 * t146, t141 ^ 2 / 0.2e1, -t141 * t140, t141 * t144, -t140 * t144, t144 ^ 2 / 0.2e1, -t135 * t144 + t137 * t140, t137 * t141 - t172 * t144, t134 ^ 2 / 0.2e1, -t134 * t133, t134 * t139, -t133 * t139, t139 ^ 2 / 0.2e1 (t163 * t131 - t158 * t132) * t139 + t135 * t133 -(t158 * t131 + t163 * t132) * t139 + t135 * t134;];
T_reg  = t1;
