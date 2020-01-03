% Calculate minimal parameter regressor of potential energy for
% S5RRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% U_reg [1x23]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRP10_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP10_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:11:05
% EndTime: 2019-12-31 20:11:05
% DurationCPUTime: 0.09s
% Computational Cost: add. (54->40), mult. (102->50), div. (0->0), fcn. (100->6), ass. (0->25)
t151 = sin(qJ(2));
t169 = qJ(3) * t151 + pkin(1);
t150 = sin(qJ(4));
t168 = pkin(4) * t150;
t154 = cos(qJ(2));
t167 = g(3) * t154;
t166 = t151 * pkin(2) + pkin(5);
t152 = sin(qJ(1));
t164 = t152 * t150;
t153 = cos(qJ(4));
t163 = t152 * t153;
t162 = t152 * t154;
t155 = cos(qJ(1));
t161 = t155 * t150;
t160 = t155 * t153;
t159 = t151 * t164;
t158 = pkin(2) * t162 + t169 * t152;
t157 = t152 * pkin(6) + (pkin(2) * t154 + t169) * t155;
t156 = g(1) * t155 + g(2) * t152;
t149 = -qJ(5) - pkin(7);
t144 = t153 * pkin(4) + pkin(3);
t139 = g(1) * t152 - g(2) * t155;
t138 = g(3) * t151 + t156 * t154;
t137 = t156 * t151 - t167;
t1 = [0, -t156, t139, 0, 0, 0, 0, 0, -t138, t137, -t139, t138, -t137, -g(1) * t157 - g(2) * (-t155 * pkin(6) + t158) - g(3) * (-t154 * qJ(3) + t166), 0, 0, 0, 0, 0, -g(1) * (t151 * t161 + t163) - g(2) * (t159 - t160) + t150 * t167, -g(1) * (t151 * t160 - t164) - g(2) * (t151 * t163 + t161) + t153 * t167, -t138, -g(1) * (t152 * t144 + t157) - g(2) * (pkin(4) * t159 - t149 * t162 + t158) - g(3) * (-t151 * t149 + (-qJ(3) - t168) * t154 + t166) + (-g(1) * (-t149 * t154 + t151 * t168) - g(2) * (-pkin(6) - t144)) * t155;];
U_reg = t1;
