% Calculate minimal parameter regressor of potential energy for
% S5RRPRP8
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
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRP8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:04:26
% EndTime: 2019-12-31 20:04:26
% DurationCPUTime: 0.08s
% Computational Cost: add. (55->34), mult. (106->44), div. (0->0), fcn. (106->6), ass. (0->22)
t152 = sin(qJ(2));
t169 = qJ(3) * t152 + pkin(1);
t153 = sin(qJ(1));
t156 = cos(qJ(1));
t159 = g(1) * t156 + g(2) * t153;
t166 = t152 * pkin(2) + pkin(5);
t151 = sin(qJ(4));
t164 = t152 * t151;
t155 = cos(qJ(2));
t163 = t153 * t155;
t162 = pkin(4) * t164;
t161 = pkin(2) * t163 + t169 * t153;
t160 = t153 * pkin(6) + (pkin(2) * t155 + t169) * t156;
t154 = cos(qJ(4));
t158 = t155 * t151 - t152 * t154;
t157 = t155 * t154 + t164;
t150 = -qJ(5) - pkin(7);
t145 = t154 * pkin(4) + pkin(3);
t140 = g(1) * t153 - g(2) * t156;
t139 = -g(3) * t152 - t159 * t155;
t138 = -g(3) * t155 + t159 * t152;
t1 = [0, -t159, t140, 0, 0, 0, 0, 0, t139, t138, t139, -t140, -t138, -g(1) * t160 - g(2) * (-t156 * pkin(6) + t161) - g(3) * (-t155 * qJ(3) + t166), 0, 0, 0, 0, 0, g(3) * t158 - t159 * t157, g(3) * t157 + t159 * t158, t140, -g(1) * (t153 * t150 + t160) - g(2) * (t145 * t163 + t153 * t162 + t161) - g(3) * (t152 * t145 + (-pkin(4) * t151 - qJ(3)) * t155 + t166) + (-g(1) * (t145 * t155 + t162) - g(2) * (-pkin(6) - t150)) * t156;];
U_reg = t1;
