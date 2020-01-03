% Calculate minimal parameter regressor of potential energy for
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRP5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:42
% EndTime: 2019-12-31 17:53:42
% DurationCPUTime: 0.08s
% Computational Cost: add. (74->44), mult. (156->58), div. (0->0), fcn. (165->6), ass. (0->29)
t155 = sin(pkin(7));
t158 = sin(qJ(1));
t156 = cos(pkin(7));
t168 = t156 * t158;
t171 = t158 * t155 * qJ(3) + pkin(2) * t168;
t159 = cos(qJ(4));
t170 = t155 * t159;
t160 = cos(qJ(1));
t169 = t155 * t160;
t167 = t156 * t160;
t166 = t160 * pkin(1) + t158 * qJ(2);
t152 = t158 * pkin(1);
t165 = -t160 * qJ(2) + t152;
t164 = pkin(2) * t167 + qJ(3) * t169 + t166;
t163 = t155 * pkin(2) - t156 * qJ(3) + pkin(5);
t162 = g(1) * t160 + g(2) * t158;
t157 = sin(qJ(4));
t141 = t155 * t157 + t156 * t159;
t137 = t157 * t168 - t158 * t170;
t139 = t157 * t167 - t159 * t169;
t161 = g(1) * t139 + g(2) * t137 + g(3) * t141;
t143 = g(1) * t158 - g(2) * t160;
t142 = -t156 * t157 + t170;
t140 = t141 * t160;
t138 = t141 * t158;
t136 = -g(3) * t155 - t162 * t156;
t135 = -g(3) * t156 + t162 * t155;
t134 = -g(1) * t140 - g(2) * t138 - g(3) * t142;
t1 = [0, -t162, t143, t136, t135, -t143, -g(3) * pkin(5) - g(1) * t166 - g(2) * t165, t136, -t143, -t135, -g(1) * t164 - g(2) * (t165 + t171) - g(3) * t163, 0, 0, 0, 0, 0, t134, t161, t134, t143, -t161, -g(1) * (pkin(3) * t167 + t140 * pkin(4) - t158 * pkin(6) + t139 * qJ(5) + t164) - g(2) * (pkin(3) * t168 + t138 * pkin(4) + t137 * qJ(5) + t152 + (pkin(6) - qJ(2)) * t160 + t171) - g(3) * (t155 * pkin(3) + t142 * pkin(4) + t141 * qJ(5) + t163);];
U_reg = t1;
