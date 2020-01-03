% Calculate minimal parameter regressor of potential energy for
% S5RPRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR15_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR15_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:37:23
% EndTime: 2019-12-31 18:37:23
% DurationCPUTime: 0.08s
% Computational Cost: add. (53->39), mult. (89->55), div. (0->0), fcn. (91->8), ass. (0->25)
t142 = sin(qJ(3));
t144 = cos(qJ(3));
t158 = pkin(3) * t142 - qJ(4) * t144;
t156 = g(3) * t144;
t139 = pkin(8) + qJ(5);
t134 = sin(t139);
t143 = sin(qJ(1));
t154 = t143 * t134;
t135 = cos(t139);
t153 = t143 * t135;
t140 = sin(pkin(8));
t152 = t143 * t140;
t141 = cos(pkin(8));
t151 = t143 * t141;
t145 = cos(qJ(1));
t150 = t145 * t134;
t149 = t145 * t135;
t148 = t145 * t140;
t147 = t145 * t141;
t146 = t145 * pkin(1) + t143 * qJ(2);
t132 = g(1) * t143 - g(2) * t145;
t137 = t143 * pkin(1);
t133 = g(1) * t145 + g(2) * t143;
t131 = -g(3) * t142 + t132 * t144;
t1 = [0, -t133, t132, t133, -t132, -g(1) * t146 - g(2) * (-t145 * qJ(2) + t137) - g(3) * pkin(5), 0, 0, 0, 0, 0, -t132 * t142 - t156, -t131, -g(1) * (t142 * t151 + t148) - g(2) * (-t142 * t147 + t152) - t141 * t156, -g(1) * (-t142 * t152 + t147) - g(2) * (t142 * t148 + t151) + t140 * t156, t131, -g(1) * (t158 * t143 + t146) - g(2) * (t143 * pkin(6) + t137) - g(3) * (t144 * pkin(3) + t142 * qJ(4) + pkin(2) + pkin(5)) + (-g(1) * pkin(6) - g(2) * (-qJ(2) - t158)) * t145, 0, 0, 0, 0, 0, -g(1) * (t142 * t153 + t150) - g(2) * (-t142 * t149 + t154) - t135 * t156, -g(1) * (-t142 * t154 + t149) - g(2) * (t142 * t150 + t153) + t134 * t156;];
U_reg = t1;
