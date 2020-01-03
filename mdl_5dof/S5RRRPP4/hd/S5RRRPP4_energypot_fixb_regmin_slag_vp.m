% Calculate minimal parameter regressor of potential energy for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% U_reg [1x23]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:55:45
% EndTime: 2019-12-31 20:55:45
% DurationCPUTime: 0.05s
% Computational Cost: add. (86->29), mult. (74->36), div. (0->0), fcn. (68->8), ass. (0->19)
t157 = qJ(2) + qJ(3);
t153 = cos(t157);
t160 = cos(qJ(2));
t145 = t160 * pkin(2) + pkin(3) * t153 + pkin(1);
t156 = -qJ(4) - pkin(7) - pkin(6);
t159 = sin(qJ(1));
t161 = cos(qJ(1));
t166 = t159 * t145 + t161 * t156;
t152 = sin(t157);
t158 = sin(qJ(2));
t165 = t158 * pkin(2) + pkin(3) * t152 + pkin(5);
t164 = t161 * t145 - t159 * t156;
t163 = g(1) * t161 + g(2) * t159;
t151 = pkin(8) + t157;
t148 = sin(t151);
t149 = cos(t151);
t162 = pkin(4) * t149 + qJ(5) * t148;
t146 = g(1) * t159 - g(2) * t161;
t1 = [0, -t163, t146, 0, 0, 0, 0, 0, -g(3) * t158 - t163 * t160, -g(3) * t160 + t163 * t158, 0, 0, 0, 0, 0, -g(3) * t152 - t163 * t153, -g(3) * t153 + t163 * t152, -t146, -g(1) * t164 - g(2) * t166 - g(3) * t165, -g(3) * t148 - t163 * t149, -t146, g(3) * t149 - t163 * t148, -g(1) * (t162 * t161 + t164) - g(2) * (t162 * t159 + t166) - g(3) * (t148 * pkin(4) - t149 * qJ(5) + t165);];
U_reg = t1;
