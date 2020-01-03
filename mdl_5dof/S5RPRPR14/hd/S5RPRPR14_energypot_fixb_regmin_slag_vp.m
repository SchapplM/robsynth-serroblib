% Calculate minimal parameter regressor of potential energy for
% S5RPRPR14
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
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR14_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR14_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:35:15
% EndTime: 2019-12-31 18:35:16
% DurationCPUTime: 0.05s
% Computational Cost: add. (39->30), mult. (61->39), div. (0->0), fcn. (59->8), ass. (0->20)
t141 = sin(qJ(3));
t152 = pkin(3) * t141;
t138 = qJ(3) + pkin(8);
t151 = g(3) * cos(t138);
t140 = sin(qJ(5));
t142 = sin(qJ(1));
t150 = t142 * t140;
t143 = cos(qJ(5));
t149 = t142 * t143;
t145 = cos(qJ(1));
t148 = t145 * t140;
t147 = t145 * t143;
t146 = t145 * pkin(1) + t142 * qJ(2);
t131 = g(1) * t142 - g(2) * t145;
t144 = cos(qJ(3));
t139 = -qJ(4) - pkin(6);
t136 = t142 * pkin(1);
t133 = sin(t138);
t132 = g(1) * t145 + g(2) * t142;
t1 = [0, -t132, t131, t132, -t131, -g(1) * t146 - g(2) * (-t145 * qJ(2) + t136) - g(3) * pkin(5), 0, 0, 0, 0, 0, -g(3) * t144 - t131 * t141, g(3) * t141 - t131 * t144, -t132, -g(1) * (-t145 * t139 + t142 * t152 + t146) - g(2) * (-t142 * t139 + t136 + (-qJ(2) - t152) * t145) - g(3) * (t144 * pkin(3) + pkin(2) + pkin(5)), 0, 0, 0, 0, 0, -g(1) * (t133 * t149 + t148) - g(2) * (-t133 * t147 + t150) - t143 * t151, -g(1) * (-t133 * t150 + t147) - g(2) * (t133 * t148 + t149) + t140 * t151;];
U_reg = t1;
