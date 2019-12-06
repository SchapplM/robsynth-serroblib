% Calculate minimal parameter regressor of potential energy for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% U_reg [1x17]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRPP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:10:08
% EndTime: 2019-12-05 16:10:08
% DurationCPUTime: 0.11s
% Computational Cost: add. (91->44), mult. (129->64), div. (0->0), fcn. (134->8), ass. (0->28)
t150 = -qJ(4) - pkin(6);
t152 = sin(qJ(2));
t166 = -t150 * t152 + pkin(1);
t165 = g(3) * t152;
t148 = sin(pkin(7));
t151 = sin(qJ(3));
t164 = t148 * t151;
t154 = cos(qJ(2));
t163 = t148 * t154;
t149 = cos(pkin(7));
t162 = t149 * t154;
t160 = t151 * t154;
t153 = cos(qJ(3));
t159 = t153 * t154;
t140 = t153 * pkin(3) + pkin(2);
t158 = t152 * t140 + t154 * t150 + qJ(1);
t157 = g(1) * t149 + g(2) * t148;
t156 = pkin(3) * t164 + t148 * pkin(5) + t140 * t162 + t166 * t149;
t155 = t140 * t163 + (-pkin(3) * t151 - pkin(5)) * t149 + t166 * t148;
t147 = qJ(3) + pkin(8);
t142 = cos(t147);
t141 = sin(t147);
t134 = -g(3) * t154 + t157 * t152;
t133 = t148 * t141 + t142 * t162;
t132 = t141 * t162 - t148 * t142;
t131 = -t149 * t141 + t142 * t163;
t130 = t141 * t163 + t149 * t142;
t1 = [-g(3) * qJ(1), 0, -t157 * t154 - t165, t134, 0, 0, 0, 0, 0, -g(1) * (t149 * t159 + t164) - g(2) * (t148 * t159 - t149 * t151) - t153 * t165, -g(1) * (t148 * t153 - t149 * t160) - g(2) * (-t148 * t160 - t149 * t153) + t151 * t165, -t134, -g(1) * t156 - g(2) * t155 - g(3) * t158, -g(1) * t133 - g(2) * t131 - t142 * t165, -t134, -g(1) * t132 - g(2) * t130 - t141 * t165, -g(1) * (t133 * pkin(4) + t132 * qJ(5) + t156) - g(2) * (t131 * pkin(4) + t130 * qJ(5) + t155) - g(3) * ((pkin(4) * t142 + qJ(5) * t141) * t152 + t158);];
U_reg = t1;
