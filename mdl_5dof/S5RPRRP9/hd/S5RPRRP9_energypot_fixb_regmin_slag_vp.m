% Calculate minimal parameter regressor of potential energy for
% S5RPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP9_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP9_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:49:10
% EndTime: 2019-12-31 18:49:11
% DurationCPUTime: 0.05s
% Computational Cost: add. (89->31), mult. (78->38), div. (0->0), fcn. (72->8), ass. (0->17)
t146 = pkin(8) + qJ(3);
t149 = sin(qJ(1));
t150 = cos(qJ(1));
t152 = g(1) * t150 + g(2) * t149;
t143 = qJ(4) + t146;
t139 = sin(t143);
t140 = cos(t143);
t142 = cos(t146);
t148 = cos(pkin(8));
t151 = t148 * pkin(2) + pkin(3) * t142 + pkin(4) * t140 + qJ(5) * t139 + pkin(1);
t147 = sin(pkin(8));
t145 = -pkin(7) - pkin(6) - qJ(2);
t141 = sin(t146);
t138 = g(1) * t149 - g(2) * t150;
t136 = -g(3) * t139 - t152 * t140;
t135 = -g(3) * t140 + t152 * t139;
t1 = [0, -t152, t138, -g(3) * t147 - t152 * t148, -g(3) * t148 + t152 * t147, -t138, -g(1) * (t150 * pkin(1) + t149 * qJ(2)) - g(2) * (t149 * pkin(1) - t150 * qJ(2)) - g(3) * pkin(5), 0, 0, 0, 0, 0, -g(3) * t141 - t152 * t142, -g(3) * t142 + t152 * t141, 0, 0, 0, 0, 0, t136, t135, t136, -t138, -t135, -g(3) * (t147 * pkin(2) + pkin(3) * t141 + t139 * pkin(4) - t140 * qJ(5) + pkin(5)) + (-g(1) * t151 - g(2) * t145) * t150 + (g(1) * t145 - g(2) * t151) * t149;];
U_reg = t1;
