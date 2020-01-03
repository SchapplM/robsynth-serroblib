% Calculate minimal parameter regressor of potential energy for
% S5RRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPP5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:58:27
% EndTime: 2019-12-31 20:58:28
% DurationCPUTime: 0.06s
% Computational Cost: add. (88->31), mult. (95->36), div. (0->0), fcn. (89->6), ass. (0->19)
t143 = qJ(2) + qJ(3);
t140 = sin(t143);
t146 = cos(qJ(2));
t156 = t146 * pkin(2) + qJ(4) * t140 + pkin(1);
t141 = cos(t143);
t145 = sin(qJ(1));
t154 = t141 * t145;
t147 = cos(qJ(1));
t153 = t141 * t147;
t152 = pkin(3) * t153 + t156 * t147;
t148 = -pkin(7) - pkin(6);
t151 = pkin(3) * t154 + t156 * t145 + t147 * t148;
t150 = g(1) * t147 + g(2) * t145;
t144 = sin(qJ(2));
t149 = t144 * pkin(2) + t140 * pkin(3) - t141 * qJ(4) + pkin(5);
t129 = g(1) * t145 - g(2) * t147;
t128 = -g(3) * t140 - t150 * t141;
t127 = -g(3) * t141 + t150 * t140;
t1 = [0, -t150, t129, 0, 0, 0, 0, 0, -g(3) * t144 - t150 * t146, -g(3) * t146 + t150 * t144, 0, 0, 0, 0, 0, t128, t127, t128, -t129, -t127, -g(1) * (-t145 * t148 + t152) - g(2) * t151 - g(3) * t149, t128, -t127, t129, -g(1) * (pkin(4) * t153 + (-qJ(5) - t148) * t145 + t152) - g(2) * (pkin(4) * t154 + t147 * qJ(5) + t151) - g(3) * (t140 * pkin(4) + t149);];
U_reg = t1;
