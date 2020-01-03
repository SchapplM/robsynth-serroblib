% Calculate minimal parameter regressor of potential energy for
% S4RRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% U_reg [1x14]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRPR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:33
% EndTime: 2019-12-31 17:01:33
% DurationCPUTime: 0.05s
% Computational Cost: add. (28->15), mult. (25->21), div. (0->0), fcn. (22->8), ass. (0->10)
t79 = qJ(1) + qJ(2);
t76 = pkin(7) + t79;
t84 = g(1) * cos(t76) + g(2) * sin(t76);
t83 = cos(qJ(1));
t82 = cos(qJ(4));
t81 = sin(qJ(1));
t80 = sin(qJ(4));
t78 = cos(t79);
t77 = sin(t79);
t1 = [0, -g(1) * t83 - g(2) * t81, g(1) * t81 - g(2) * t83, 0, -g(1) * t78 - g(2) * t77, g(1) * t77 - g(2) * t78, -g(1) * (t83 * pkin(1) + pkin(2) * t78) - g(2) * (pkin(1) * t81 + pkin(2) * t77) - g(3) * (qJ(3) + pkin(5) + pkin(4)), 0, 0, 0, 0, 0, -g(3) * t80 - t84 * t82, -g(3) * t82 + t84 * t80;];
U_reg = t1;
