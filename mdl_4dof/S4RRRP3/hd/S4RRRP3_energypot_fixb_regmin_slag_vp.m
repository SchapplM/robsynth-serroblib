% Calculate minimal parameter regressor of potential energy for
% S4RRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% U_reg [1x17]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRRP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:14
% EndTime: 2019-12-31 17:14:14
% DurationCPUTime: 0.04s
% Computational Cost: add. (48->20), mult. (49->24), div. (0->0), fcn. (46->6), ass. (0->14)
t81 = qJ(1) + qJ(2);
t79 = sin(t81);
t80 = cos(t81);
t88 = g(1) * t80 + g(2) * t79;
t83 = sin(qJ(1));
t85 = cos(qJ(1));
t87 = -g(1) * t85 - g(2) * t83;
t82 = sin(qJ(3));
t84 = cos(qJ(3));
t86 = pkin(3) * t84 + qJ(4) * t82 + pkin(2);
t78 = g(1) * t79 - g(2) * t80;
t77 = -g(3) * t82 - t88 * t84;
t76 = -g(3) * t84 + t88 * t82;
t1 = [0, t87, g(1) * t83 - g(2) * t85, 0, -t88, t78, 0, 0, 0, 0, 0, t77, t76, t77, -t78, -t76, -g(3) * (t82 * pkin(3) - t84 * qJ(4) + pkin(4) + pkin(5)) + t87 * pkin(1) + (g(2) * pkin(6) - g(1) * t86) * t80 + (-g(1) * pkin(6) - g(2) * t86) * t79;];
U_reg = t1;
