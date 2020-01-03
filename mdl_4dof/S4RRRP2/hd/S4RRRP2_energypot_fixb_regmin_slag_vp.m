% Calculate minimal parameter regressor of potential energy for
% S4RRRP2
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
% U_reg [1x15]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRRP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:13:09
% EndTime: 2019-12-31 17:13:09
% DurationCPUTime: 0.03s
% Computational Cost: add. (35->19), mult. (32->23), div. (0->0), fcn. (29->6), ass. (0->12)
t77 = qJ(1) + qJ(2);
t75 = sin(t77);
t76 = cos(t77);
t83 = g(1) * t76 + g(2) * t75;
t82 = cos(qJ(1));
t81 = cos(qJ(3));
t80 = sin(qJ(1));
t79 = sin(qJ(3));
t78 = -qJ(4) - pkin(6);
t74 = t81 * pkin(3) + pkin(2);
t73 = g(1) * t75 - g(2) * t76;
t1 = [0, -g(1) * t82 - g(2) * t80, g(1) * t80 - g(2) * t82, 0, -t83, t73, 0, 0, 0, 0, 0, -g(3) * t79 - t83 * t81, -g(3) * t81 + t83 * t79, -t73, -g(1) * (t82 * pkin(1) + t76 * t74 - t75 * t78) - g(2) * (t80 * pkin(1) + t75 * t74 + t76 * t78) - g(3) * (t79 * pkin(3) + pkin(4) + pkin(5));];
U_reg = t1;
