% Calculate minimal parameter regressor of potential energy for
% S5PPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
% 
% Output:
% U_reg [1x15]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PPRPR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:29
% EndTime: 2019-12-31 17:33:29
% DurationCPUTime: 0.06s
% Computational Cost: add. (40->22), mult. (71->28), div. (0->0), fcn. (80->6), ass. (0->14)
t85 = cos(qJ(3));
t84 = sin(qJ(3));
t83 = g(3) * qJ(1);
t75 = cos(pkin(7));
t81 = sin(pkin(7));
t82 = t75 * pkin(1) + t81 * qJ(2);
t80 = t81 * pkin(1) - t75 * qJ(2);
t67 = -t75 * t85 - t81 * t84;
t68 = t75 * t84 - t81 * t85;
t79 = g(1) * t68 - g(2) * t67;
t78 = g(1) * t67 + g(2) * t68;
t77 = cos(qJ(5));
t76 = sin(qJ(5));
t1 = [-t83, -g(1) * t82 - g(2) * t80 - t83, 0, t78, t79, -t78, -t79, -g(1) * (pkin(2) * t75 - pkin(3) * t67 + qJ(4) * t68 + t82) - g(2) * (t81 * pkin(2) - t68 * pkin(3) - t67 * qJ(4) + t80) - g(3) * (-pkin(5) + qJ(1)), 0, 0, 0, 0, 0, g(3) * t77 - t79 * t76, -g(3) * t76 - t79 * t77;];
U_reg = t1;
