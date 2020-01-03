% Calculate minimal parameter regressor of potential energy for
% S4RPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% U_reg [1x16]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPPR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:48
% EndTime: 2019-12-31 16:39:48
% DurationCPUTime: 0.03s
% Computational Cost: add. (29->18), mult. (54->28), div. (0->0), fcn. (56->6), ass. (0->14)
t76 = sin(qJ(1));
t78 = cos(qJ(1));
t83 = t78 * pkin(1) + t76 * qJ(2);
t82 = cos(pkin(6));
t81 = sin(pkin(6));
t80 = t76 * pkin(1) - t78 * qJ(2);
t66 = -t76 * t81 - t78 * t82;
t67 = t76 * t82 - t78 * t81;
t79 = g(1) * t66 - g(2) * t67;
t77 = cos(qJ(4));
t75 = sin(qJ(4));
t69 = -g(1) * t78 - g(2) * t76;
t68 = g(1) * t76 - g(2) * t78;
t1 = [0, t69, t68, t69, -t68, -g(3) * pkin(4) - g(1) * t83 - g(2) * t80, t79, -g(1) * t67 - g(2) * t66, -g(1) * (t78 * pkin(2) + t83) - g(2) * (t76 * pkin(2) + t80) - g(3) * (-qJ(3) + pkin(4)), 0, 0, 0, 0, 0, g(3) * t75 + t79 * t77, g(3) * t77 - t79 * t75;];
U_reg = t1;
