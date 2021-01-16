% Calculate minimal parameter regressor of potential energy for
% S4PRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% U_reg [1x15]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 22:36
% Revision: beb2ba9bd8c5bd556f42a244985f3dab86917626 (2021-01-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PRRP5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 22:36:07
% EndTime: 2021-01-14 22:36:07
% DurationCPUTime: 0.09s
% Computational Cost: add. (41->24), mult. (80->36), div. (0->0), fcn. (84->6), ass. (0->18)
t83 = sin(qJ(2));
t91 = g(3) * t83;
t82 = sin(qJ(3));
t85 = cos(qJ(2));
t90 = t82 * t85;
t84 = cos(qJ(3));
t89 = t84 * t85;
t88 = pkin(3) * t82 + pkin(4);
t79 = sin(pkin(6));
t80 = cos(pkin(6));
t87 = g(1) * t80 + g(2) * t79;
t78 = t84 * pkin(3) + pkin(2);
t81 = qJ(4) + pkin(5);
t86 = t78 * t85 + t81 * t83 + pkin(1);
t77 = -g(3) * t85 + t87 * t83;
t76 = -g(1) * (t79 * t82 + t80 * t89) - g(2) * (t79 * t89 - t80 * t82) - t84 * t91;
t75 = -g(1) * (t79 * t84 - t80 * t90) - g(2) * (-t79 * t90 - t80 * t84) + t82 * t91;
t1 = [-g(3) * qJ(1), 0, -t87 * t85 - t91, t77, 0, 0, 0, 0, 0, t76, t75, t76, t75, -t77, -g(3) * (t83 * t78 - t85 * t81 + qJ(1)) + (-g(1) * t86 + g(2) * t88) * t80 + (-g(1) * t88 - g(2) * t86) * t79;];
U_reg = t1;
