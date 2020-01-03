% Calculate minimal parameter regressor of potential energy for
% S4PRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% U_reg [1x14]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PRRR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:45
% EndTime: 2019-12-31 16:33:45
% DurationCPUTime: 0.03s
% Computational Cost: add. (29->15), mult. (41->25), div. (0->0), fcn. (44->8), ass. (0->16)
t80 = qJ(2) + qJ(3);
t78 = sin(t80);
t92 = g(3) * t78;
t81 = sin(pkin(7));
t83 = sin(qJ(4));
t91 = t81 * t83;
t85 = cos(qJ(4));
t90 = t81 * t85;
t82 = cos(pkin(7));
t89 = t82 * t83;
t88 = t82 * t85;
t87 = g(1) * t82 + g(2) * t81;
t86 = cos(qJ(2));
t84 = sin(qJ(2));
t79 = cos(t80);
t1 = [-g(3) * qJ(1), 0, -g(3) * t84 - t87 * t86, -g(3) * t86 + t87 * t84, 0, -t87 * t79 - t92, -g(3) * t79 + t87 * t78, 0, 0, 0, 0, 0, -g(1) * (t79 * t88 + t91) - g(2) * (t79 * t90 - t89) - t85 * t92, -g(1) * (-t79 * t89 + t90) - g(2) * (-t79 * t91 - t88) + t83 * t92;];
U_reg = t1;
