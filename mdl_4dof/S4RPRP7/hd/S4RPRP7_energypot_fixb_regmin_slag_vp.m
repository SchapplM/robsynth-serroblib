% Calculate minimal parameter regressor of potential energy for
% S4RPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% U_reg [1x17]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPRP7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_energypot_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:16
% EndTime: 2019-12-31 16:47:16
% DurationCPUTime: 0.04s
% Computational Cost: add. (30->23), mult. (56->27), div. (0->0), fcn. (50->4), ass. (0->12)
t66 = sin(qJ(3));
t68 = cos(qJ(3));
t73 = pkin(3) * t66 - qJ(4) * t68;
t67 = sin(qJ(1));
t69 = cos(qJ(1));
t71 = t69 * pkin(1) + t67 * qJ(2);
t61 = g(1) * t67 - g(2) * t69;
t64 = t67 * pkin(1);
t62 = g(1) * t69 + g(2) * t67;
t60 = -g(3) * t66 + t61 * t68;
t59 = -g(3) * t68 - t61 * t66;
t1 = [0, -t62, t61, t62, -t61, -g(1) * t71 - g(2) * (-t69 * qJ(2) + t64) - g(3) * pkin(4), 0, 0, 0, 0, 0, t59, -t60, t59, -t62, t60, -g(1) * (t73 * t67 + t71) - g(2) * (t67 * pkin(5) + t64) - g(3) * (t68 * pkin(3) + t66 * qJ(4) + pkin(2) + pkin(4)) + (-g(1) * pkin(5) - g(2) * (-qJ(2) - t73)) * t69;];
U_reg = t1;
