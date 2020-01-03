% Calculate minimal parameter regressor of potential energy for
% S4PRRP4
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
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PRRP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:27:58
% EndTime: 2019-12-31 16:27:58
% DurationCPUTime: 0.04s
% Computational Cost: add. (47->20), mult. (46->23), div. (0->0), fcn. (42->6), ass. (0->11)
t66 = pkin(6) + qJ(2);
t64 = sin(t66);
t65 = cos(t66);
t70 = g(1) * t65 + g(2) * t64;
t67 = sin(qJ(3));
t68 = cos(qJ(3));
t69 = pkin(3) * t68 + qJ(4) * t67 + pkin(2);
t63 = g(1) * t64 - g(2) * t65;
t62 = -g(3) * t67 - t70 * t68;
t61 = -g(3) * t68 + t70 * t67;
t1 = [-g(3) * qJ(1), 0, -t70, t63, 0, 0, 0, 0, 0, t62, t61, t62, -t63, -t61, -g(3) * (t67 * pkin(3) - t68 * qJ(4) + pkin(4) + qJ(1)) + (-g(1) * cos(pkin(6)) - g(2) * sin(pkin(6))) * pkin(1) + (g(2) * pkin(5) - g(1) * t69) * t65 + (-g(1) * pkin(5) - g(2) * t69) * t64;];
U_reg = t1;
