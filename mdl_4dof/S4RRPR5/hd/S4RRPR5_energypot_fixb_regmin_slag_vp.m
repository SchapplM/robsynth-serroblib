% Calculate minimal parameter regressor of potential energy for
% S4RRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% U_reg [1x16]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRPR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:32
% EndTime: 2019-12-31 17:03:32
% DurationCPUTime: 0.03s
% Computational Cost: add. (33->16), mult. (31->21), div. (0->0), fcn. (28->6), ass. (0->10)
t69 = qJ(1) + qJ(2);
t67 = sin(t69);
t68 = cos(t69);
t65 = g(1) * t67 - g(2) * t68;
t73 = cos(qJ(1));
t72 = cos(qJ(4));
t71 = sin(qJ(1));
t70 = sin(qJ(4));
t66 = g(1) * t68 + g(2) * t67;
t1 = [0, -g(1) * t73 - g(2) * t71, g(1) * t71 - g(2) * t73, 0, -t66, t65, t66, -t65, -g(1) * (t73 * pkin(1) + t68 * pkin(2) + t67 * qJ(3)) - g(2) * (t71 * pkin(1) + t67 * pkin(2) - t68 * qJ(3)) - g(3) * (pkin(5) + pkin(4)), 0, 0, 0, 0, 0, -g(3) * t72 - t65 * t70, g(3) * t70 - t65 * t72;];
U_reg = t1;
