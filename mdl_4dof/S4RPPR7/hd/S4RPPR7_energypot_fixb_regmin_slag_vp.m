% Calculate minimal parameter regressor of potential energy for
% S4RPPR7
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
% U_reg [1x17]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPPR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:40
% EndTime: 2019-12-31 16:41:40
% DurationCPUTime: 0.05s
% Computational Cost: add. (30->19), mult. (46->24), div. (0->0), fcn. (40->6), ass. (0->12)
t77 = sin(qJ(1));
t78 = cos(qJ(1));
t80 = t78 * pkin(1) + t77 * qJ(2);
t79 = t77 * pkin(1) - qJ(2) * t78;
t67 = g(1) * t77 - g(2) * t78;
t76 = cos(pkin(6));
t75 = sin(pkin(6));
t74 = pkin(6) + qJ(4);
t70 = cos(t74);
t69 = sin(t74);
t68 = g(1) * t78 + g(2) * t77;
t1 = [0, -t68, t67, t68, -t67, -g(3) * pkin(4) - g(1) * t80 - g(2) * t79, -g(3) * t76 - t67 * t75, g(3) * t75 - t67 * t76, -t68, -g(1) * (qJ(3) * t78 + t80) - g(2) * (t77 * qJ(3) + t79) - g(3) * (pkin(2) + pkin(4)), 0, 0, 0, 0, 0, -g(3) * t70 - t67 * t69, g(3) * t69 - t67 * t70;];
U_reg = t1;
