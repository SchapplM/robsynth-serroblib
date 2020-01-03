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
% U_reg [1x13]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2019-12-31 16:29:21
% EndTime: 2019-12-31 16:29:21
% DurationCPUTime: 0.05s
% Computational Cost: add. (33->24), mult. (60->36), div. (0->0), fcn. (60->6), ass. (0->16)
t80 = sin(qJ(2));
t88 = g(3) * t80;
t79 = sin(qJ(3));
t82 = cos(qJ(2));
t87 = t79 * t82;
t81 = cos(qJ(3));
t86 = t81 * t82;
t85 = pkin(3) * t79 + pkin(4);
t76 = sin(pkin(6));
t77 = cos(pkin(6));
t84 = g(1) * t77 + g(2) * t76;
t75 = t81 * pkin(3) + pkin(2);
t78 = -qJ(4) - pkin(5);
t83 = t75 * t82 - t78 * t80 + pkin(1);
t74 = -g(3) * t82 + t84 * t80;
t1 = [-g(3) * qJ(1), 0, -t84 * t82 - t88, t74, 0, 0, 0, 0, 0, -g(1) * (t76 * t79 + t77 * t86) - g(2) * (t76 * t86 - t77 * t79) - t81 * t88, -g(1) * (t76 * t81 - t77 * t87) - g(2) * (-t76 * t87 - t77 * t81) + t79 * t88, -t74, -g(3) * (t80 * t75 + t82 * t78 + qJ(1)) + (-g(1) * t83 + g(2) * t85) * t77 + (-g(1) * t85 - g(2) * t83) * t76;];
U_reg = t1;
