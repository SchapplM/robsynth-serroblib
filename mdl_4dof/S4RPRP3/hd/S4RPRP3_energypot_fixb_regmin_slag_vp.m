% Calculate minimal parameter regressor of potential energy for
% S4RPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% U_reg [1x13]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPRP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:47
% EndTime: 2019-12-31 16:42:47
% DurationCPUTime: 0.03s
% Computational Cost: add. (32->18), mult. (33->25), div. (0->0), fcn. (27->6), ass. (0->13)
t86 = qJ(2) + pkin(4);
t78 = qJ(1) + pkin(6);
t76 = sin(t78);
t77 = cos(t78);
t85 = g(1) * t77 + g(2) * t76;
t81 = sin(qJ(1));
t83 = cos(qJ(1));
t84 = -g(1) * t83 - g(2) * t81;
t82 = cos(qJ(3));
t80 = sin(qJ(3));
t79 = -qJ(4) - pkin(5);
t75 = t82 * pkin(3) + pkin(2);
t1 = [0, t84, g(1) * t81 - g(2) * t83, t84 * pkin(1) - g(3) * t86, 0, 0, 0, 0, 0, -g(3) * t80 - t85 * t82, -g(3) * t82 + t85 * t80, -g(1) * t76 + g(2) * t77, -g(1) * (t83 * pkin(1) + t77 * t75 - t76 * t79) - g(2) * (t81 * pkin(1) + t76 * t75 + t77 * t79) - g(3) * (t80 * pkin(3) + t86);];
U_reg = t1;
