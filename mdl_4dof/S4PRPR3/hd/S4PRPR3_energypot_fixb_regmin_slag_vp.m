% Calculate minimal parameter regressor of potential energy for
% S4PRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% U_reg [1x15]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PRPR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:54
% EndTime: 2019-12-31 16:20:54
% DurationCPUTime: 0.05s
% Computational Cost: add. (43->18), mult. (36->22), div. (0->0), fcn. (32->8), ass. (0->11)
t76 = pkin(6) + qJ(2);
t72 = sin(t76);
t74 = cos(t76);
t79 = g(1) * t74 + g(2) * t72;
t78 = cos(pkin(7));
t77 = sin(pkin(7));
t75 = pkin(7) + qJ(4);
t73 = cos(t75);
t71 = sin(t75);
t70 = g(1) * t72 - g(2) * t74;
t1 = [-g(3) * qJ(1), 0, -t79, t70, -g(3) * t77 - t79 * t78, -g(3) * t78 + t79 * t77, -t70, -g(1) * (t74 * pkin(2) + t72 * qJ(3) + cos(pkin(6)) * pkin(1)) - g(2) * (t72 * pkin(2) - t74 * qJ(3) + sin(pkin(6)) * pkin(1)) - g(3) * (pkin(4) + qJ(1)), 0, 0, 0, 0, 0, -g(3) * t71 - t79 * t73, -g(3) * t73 + t79 * t71;];
U_reg = t1;
