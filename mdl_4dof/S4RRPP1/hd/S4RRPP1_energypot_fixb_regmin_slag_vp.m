% Calculate minimal parameter regressor of potential energy for
% S4RRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% 
% Output:
% U_reg [1x10]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRPP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:33:09
% EndTime: 2019-03-08 18:33:09
% DurationCPUTime: 0.05s
% Computational Cost: add. (46->20), mult. (30->25), div. (0->0), fcn. (24->6), ass. (0->12)
t74 = g(3) * (qJ(3) + pkin(5) + pkin(4));
t69 = qJ(1) + qJ(2);
t64 = sin(t69);
t70 = sin(qJ(1));
t73 = t70 * pkin(1) + pkin(2) * t64;
t65 = cos(t69);
t71 = cos(qJ(1));
t72 = t71 * pkin(1) + pkin(2) * t65;
t63 = pkin(6) + t69;
t60 = cos(t63);
t59 = sin(t63);
t1 = [0, -g(1) * t71 - g(2) * t70, g(1) * t70 - g(2) * t71, 0, -g(1) * t65 - g(2) * t64, g(1) * t64 - g(2) * t65, -g(1) * t72 - g(2) * t73 - t74, -g(1) * t60 - g(2) * t59, -g(1) * t59 + g(2) * t60, -g(1) * (t60 * pkin(3) + t59 * qJ(4) + t72) - g(2) * (t59 * pkin(3) - t60 * qJ(4) + t73) - t74;];
U_reg  = t1;
