% Calculate minimal parameter regressor of potential energy for
% S4RPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% U_reg [1x14]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPPR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:53
% EndTime: 2019-12-31 16:38:54
% DurationCPUTime: 0.05s
% Computational Cost: add. (30->15), mult. (32->22), div. (0->0), fcn. (26->6), ass. (0->11)
t68 = g(3) * (qJ(2) + pkin(4));
t60 = qJ(1) + pkin(6);
t58 = sin(t60);
t59 = cos(t60);
t67 = -g(1) * t58 + g(2) * t59;
t63 = sin(qJ(1));
t65 = cos(qJ(1));
t66 = -g(1) * t65 - g(2) * t63;
t64 = cos(qJ(4));
t62 = sin(qJ(4));
t1 = [0, t66, g(1) * t63 - g(2) * t65, t66 * pkin(1) - t68, g(1) * t59 + g(2) * t58, t67, -g(1) * (t65 * pkin(1) + t59 * pkin(2) + t58 * qJ(3)) - g(2) * (t63 * pkin(1) + t58 * pkin(2) - t59 * qJ(3)) - t68, 0, 0, 0, 0, 0, -g(3) * t64 + t67 * t62, g(3) * t62 + t67 * t64;];
U_reg = t1;
