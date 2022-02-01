% Calculate minimal parameter regressor of potential energy for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x17]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:34
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:34:19
% EndTime: 2022-01-20 10:34:19
% DurationCPUTime: 0.10s
% Computational Cost: add. (46->17), mult. (29->23), div. (0->0), fcn. (26->8), ass. (0->12)
t39 = qJ(1) + qJ(2);
t36 = pkin(9) + qJ(4) + t39;
t34 = sin(t36);
t35 = cos(t36);
t44 = g(1) * t35 + g(2) * t34;
t43 = cos(qJ(1));
t42 = cos(qJ(5));
t41 = sin(qJ(1));
t40 = sin(qJ(5));
t38 = cos(t39);
t37 = sin(t39);
t1 = [0, -g(1) * t43 - g(2) * t41, g(1) * t41 - g(2) * t43, 0, -g(1) * t38 - g(2) * t37, g(1) * t37 - g(2) * t38, -g(1) * (t43 * pkin(1) + pkin(2) * t38) - g(2) * (t41 * pkin(1) + pkin(2) * t37) - g(3) * (qJ(3) + pkin(6) + pkin(5)), 0, -t44, g(1) * t34 - g(2) * t35, 0, 0, 0, 0, 0, -g(3) * t40 - t44 * t42, -g(3) * t42 + t44 * t40;];
U_reg = t1;
