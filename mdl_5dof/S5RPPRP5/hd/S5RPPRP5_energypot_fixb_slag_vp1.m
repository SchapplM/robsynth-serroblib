% Calculate potential energy for
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRP5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_energypot_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP5_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP5_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:17
% EndTime: 2019-12-31 17:53:17
% DurationCPUTime: 0.23s
% Computational Cost: add. (112->76), mult. (194->91), div. (0->0), fcn. (188->6), ass. (0->28)
t82 = rSges(6,1) + pkin(4);
t63 = sin(pkin(7));
t81 = t63 * pkin(2) + pkin(5);
t66 = sin(qJ(1));
t80 = t63 * t66;
t67 = cos(qJ(4));
t79 = t63 * t67;
t64 = cos(pkin(7));
t78 = t64 * t66;
t68 = cos(qJ(1));
t77 = t64 * t68;
t76 = t68 * pkin(1) + t66 * qJ(2);
t75 = qJ(3) * t63;
t74 = rSges(6,3) + qJ(5);
t60 = t66 * pkin(1);
t73 = pkin(2) * t78 + t66 * t75 + t60;
t72 = pkin(2) * t77 + t68 * t75 + t76;
t71 = pkin(3) * t78 + t68 * pkin(6) + t73;
t70 = pkin(3) * t77 + t72;
t65 = sin(qJ(4));
t47 = t63 * t65 + t64 * t67;
t69 = t63 * pkin(3) - qJ(3) * t64 + t81;
t48 = -t64 * t65 + t79;
t46 = t47 * t68;
t45 = t65 * t77 - t68 * t79;
t44 = t47 * t66;
t43 = t65 * t78 - t66 * t79;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t68 - t66 * rSges(2,2)) + g(2) * (t66 * rSges(2,1) + rSges(2,2) * t68) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t66 * rSges(3,3) + t76) + g(2) * (rSges(3,1) * t78 - rSges(3,2) * t80 + t60) + g(3) * (rSges(3,1) * t63 + rSges(3,2) * t64 + pkin(5)) + (g(1) * (rSges(3,1) * t64 - rSges(3,2) * t63) + g(2) * (-rSges(3,3) - qJ(2))) * t68) - m(4) * (g(1) * (t66 * rSges(4,2) + t72) + g(2) * (rSges(4,1) * t78 + rSges(4,3) * t80 + t73) + g(3) * (rSges(4,1) * t63 + (-rSges(4,3) - qJ(3)) * t64 + t81) + (g(1) * (rSges(4,1) * t64 + rSges(4,3) * t63) + g(2) * (-rSges(4,2) - qJ(2))) * t68) - m(5) * (g(1) * (rSges(5,1) * t46 - rSges(5,2) * t45 + (-rSges(5,3) - pkin(6)) * t66 + t70) + g(2) * (t44 * rSges(5,1) - t43 * rSges(5,2) + (rSges(5,3) - qJ(2)) * t68 + t71) + g(3) * (rSges(5,1) * t48 - rSges(5,2) * t47 + t69)) - m(6) * (g(1) * ((-rSges(6,2) - pkin(6)) * t66 + t82 * t46 + t74 * t45 + t70) + g(2) * ((rSges(6,2) - qJ(2)) * t68 + t82 * t44 + t74 * t43 + t71) + g(3) * (t74 * t47 + t82 * t48 + t69));
U = t1;
