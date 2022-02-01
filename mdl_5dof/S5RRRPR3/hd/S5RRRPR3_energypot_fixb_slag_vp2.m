% Calculate potential energy for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:44
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPR3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR3_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:42:23
% EndTime: 2022-01-20 11:42:23
% DurationCPUTime: 0.23s
% Computational Cost: add. (135->51), mult. (104->38), div. (0->0), fcn. (73->10), ass. (0->22)
t64 = -m(3) - m(4);
t65 = -mrSges(2,1) + pkin(1) * (-m(5) - m(6) + t64);
t49 = qJ(3) + pkin(9);
t41 = qJ(5) + t49;
t36 = sin(t41);
t37 = cos(t41);
t54 = cos(qJ(3));
t38 = t54 * pkin(3) + pkin(2);
t39 = sin(t49);
t40 = cos(t49);
t52 = sin(qJ(3));
t63 = -mrSges(3,1) - m(4) * pkin(2) - t54 * mrSges(4,1) + t52 * mrSges(4,2) - m(5) * t38 - t40 * mrSges(5,1) + t39 * mrSges(5,2) - m(6) * (pkin(4) * t40 + t38) - t37 * mrSges(6,1) + t36 * mrSges(6,2);
t51 = -qJ(4) - pkin(7);
t61 = m(4) * pkin(7) - m(5) * t51 + m(6) * (pkin(8) - t51) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t56 = pkin(6) + pkin(5);
t60 = t52 * pkin(3) + t56;
t55 = cos(qJ(1));
t53 = sin(qJ(1));
t50 = qJ(1) + qJ(2);
t43 = cos(t50);
t42 = sin(t50);
t1 = (-mrSges(1,3) - m(2) * pkin(5) - mrSges(2,3) - mrSges(3,3) - t52 * mrSges(4,1) - t54 * mrSges(4,2) - m(5) * t60 - t39 * mrSges(5,1) - t40 * mrSges(5,2) - m(6) * (pkin(4) * t39 + t60) - t36 * mrSges(6,1) - t37 * mrSges(6,2) + t64 * t56) * g(3) + (-mrSges(2,2) * t55 + t63 * t42 + t61 * t43 + t65 * t53 - mrSges(1,2)) * g(2) + (mrSges(2,2) * t53 - t61 * t42 + t63 * t43 + t65 * t55 - mrSges(1,1)) * g(1);
U = t1;
