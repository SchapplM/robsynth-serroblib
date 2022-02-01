% Calculate potential energy for
% S5RRRPR2
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
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPR2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR2_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:30:30
% EndTime: 2022-01-20 11:30:30
% DurationCPUTime: 0.16s
% Computational Cost: add. (134->50), mult. (81->41), div. (0->0), fcn. (50->10), ass. (0->23)
t64 = -m(5) - m(6);
t51 = sin(qJ(5));
t53 = cos(qJ(5));
t63 = -m(6) * pkin(4) - t53 * mrSges(6,1) + t51 * mrSges(6,2) - mrSges(5,1);
t62 = m(6) * pkin(8) - mrSges(5,2) + mrSges(6,3);
t61 = pkin(6) + pkin(5);
t52 = sin(qJ(1));
t48 = t52 * pkin(1);
t54 = cos(qJ(1));
t49 = t54 * pkin(1);
t50 = qJ(1) + qJ(2);
t44 = sin(t50);
t60 = pkin(2) * t44 + t48;
t45 = cos(t50);
t59 = pkin(2) * t45 + t49;
t58 = pkin(7) + t61;
t47 = qJ(3) + t50;
t43 = cos(t47);
t42 = sin(t47);
t41 = pkin(9) + t47;
t36 = cos(t41);
t35 = sin(t41);
t1 = (-m(2) * pkin(5) - m(3) * t61 - m(4) * t58 - t51 * mrSges(6,1) - t53 * mrSges(6,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) + t64 * (qJ(4) + t58)) * g(3) + (-m(3) * t48 - m(4) * t60 - t52 * mrSges(2,1) - t44 * mrSges(3,1) - t42 * mrSges(4,1) - t54 * mrSges(2,2) - t45 * mrSges(3,2) - t43 * mrSges(4,2) - mrSges(1,2) + t64 * (pkin(3) * t42 + t60) + t62 * t36 + t63 * t35) * g(2) + (-m(3) * t49 - m(4) * t59 - t54 * mrSges(2,1) - t45 * mrSges(3,1) - t43 * mrSges(4,1) + t52 * mrSges(2,2) + t44 * mrSges(3,2) + t42 * mrSges(4,2) - mrSges(1,1) + t64 * (pkin(3) * t43 + t59) + t63 * t36 - t62 * t35) * g(1);
U = t1;
