% Calculate potential energy for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2022-01-23 09:28
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP2_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:27:47
% EndTime: 2022-01-23 09:27:48
% DurationCPUTime: 0.24s
% Computational Cost: add. (130->45), mult. (92->34), div. (0->0), fcn. (61->8), ass. (0->20)
t60 = -m(5) - m(6);
t59 = mrSges(5,2) + mrSges(6,2);
t58 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t57 = -m(4) + t60;
t46 = sin(qJ(4));
t48 = cos(qJ(4));
t56 = pkin(3) * t60 + t59 * t46 + t58 * t48 - mrSges(4,1);
t55 = m(5) * pkin(7) - m(6) * (-qJ(5) - pkin(7)) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t47 = sin(qJ(1));
t41 = t47 * pkin(1);
t49 = cos(qJ(1));
t42 = t49 * pkin(1);
t54 = qJ(2) + pkin(5);
t44 = qJ(1) + pkin(8);
t40 = qJ(3) + t44;
t39 = cos(t44);
t38 = sin(t44);
t36 = cos(t40);
t35 = sin(t40);
t1 = (-m(2) * pkin(5) - m(3) * t54 - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3) - t59 * t48 + t58 * t46 + t57 * (pkin(6) + t54)) * g(3) + (-m(3) * t41 - t47 * mrSges(2,1) - t38 * mrSges(3,1) - t49 * mrSges(2,2) - t39 * mrSges(3,2) - mrSges(1,2) + t57 * (pkin(2) * t38 + t41) + t55 * t36 + t56 * t35) * g(2) + (-m(3) * t42 - t49 * mrSges(2,1) - t39 * mrSges(3,1) + t47 * mrSges(2,2) + t38 * mrSges(3,2) - mrSges(1,1) + t57 * (pkin(2) * t39 + t42) + t56 * t36 - t55 * t35) * g(1);
U = t1;
