% Calculate potential energy for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% m_mdh [6x1]
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
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR3_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:35:52
% EndTime: 2020-01-03 11:35:53
% DurationCPUTime: 0.36s
% Computational Cost: add. (152->52), mult. (112->44), div. (0->0), fcn. (85->10), ass. (0->23)
t65 = m(5) + m(6);
t47 = sin(qJ(5));
t49 = cos(qJ(5));
t64 = -m(6) * pkin(4) - t49 * mrSges(6,1) + t47 * mrSges(6,2) - mrSges(5,1);
t62 = -m(4) - t65;
t61 = m(6) * pkin(7) + mrSges(6,3);
t45 = sin(pkin(9));
t46 = cos(pkin(9));
t60 = t45 * mrSges(5,2) + t64 * t46 - mrSges(4,1);
t59 = t47 * mrSges(6,1) + t49 * mrSges(6,2) + t65 * qJ(4) - mrSges(4,2) + mrSges(5,3);
t50 = cos(qJ(1));
t58 = pkin(1) * t50;
t48 = sin(qJ(1));
t42 = t48 * pkin(1);
t57 = qJ(2) + pkin(5);
t44 = qJ(1) + pkin(8);
t39 = sin(t44);
t56 = pkin(2) * t39 + t42;
t41 = qJ(3) + t44;
t40 = cos(t44);
t38 = cos(t41);
t37 = sin(t41);
t1 = (m(3) * t58 + t50 * mrSges(2,1) + t40 * mrSges(3,1) - t48 * mrSges(2,2) - t39 * mrSges(3,2) - mrSges(1,3) + t62 * (-pkin(2) * t40 - t58) + (m(5) * pkin(3) - m(6) * (-pkin(7) * t45 - pkin(3)) + t45 * mrSges(6,3) - t60) * t38 + t59 * t37) * g(3) + (-m(3) * t42 - m(4) * t56 - t48 * mrSges(2,1) - t39 * mrSges(3,1) - t50 * mrSges(2,2) - t40 * mrSges(3,2) - mrSges(1,2) - t65 * (t37 * pkin(3) + t56) + t59 * t38 + (-t61 * t45 + t60) * t37) * g(2) + (-m(2) * pkin(5) - m(3) * t57 - mrSges(1,1) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3) + t64 * t45 + (-mrSges(5,2) + t61) * t46 + t62 * (pkin(6) + t57)) * g(1);
U = t1;
