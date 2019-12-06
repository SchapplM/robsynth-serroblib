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
% Datum: 2019-12-05 17:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:51:10
% EndTime: 2019-12-05 17:51:11
% DurationCPUTime: 0.30s
% Computational Cost: add. (152->54), mult. (112->46), div. (0->0), fcn. (85->10), ass. (0->23)
t50 = sin(qJ(5));
t52 = cos(qJ(5));
t67 = -m(6) * pkin(4) - t52 * mrSges(6,1) + t50 * mrSges(6,2) - mrSges(5,1);
t66 = -m(5) - m(6);
t65 = -t50 * mrSges(6,1) - t52 * mrSges(6,2) + mrSges(4,2) - mrSges(5,3);
t64 = m(6) * pkin(7) + mrSges(6,3);
t48 = sin(pkin(9));
t49 = cos(pkin(9));
t63 = t48 * mrSges(5,2) + t67 * t49 - mrSges(4,1);
t51 = sin(qJ(1));
t62 = pkin(1) * t51;
t53 = cos(qJ(1));
t45 = t53 * pkin(1);
t61 = qJ(2) + pkin(5);
t47 = qJ(1) + pkin(8);
t43 = cos(t47);
t60 = pkin(2) * t43 + t45;
t42 = sin(t47);
t58 = -pkin(2) * t42 - t62;
t44 = qJ(3) + t47;
t41 = cos(t44);
t40 = sin(t44);
t1 = (-m(3) * t45 - m(4) * t60 - t53 * mrSges(2,1) - t43 * mrSges(3,1) + t51 * mrSges(2,2) + t42 * mrSges(3,2) - mrSges(1,3) + t66 * (t41 * pkin(3) + t40 * qJ(4) + t60) + (-t64 * t48 + t63) * t41 + t65 * t40) * g(3) + (m(3) * t62 - m(4) * t58 + t51 * mrSges(2,1) + t42 * mrSges(3,1) + t53 * mrSges(2,2) + t43 * mrSges(3,2) - mrSges(1,2) + t66 * (t41 * qJ(4) + t58) + t65 * t41 + (m(5) * pkin(3) - m(6) * (-pkin(7) * t48 - pkin(3)) + t48 * mrSges(6,3) - t63) * t40) * g(2) + (-m(2) * pkin(5) - m(3) * t61 - mrSges(1,1) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3) + t67 * t48 + (-mrSges(5,2) + t64) * t49 + (-m(4) + t66) * (pkin(6) + t61)) * g(1);
U = t1;
