% Calculate potential energy for
% S5RRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPPR3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR3_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR3_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:26
% EndTime: 2019-12-31 19:26:26
% DurationCPUTime: 0.18s
% Computational Cost: add. (125->48), mult. (83->38), div. (0->0), fcn. (52->8), ass. (0->19)
t58 = m(5) + m(6);
t46 = sin(qJ(5));
t48 = cos(qJ(5));
t57 = -t46 * mrSges(6,1) - t48 * mrSges(6,2) + mrSges(4,2) - mrSges(5,3);
t56 = -m(6) * pkin(7) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t55 = pkin(6) + pkin(5);
t47 = sin(qJ(1));
t42 = t47 * pkin(1);
t49 = cos(qJ(1));
t43 = t49 * pkin(1);
t45 = qJ(1) + qJ(2);
t40 = sin(t45);
t54 = pkin(2) * t40 + t42;
t41 = cos(t45);
t53 = pkin(2) * t41 + t43;
t39 = pkin(8) + t45;
t36 = cos(t39);
t35 = sin(t39);
t1 = (-m(2) * pkin(5) - m(3) * t55 - m(6) * pkin(4) - t48 * mrSges(6,1) + t46 * mrSges(6,2) - mrSges(5,1) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3) + (-m(4) - t58) * (qJ(3) + t55)) * g(3) + (-m(3) * t42 - m(4) * t54 - t47 * mrSges(2,1) - t40 * mrSges(3,1) - t49 * mrSges(2,2) - t41 * mrSges(3,2) - mrSges(1,2) - t58 * (t35 * pkin(3) + t54) + (t58 * qJ(4) - t57) * t36 + t56 * t35) * g(2) + (-m(3) * t43 - m(4) * t53 - t49 * mrSges(2,1) - t41 * mrSges(3,1) + t47 * mrSges(2,2) + t40 * mrSges(3,2) - mrSges(1,1) - t58 * (t36 * pkin(3) + t35 * qJ(4) + t53) + t56 * t36 + t57 * t35) * g(1);
U = t1;
