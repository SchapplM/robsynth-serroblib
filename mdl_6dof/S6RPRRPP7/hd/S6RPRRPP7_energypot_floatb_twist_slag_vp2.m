% Calculate potential energy for
% S6RPRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPP7_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRPP7_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP7_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP7_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP7_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:50:01
% EndTime: 2019-03-09 04:50:01
% DurationCPUTime: 0.57s
% Computational Cost: add. (164->77), mult. (250->67), div. (0->0), fcn. (230->6), ass. (0->36)
t56 = -mrSges(4,2) + mrSges(6,2) + mrSges(5,3);
t53 = -m(6) - m(7);
t55 = -m(5) + t53;
t54 = -m(1) - m(2);
t52 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t51 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t50 = m(7) * qJ(6) + mrSges(7,3);
t21 = sin(qJ(3));
t24 = cos(qJ(3));
t49 = -t21 * mrSges(4,1) + t56 * t24 + mrSges(2,2) - mrSges(3,3);
t48 = -m(7) * pkin(5) - mrSges(5,1) - mrSges(6,1) - mrSges(7,1);
t47 = pkin(3) * t21;
t20 = sin(qJ(4));
t25 = cos(qJ(1));
t46 = t20 * t25;
t22 = sin(qJ(1));
t45 = t22 * t20;
t23 = cos(qJ(4));
t44 = t22 * t23;
t43 = t22 * t24;
t40 = t25 * t23;
t19 = pkin(6) + r_base(3);
t39 = pkin(8) * t43;
t38 = t22 * pkin(1) + r_base(2);
t37 = pkin(2) + t19;
t35 = t22 * pkin(7) + t38;
t34 = t25 * pkin(1) + t22 * qJ(2) + r_base(1);
t33 = t25 * t24 * pkin(8) + t35;
t32 = t25 * pkin(7) + t34;
t29 = t22 * t47 + t32;
t3 = t21 * t45 - t40;
t4 = t21 * t44 + t46;
t26 = t4 * pkin(4) + t3 * qJ(5) + t29;
t6 = -t21 * t40 + t45;
t5 = t21 * t46 + t44;
t1 = (-m(1) * r_base(3) - m(4) * t37 - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * t19 + (t50 - t56) * t21 + t55 * (t24 * pkin(3) + t21 * pkin(8) + t37) + (t53 * (pkin(4) * t23 + qJ(5) * t20) + t51 * t20 + t48 * t23 - mrSges(4,1)) * t24) * g(3) + (-m(3) * t38 - m(4) * t35 - m(5) * t33 - mrSges(1,2) + t54 * r_base(2) + t53 * (t6 * pkin(4) - t5 * qJ(5) + t33) + t48 * t6 - t51 * t5 + t52 * t22 + (t50 * t24 + t55 * (-qJ(2) - t47) + (m(3) + m(4)) * qJ(2) - t49) * t25) * g(2) + (-mrSges(1,1) - m(3) * t34 - m(4) * t32 - m(5) * (t29 - t39) - m(6) * (t26 - t39) - m(7) * t26 - (m(7) * (-pkin(8) + qJ(6)) + mrSges(7,3)) * t43 + t54 * r_base(1) + t48 * t4 + t51 * t3 + t52 * t25 + t49 * t22) * g(1);
U  = t1;
