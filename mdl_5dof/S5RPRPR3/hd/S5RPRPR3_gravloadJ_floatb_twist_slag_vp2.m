% Calculate Gravitation load on the joints for
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
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:20:36
% EndTime: 2022-01-23 09:20:37
% DurationCPUTime: 0.28s
% Computational Cost: add. (261->59), mult. (185->70), div. (0->0), fcn. (158->10), ass. (0->35)
t62 = -mrSges(5,3) + mrSges(4,2);
t28 = sin(pkin(9));
t29 = cos(pkin(9));
t61 = pkin(4) * t29 + pkin(7) * t28;
t60 = -t29 * mrSges(5,1) - mrSges(4,1) + (mrSges(5,2) - mrSges(6,3)) * t28;
t58 = m(5) + m(6);
t27 = qJ(1) + pkin(8);
t25 = qJ(3) + t27;
t22 = cos(t25);
t57 = t61 * t22;
t21 = sin(t25);
t32 = cos(qJ(5));
t30 = sin(qJ(5));
t47 = t29 * t30;
t5 = t21 * t47 + t22 * t32;
t46 = t29 * t32;
t6 = -t21 * t46 + t22 * t30;
t56 = -t6 * mrSges(6,1) - t5 * mrSges(6,2) + t62 * t22 + (-m(6) * (-pkin(3) - t61) - t60) * t21;
t7 = t21 * t32 - t22 * t47;
t8 = t21 * t30 + t22 * t46;
t55 = -t8 * mrSges(6,1) - t7 * mrSges(6,2) + t62 * t21 + t60 * t22;
t31 = sin(qJ(1));
t54 = pkin(1) * t31;
t53 = pkin(3) * t21;
t33 = cos(qJ(1));
t26 = t33 * pkin(1);
t45 = t22 * pkin(3) + t21 * qJ(4);
t24 = cos(t27);
t44 = pkin(2) * t24 + t26;
t41 = t44 + t45;
t23 = sin(t27);
t40 = -pkin(2) * t23 - t54;
t15 = t22 * qJ(4);
t35 = t15 + t40;
t1 = [(-mrSges(2,1) * t33 + t31 * mrSges(2,2) - m(3) * t26 - t24 * mrSges(3,1) + t23 * mrSges(3,2) - m(4) * t44 - m(5) * t41 - m(6) * (t41 + t57) + t55) * g(2) + (t31 * mrSges(2,1) + mrSges(2,2) * t33 + m(3) * t54 + mrSges(3,1) * t23 + mrSges(3,2) * t24 - m(4) * t40 - m(5) * (t35 - t53) - m(6) * t35 + t56) * g(1), (-m(3) - m(4) - t58) * g(3), (-m(5) * t45 - m(6) * (t45 + t57) + t55) * g(2) + (-m(5) * (t15 - t53) - m(6) * t15 + t56) * g(1), t58 * (-g(1) * t21 + g(2) * t22), -g(1) * (mrSges(6,1) * t7 - mrSges(6,2) * t8) - g(2) * (-mrSges(6,1) * t5 + mrSges(6,2) * t6) - g(3) * (-mrSges(6,1) * t30 - mrSges(6,2) * t32) * t28];
taug = t1(:);
