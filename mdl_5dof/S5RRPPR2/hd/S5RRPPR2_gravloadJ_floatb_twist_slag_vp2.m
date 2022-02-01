% Calculate Gravitation load on the joints for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:05:27
% EndTime: 2022-01-20 10:05:28
% DurationCPUTime: 0.30s
% Computational Cost: add. (280->62), mult. (198->72), div. (0->0), fcn. (168->10), ass. (0->35)
t64 = -mrSges(5,3) + mrSges(4,2);
t29 = sin(pkin(9));
t30 = cos(pkin(9));
t63 = -t30 * mrSges(5,1) - mrSges(4,1) + (mrSges(5,2) - mrSges(6,3)) * t29;
t61 = pkin(4) * t30 + pkin(7) * t29;
t60 = m(5) + m(6);
t28 = qJ(1) + qJ(2);
t24 = pkin(8) + t28;
t21 = sin(t24);
t22 = cos(t24);
t25 = sin(t28);
t26 = cos(t28);
t33 = cos(qJ(5));
t31 = sin(qJ(5));
t50 = t30 * t31;
t5 = t21 * t50 + t22 * t33;
t49 = t30 * t33;
t6 = -t21 * t49 + t22 * t31;
t59 = t25 * mrSges(3,1) - t6 * mrSges(6,1) + t26 * mrSges(3,2) - t5 * mrSges(6,2) + t64 * t22 + (-m(6) * (-pkin(3) - t61) - t63) * t21;
t7 = t21 * t33 - t22 * t50;
t8 = t21 * t31 + t22 * t49;
t58 = -t26 * mrSges(3,1) - t8 * mrSges(6,1) + t25 * mrSges(3,2) - t7 * mrSges(6,2) + t64 * t21 + t63 * t22;
t57 = pkin(2) * t25;
t23 = pkin(2) * t26;
t32 = sin(qJ(1));
t54 = t32 * pkin(1);
t34 = cos(qJ(1));
t27 = t34 * pkin(1);
t47 = t22 * pkin(3) + t21 * qJ(4) + t23;
t15 = t22 * qJ(4);
t46 = t15 - t57;
t43 = t61 * t22 + t47;
t42 = -t54 - t57;
t36 = -t21 * pkin(3) + t46;
t1 = [(-t34 * mrSges(2,1) + t32 * mrSges(2,2) - m(3) * t27 - m(4) * (t23 + t27) - m(5) * (t27 + t47) - m(6) * (t27 + t43) + t58) * g(2) + (t32 * mrSges(2,1) + t34 * mrSges(2,2) + m(3) * t54 - m(4) * t42 - m(5) * (t36 - t54) - m(6) * (t15 + t42) + t59) * g(1), (-m(4) * t23 - m(5) * t47 - m(6) * t43 + t58) * g(2) + (m(4) * t57 - m(5) * t36 - m(6) * t46 + t59) * g(1), (-m(4) - t60) * g(3), t60 * (-g(1) * t21 + g(2) * t22), -g(1) * (mrSges(6,1) * t7 - mrSges(6,2) * t8) - g(2) * (-mrSges(6,1) * t5 + mrSges(6,2) * t6) - g(3) * (-mrSges(6,1) * t31 - mrSges(6,2) * t33) * t29];
taug = t1(:);
