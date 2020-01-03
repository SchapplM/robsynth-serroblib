% Calculate joint inertia matrix for
% S4RRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP2_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP2_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP2_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:12:53
% EndTime: 2019-12-31 17:12:54
% DurationCPUTime: 0.23s
% Computational Cost: add. (168->77), mult. (321->96), div. (0->0), fcn. (182->4), ass. (0->28)
t26 = cos(qJ(3));
t41 = t26 ^ 2;
t24 = sin(qJ(3));
t9 = -t26 * mrSges(5,1) + t24 * mrSges(5,2);
t40 = 0.2e1 * t9;
t27 = cos(qJ(2));
t39 = t27 * pkin(1);
t38 = t24 * mrSges(5,3);
t37 = t24 ^ 2 + t41;
t36 = 0.2e1 * mrSges(5,3);
t15 = -t26 * pkin(3) - pkin(2);
t25 = sin(qJ(2));
t13 = t25 * pkin(1) + pkin(6);
t35 = t37 * t13;
t34 = (Ifges(4,6) + Ifges(5,6)) * t26 + (Ifges(4,5) + Ifges(5,5)) * t24;
t33 = Ifges(3,3) + (Ifges(5,2) + Ifges(4,2)) * t41 + ((Ifges(5,1) + Ifges(4,1)) * t24 + 0.2e1 * (Ifges(4,4) + Ifges(5,4)) * t26) * t24;
t32 = -mrSges(4,1) * t24 - mrSges(4,2) * t26;
t31 = 0.2e1 * t37 * mrSges(4,3);
t30 = (t27 * mrSges(3,1) - t25 * mrSges(3,2)) * pkin(1);
t16 = t26 * qJ(4);
t14 = -pkin(2) - t39;
t11 = t26 * pkin(6) + t16;
t10 = -t26 * mrSges(4,1) + t24 * mrSges(4,2);
t8 = (-qJ(4) - pkin(6)) * t24;
t7 = t15 - t39;
t2 = t26 * t13 + t16;
t1 = (-qJ(4) - t13) * t24;
t3 = [0.2e1 * t14 * t10 + t7 * t40 + Ifges(2,3) + 0.2e1 * t30 + (-t1 * t24 + t2 * t26) * t36 + t13 * t31 + m(5) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) + m(4) * (t37 * t13 ^ 2 + t14 ^ 2) + m(3) * (t25 ^ 2 + t27 ^ 2) * pkin(1) ^ 2 + t33; (t7 + t15) * t9 + (t14 - pkin(2)) * t10 + t30 + m(5) * (t8 * t1 + t11 * t2 + t15 * t7) + m(4) * (-pkin(2) * t14 + pkin(6) * t35) + ((t11 + t2) * t26 + (-t1 - t8) * t24) * mrSges(5,3) + (t37 * pkin(6) + t35) * mrSges(4,3) + t33; -0.2e1 * pkin(2) * t10 + t15 * t40 + (t11 * t26 - t8 * t24) * t36 + pkin(6) * t31 + m(5) * (t11 ^ 2 + t15 ^ 2 + t8 ^ 2) + m(4) * (t37 * pkin(6) ^ 2 + pkin(2) ^ 2) + t33; t1 * mrSges(5,1) - t2 * mrSges(5,2) + t32 * t13 + (m(5) * t1 - t38) * pkin(3) + t34; t8 * mrSges(5,1) - t11 * mrSges(5,2) + t32 * pkin(6) + (m(5) * t8 - t38) * pkin(3) + t34; Ifges(4,3) + Ifges(5,3) + (m(5) * pkin(3) + 0.2e1 * mrSges(5,1)) * pkin(3); m(5) * t7 + t9; m(5) * t15 + t9; 0; m(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t3(1), t3(2), t3(4), t3(7); t3(2), t3(3), t3(5), t3(8); t3(4), t3(5), t3(6), t3(9); t3(7), t3(8), t3(9), t3(10);];
Mq = res;
