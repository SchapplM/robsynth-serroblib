% Calculate joint inertia matrix for
% S4RRRP3
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
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP3_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP3_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP3_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:01
% EndTime: 2019-12-31 17:14:02
% DurationCPUTime: 0.23s
% Computational Cost: add. (158->64), mult. (324->75), div. (0->0), fcn. (161->4), ass. (0->28)
t27 = cos(qJ(3));
t24 = t27 ^ 2;
t25 = sin(qJ(3));
t49 = t25 ^ 2 + t24;
t48 = -m(5) * pkin(3) - mrSges(5,1);
t47 = (mrSges(4,3) + mrSges(5,2)) * t49;
t8 = -t27 * mrSges(5,1) - t25 * mrSges(5,3);
t46 = 0.2e1 * t8;
t26 = sin(qJ(2));
t14 = t26 * pkin(1) + pkin(6);
t45 = t49 * pkin(6) * t14;
t44 = m(5) * t25;
t28 = cos(qJ(2));
t42 = t28 * pkin(1);
t17 = t25 * mrSges(5,2);
t41 = t49 * t14 ^ 2;
t40 = t49 * pkin(6) ^ 2;
t38 = qJ(4) * t27;
t5 = -t27 * pkin(3) - t25 * qJ(4) - pkin(2);
t36 = Ifges(3,3) + (Ifges(5,3) + Ifges(4,2)) * t24 + ((Ifges(5,1) + Ifges(4,1)) * t25 + 0.2e1 * (Ifges(4,4) - Ifges(5,5)) * t27) * t25;
t35 = (t28 * mrSges(3,1) - t26 * mrSges(3,2)) * pkin(1);
t33 = mrSges(5,2) * t38 - pkin(3) * t17 + (Ifges(4,6) - Ifges(5,6)) * t27 + (Ifges(5,4) + Ifges(4,5)) * t25;
t32 = 0.2e1 * t47;
t31 = m(5) * t38 + (mrSges(5,3) - mrSges(4,2)) * t27 + (-mrSges(4,1) + t48) * t25;
t15 = -pkin(2) - t42;
t9 = -t27 * mrSges(4,1) + t25 * mrSges(4,2);
t1 = t5 - t42;
t2 = [t1 * t46 + 0.2e1 * t15 * t9 + Ifges(2,3) + 0.2e1 * t35 + m(4) * (t15 ^ 2 + t41) + m(5) * (t1 ^ 2 + t41) + m(3) * (t26 ^ 2 + t28 ^ 2) * pkin(1) ^ 2 + t32 * t14 + t36; (t15 - pkin(2)) * t9 + (t1 + t5) * t8 + t35 + m(4) * (-pkin(2) * t15 + t45) + m(5) * (t5 * t1 + t45) + t36 + (pkin(6) + t14) * t47; -0.2e1 * pkin(2) * t9 + t5 * t46 + m(5) * (t5 ^ 2 + t40) + m(4) * (pkin(2) ^ 2 + t40) + t32 * pkin(6) + t36; t14 * t31 + t33; pkin(6) * t31 + t33; Ifges(5,2) + Ifges(4,3) + 0.2e1 * pkin(3) * mrSges(5,1) + 0.2e1 * qJ(4) * mrSges(5,3) + m(5) * (pkin(3) ^ 2 + qJ(4) ^ 2); t14 * t44 + t17; pkin(6) * t44 + t17; t48; m(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t2(1), t2(2), t2(4), t2(7); t2(2), t2(3), t2(5), t2(8); t2(4), t2(5), t2(6), t2(9); t2(7), t2(8), t2(9), t2(10);];
Mq = res;
