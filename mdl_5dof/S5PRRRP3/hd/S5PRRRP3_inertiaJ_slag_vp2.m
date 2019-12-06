% Calculate joint inertia matrix for
% S5PRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP3_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP3_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP3_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:43:32
% EndTime: 2019-12-05 16:43:33
% DurationCPUTime: 0.30s
% Computational Cost: add. (281->92), mult. (559->127), div. (0->0), fcn. (488->4), ass. (0->35)
t32 = cos(qJ(3));
t51 = t32 ^ 2;
t29 = sin(qJ(4));
t31 = cos(qJ(4));
t50 = (t31 * mrSges(5,1) + (-mrSges(5,2) - mrSges(6,2)) * t29) * pkin(3);
t49 = 2 * mrSges(6,1);
t48 = m(6) * pkin(4);
t47 = -pkin(7) - pkin(6);
t30 = sin(qJ(3));
t17 = -t29 * t30 + t31 * t32;
t46 = t17 * pkin(4);
t45 = t31 * pkin(3);
t44 = t29 * t17;
t42 = Ifges(5,3) + Ifges(6,3);
t18 = t29 * t32 + t31 * t30;
t41 = -t17 * mrSges(6,1) + t18 * mrSges(6,2);
t22 = t47 * t30;
t23 = t47 * t32;
t7 = t29 * t22 - t31 * t23;
t40 = t30 ^ 2 + t51;
t25 = -t32 * pkin(3) - pkin(2);
t6 = t31 * t22 + t29 * t23;
t2 = -t18 * qJ(5) + t6;
t38 = m(6) * t2 - t18 * mrSges(6,3);
t37 = -t32 * mrSges(4,1) + t30 * mrSges(4,2);
t11 = t17 * mrSges(5,1);
t36 = -t18 * mrSges(5,2) + t11 - t41;
t3 = t17 * qJ(5) + t7;
t35 = t6 * mrSges(5,1) + t2 * mrSges(6,1) - t7 * mrSges(5,2) - t3 * mrSges(6,2) + (Ifges(5,5) + Ifges(6,5)) * t18 + (Ifges(5,6) + Ifges(6,6)) * t17;
t34 = pkin(3) ^ 2;
t26 = t29 ^ 2 * t34;
t24 = pkin(4) + t45;
t9 = t29 * pkin(3) * t18;
t8 = t25 - t46;
t1 = [m(2) + m(3) + m(4) * t40 + 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (t17 ^ 2 + t18 ^ 2); m(5) * (t6 * t17 + t7 * t18) + m(6) * (t2 * t17 + t3 * t18); Ifges(3,3) + Ifges(4,2) * t51 - 0.2e1 * t25 * t11 + 0.2e1 * t8 * t41 - 0.2e1 * pkin(2) * t37 + m(6) * (t2 ^ 2 + t3 ^ 2 + t8 ^ 2) + m(5) * (t25 ^ 2 + t6 ^ 2 + t7 ^ 2) + m(4) * (t40 * pkin(6) ^ 2 + pkin(2) ^ 2) + (0.2e1 * t25 * mrSges(5,2) - 0.2e1 * t6 * mrSges(5,3) - 0.2e1 * t2 * mrSges(6,3) + (Ifges(5,1) + Ifges(6,1)) * t18) * t18 + 0.2e1 * t40 * pkin(6) * mrSges(4,3) + (Ifges(4,1) * t30 + 0.2e1 * Ifges(4,4) * t32) * t30 + (0.2e1 * t7 * mrSges(5,3) + 0.2e1 * t3 * mrSges(6,3) + 0.2e1 * (Ifges(5,4) + Ifges(6,4)) * t18 + (Ifges(6,2) + Ifges(5,2)) * t17) * t17; m(5) * (t17 * t45 + t9) + m(6) * (t24 * t17 + t9) + t36 - t37; Ifges(4,5) * t30 + Ifges(4,6) * t32 + t38 * t24 + (-t30 * mrSges(4,1) - t32 * mrSges(4,2)) * pkin(6) + (mrSges(6,3) * t44 + (-t31 * t18 + t44) * mrSges(5,3) + m(6) * t29 * t3 + m(5) * (t29 * t7 + t31 * t6)) * pkin(3) + t35; t24 * t49 + Ifges(4,3) + m(6) * (t24 ^ 2 + t26) + m(5) * (t31 ^ 2 * t34 + t26) + 0.2e1 * t50 + t42; m(6) * t46 + t36; pkin(4) * t38 + t35; t24 * t48 + (pkin(4) + t24) * mrSges(6,1) + t50 + t42; (t49 + t48) * pkin(4) + t42; 0; m(6) * t8 + t41; 0; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
