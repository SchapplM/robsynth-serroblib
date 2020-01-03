% Calculate joint inertia matrix for
% S5RRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP2_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP2_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP2_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:51:29
% EndTime: 2019-12-31 20:51:31
% DurationCPUTime: 0.36s
% Computational Cost: add. (308->123), mult. (569->134), div. (0->0), fcn. (309->4), ass. (0->44)
t42 = cos(qJ(3));
t39 = t42 ^ 2;
t73 = -Ifges(6,4) - Ifges(5,5);
t40 = sin(qJ(3));
t72 = t40 ^ 2 + t39;
t71 = t42 * pkin(3) + t40 * qJ(4);
t70 = -m(5) * pkin(3) - mrSges(5,1);
t69 = (mrSges(4,3) + mrSges(5,2)) * t72;
t68 = -2 * mrSges(6,3);
t14 = -t42 * mrSges(5,1) - t40 * mrSges(5,3);
t67 = 0.2e1 * t14;
t15 = t42 * mrSges(6,1) + t40 * mrSges(6,2);
t66 = 0.2e1 * t15;
t64 = pkin(7) - qJ(5);
t41 = sin(qJ(2));
t23 = t41 * pkin(1) + pkin(7);
t63 = t72 * pkin(7) * t23;
t62 = (Ifges(6,2) + Ifges(5,3)) * t42 + t73 * t40;
t61 = t72 * t23 ^ 2;
t60 = t72 * pkin(7) ^ 2;
t58 = qJ(4) * t42;
t57 = -qJ(5) + t23;
t56 = t40 * t68;
t55 = pkin(2) + t71;
t43 = cos(qJ(2));
t24 = -t43 * pkin(1) - pkin(2);
t54 = Ifges(4,2) * t39 + Ifges(3,3) + ((Ifges(6,1) + Ifges(5,1) + Ifges(4,1)) * t40 + ((2 * Ifges(4,4)) + t73) * t42) * t40;
t2 = t24 - t71;
t52 = (t43 * mrSges(3,1) - t41 * mrSges(3,2)) * pkin(1);
t50 = 0.2e1 * t69;
t49 = m(5) * t58 + (mrSges(5,3) - mrSges(4,2)) * t42 + (-mrSges(4,1) + t70) * t40;
t44 = -pkin(3) - pkin(4);
t48 = mrSges(5,2) * t58 + (-qJ(4) * mrSges(6,3) + Ifges(4,6) - Ifges(5,6) + Ifges(6,6)) * t42 + (-pkin(3) * mrSges(5,2) - t44 * mrSges(6,3) + Ifges(5,4) + Ifges(4,5) - Ifges(6,5)) * t40;
t45 = qJ(4) ^ 2;
t36 = t42 * pkin(4);
t28 = t40 * mrSges(5,2);
t17 = t64 * t42;
t16 = -t42 * mrSges(4,1) + t40 * mrSges(4,2);
t13 = t64 * t40;
t5 = t57 * t42;
t4 = t57 * t40;
t3 = t36 + t55;
t1 = -t2 + t36;
t6 = [t4 * t56 + t1 * t66 + t2 * t67 + 0.2e1 * t24 * t16 + Ifges(2,3) + 0.2e1 * t52 + (t5 * t68 + t62) * t42 + t50 * t23 + m(5) * (t2 ^ 2 + t61) + m(6) * (t1 ^ 2 + t4 ^ 2 + t5 ^ 2) + m(4) * (t24 ^ 2 + t61) + m(3) * (t41 ^ 2 + t43 ^ 2) * pkin(1) ^ 2 + t54; t62 * t42 + (-pkin(2) + t24) * t16 + (t3 + t1) * t15 + (-t55 + t2) * t14 + t52 + m(5) * (-t2 * t55 + t63) + m(6) * (t3 * t1 + t13 * t4 + t17 * t5) + m(4) * (-pkin(2) * t24 + t63) + ((-t17 - t5) * t42 + (-t13 - t4) * t40) * mrSges(6,3) + t54 + (pkin(7) + t23) * t69; t13 * t56 - 0.2e1 * pkin(2) * t16 - t55 * t67 + t3 * t66 + (t17 * t68 + t62) * t42 + m(5) * (t55 ^ 2 + t60) + m(6) * (t13 ^ 2 + t17 ^ 2 + t3 ^ 2) + m(4) * (pkin(2) ^ 2 + t60) + t50 * pkin(7) + t54; -t4 * mrSges(6,1) + t5 * mrSges(6,2) + m(6) * (qJ(4) * t5 + t44 * t4) + t49 * t23 + t48; m(6) * (qJ(4) * t17 + t44 * t13) - t13 * mrSges(6,1) + t17 * mrSges(6,2) + t49 * pkin(7) + t48; 0.2e1 * pkin(3) * mrSges(5,1) - 0.2e1 * t44 * mrSges(6,1) + Ifges(5,2) + Ifges(4,3) + Ifges(6,3) + 0.2e1 * (mrSges(5,3) + mrSges(6,2)) * qJ(4) + m(5) * (pkin(3) ^ 2 + t45) + m(6) * (t44 ^ 2 + t45); t28 + m(6) * t4 + (m(5) * t23 - mrSges(6,3)) * t40; t28 + m(6) * t13 + (m(5) * pkin(7) - mrSges(6,3)) * t40; m(6) * t44 - mrSges(6,1) + t70; m(5) + m(6); m(6) * t1 + t15; m(6) * t3 + t15; 0; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t6(1), t6(2), t6(4), t6(7), t6(11); t6(2), t6(3), t6(5), t6(8), t6(12); t6(4), t6(5), t6(6), t6(9), t6(13); t6(7), t6(8), t6(9), t6(10), t6(14); t6(11), t6(12), t6(13), t6(14), t6(15);];
Mq = res;
