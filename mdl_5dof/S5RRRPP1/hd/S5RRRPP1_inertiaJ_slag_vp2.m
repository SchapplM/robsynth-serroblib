% Calculate joint inertia matrix for
% S5RRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP1_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP1_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP1_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:25
% EndTime: 2019-12-31 20:49:26
% DurationCPUTime: 0.43s
% Computational Cost: add. (571->132), mult. (1069->168), div. (0->0), fcn. (929->6), ass. (0->49)
t53 = cos(qJ(3));
t81 = t53 ^ 2;
t79 = (mrSges(5,3) + mrSges(6,2));
t80 = 2 * t79;
t49 = sin(pkin(8));
t50 = cos(pkin(8));
t51 = sin(qJ(3));
t29 = t49 * t51 - t50 * t53;
t30 = t49 * t53 + t50 * t51;
t12 = t29 * mrSges(6,1) - t30 * mrSges(6,3);
t78 = 0.2e1 * t12;
t13 = t29 * mrSges(5,1) + t30 * mrSges(5,2);
t77 = 0.2e1 * t13;
t76 = t49 * pkin(3);
t54 = cos(qJ(2));
t75 = t54 * pkin(1);
t22 = t30 * mrSges(6,2);
t74 = t51 ^ 2 + t81;
t52 = sin(qJ(2));
t41 = t52 * pkin(1) + pkin(7);
t44 = t53 * qJ(4);
t28 = t53 * t41 + t44;
t67 = (-qJ(4) - t41) * t51;
t10 = t50 * t28 + t49 * t67;
t8 = t49 * t28 - t50 * t67;
t73 = t10 ^ 2 + t8 ^ 2;
t36 = t53 * pkin(7) + t44;
t69 = (-qJ(4) - pkin(7)) * t51;
t15 = t49 * t36 - t50 * t69;
t17 = t50 * t36 + t49 * t69;
t72 = t15 ^ 2 + t17 ^ 2;
t71 = t50 * t30 * mrSges(5,3);
t43 = -t53 * pkin(3) - pkin(2);
t70 = t17 * t10 + t15 * t8;
t68 = t74 * t41;
t66 = -t51 * mrSges(4,1) - t53 * mrSges(4,2);
t65 = 0.2e1 * t74 * mrSges(4,3);
t64 = Ifges(4,2) * t81 + Ifges(3,3) + (Ifges(4,1) * t51 + 0.2e1 * Ifges(4,4) * t53) * t51 + (Ifges(6,1) + Ifges(5,1)) * t30 ^ 2 + ((Ifges(6,3) + Ifges(5,2)) * t29 + 0.2e1 * (-Ifges(5,4) + Ifges(6,5)) * t30) * t29;
t63 = (t54 * mrSges(3,1) - t52 * mrSges(3,2)) * pkin(1);
t62 = t12 + t13;
t11 = t29 * pkin(4) - t30 * qJ(5) + t43;
t37 = qJ(5) + t76;
t40 = -t50 * pkin(3) - pkin(4);
t58 = Ifges(4,5) * t51 + Ifges(4,6) * t53 + t40 * t22 + (Ifges(6,4) + Ifges(5,5)) * t30 + (-mrSges(6,2) * t37 - mrSges(5,3) * t76 - Ifges(5,6) + Ifges(6,6)) * t29;
t42 = -pkin(2) - t75;
t35 = -t53 * mrSges(4,1) + t51 * mrSges(4,2);
t34 = t43 - t75;
t6 = t11 - t75;
t1 = [t64 + m(4) * (t74 * t41 ^ 2 + t42 ^ 2) + m(3) * (t52 ^ 2 + t54 ^ 2) * pkin(1) ^ 2 + t41 * t65 + m(5) * (t34 ^ 2 + t73) + m(6) * (t6 ^ 2 + t73) + 0.2e1 * t42 * t35 + t34 * t77 + t6 * t78 + Ifges(2,3) + 0.2e1 * t63 + (-t10 * t29 + t8 * t30) * t80; t64 + t63 + m(4) * (-pkin(2) * t42 + pkin(7) * t68) + (t74 * pkin(7) + t68) * mrSges(4,3) + m(6) * (t11 * t6 + t70) + m(5) * (t43 * t34 + t70) + (t42 - pkin(2)) * t35 + (t43 + t34) * t13 + (t11 + t6) * t12 + t79 * ((t15 + t8) * t30 + (-t10 - t17) * t29); -0.2e1 * pkin(2) * t35 + t11 * t78 + t43 * t77 + pkin(7) * t65 + m(5) * (t43 ^ 2 + t72) + m(6) * (t11 ^ 2 + t72) + m(4) * (t74 * pkin(7) ^ 2 + pkin(2) ^ 2) + t64 + (t15 * t30 - t17 * t29) * t80; t58 + (-t71 + m(5) * (t10 * t49 - t50 * t8)) * pkin(3) + m(6) * (t37 * t10 + t40 * t8) + t66 * t41 - t8 * mrSges(6,1) - t8 * mrSges(5,1) + t10 * mrSges(6,3) - t10 * mrSges(5,2); t58 + (-t71 + m(5) * (-t15 * t50 + t17 * t49)) * pkin(3) + m(6) * (t40 * t15 + t37 * t17) + t66 * pkin(7) - t15 * mrSges(5,1) + t17 * mrSges(6,3) - t17 * mrSges(5,2) - t15 * mrSges(6,1); -0.2e1 * t40 * mrSges(6,1) + 0.2e1 * t37 * mrSges(6,3) + Ifges(6,2) + Ifges(4,3) + Ifges(5,3) + m(6) * (t37 ^ 2 + t40 ^ 2) + (0.2e1 * t50 * mrSges(5,1) - 0.2e1 * t49 * mrSges(5,2) + m(5) * (t49 ^ 2 + t50 ^ 2) * pkin(3)) * pkin(3); m(5) * t34 + m(6) * t6 + t62; m(5) * t43 + m(6) * t11 + t62; 0; m(5) + m(6); m(6) * t8 + t22; m(6) * t15 + t22; m(6) * t40 - mrSges(6,1); 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
