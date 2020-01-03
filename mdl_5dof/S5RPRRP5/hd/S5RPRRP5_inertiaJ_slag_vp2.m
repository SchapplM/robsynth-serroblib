% Calculate joint inertia matrix for
% S5RPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP5_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP5_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP5_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:42
% EndTime: 2019-12-31 18:40:43
% DurationCPUTime: 0.25s
% Computational Cost: add. (248->79), mult. (472->92), div. (0->0), fcn. (287->6), ass. (0->35)
t36 = cos(qJ(4));
t31 = t36 ^ 2;
t34 = sin(qJ(4));
t48 = t34 ^ 2 + t31;
t59 = -m(6) * pkin(4) - mrSges(6,1);
t58 = (mrSges(5,3) + mrSges(6,2)) * t48;
t15 = -t36 * mrSges(6,1) - t34 * mrSges(6,3);
t57 = 0.2e1 * t15;
t33 = cos(pkin(8));
t22 = t33 * pkin(1) + pkin(2);
t35 = sin(qJ(3));
t37 = cos(qJ(3));
t32 = sin(pkin(8));
t52 = pkin(1) * t32;
t10 = t35 * t22 + t37 * t52;
t8 = pkin(7) + t10;
t55 = t48 * pkin(7) * t8;
t54 = t48 * t8 ^ 2;
t53 = m(6) * t34;
t9 = t37 * t22 - t35 * t52;
t51 = t9 * mrSges(4,1);
t50 = t10 * mrSges(4,2);
t24 = t34 * mrSges(6,2);
t49 = t48 * pkin(7) ^ 2;
t47 = qJ(5) * t36;
t16 = -t36 * mrSges(5,1) + t34 * mrSges(5,2);
t45 = t36 * pkin(4) + t34 * qJ(5);
t14 = -pkin(3) - t45;
t44 = Ifges(4,3) + (Ifges(6,3) + Ifges(5,2)) * t31 + ((Ifges(6,1) + Ifges(5,1)) * t34 + 0.2e1 * (Ifges(5,4) - Ifges(6,5)) * t36) * t34;
t42 = mrSges(6,2) * t47 - pkin(4) * t24 + (Ifges(5,6) - Ifges(6,6)) * t36 + (Ifges(6,4) + Ifges(5,5)) * t34;
t41 = 0.2e1 * t58;
t40 = m(6) * t47 + (mrSges(6,3) - mrSges(5,2)) * t36 + (-mrSges(5,1) + t59) * t34;
t7 = -pkin(3) - t9;
t1 = t14 - t9;
t2 = [0.2e1 * t51 - 0.2e1 * t50 + t1 * t57 + 0.2e1 * t7 * t16 + Ifges(2,3) + Ifges(3,3) + t41 * t8 + m(6) * (t1 ^ 2 + t54) + m(5) * (t7 ^ 2 + t54) + m(4) * (t10 ^ 2 + t9 ^ 2) + t44 + (0.2e1 * t33 * mrSges(3,1) - 0.2e1 * t32 * mrSges(3,2) + m(3) * (t32 ^ 2 + t33 ^ 2) * pkin(1)) * pkin(1); 0; m(3) + m(4) + 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t48; t51 - t50 + (t7 - pkin(3)) * t16 + (t1 + t14) * t15 + m(6) * (t14 * t1 + t55) + m(5) * (-pkin(3) * t7 + t55) + t44 + (pkin(7) + t8) * t58; 0; -0.2e1 * pkin(3) * t16 + t14 * t57 + m(6) * (t14 ^ 2 + t49) + m(5) * (pkin(3) ^ 2 + t49) + t41 * pkin(7) + t44; t40 * t8 + t42; m(6) * t45 - t15 - t16; t40 * pkin(7) + t42; Ifges(6,2) + Ifges(5,3) + 0.2e1 * pkin(4) * mrSges(6,1) + 0.2e1 * qJ(5) * mrSges(6,3) + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2); t8 * t53 + t24; -m(6) * t36; pkin(7) * t53 + t24; t59; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
