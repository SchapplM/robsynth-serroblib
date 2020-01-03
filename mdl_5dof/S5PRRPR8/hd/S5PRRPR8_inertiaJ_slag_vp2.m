% Calculate joint inertia matrix for
% S5PRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR8_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR8_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR8_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:19
% EndTime: 2019-12-31 17:42:20
% DurationCPUTime: 0.25s
% Computational Cost: add. (285->81), mult. (587->116), div. (0->0), fcn. (517->8), ass. (0->40)
t36 = cos(qJ(5));
t30 = t36 ^ 2;
t33 = sin(qJ(5));
t49 = t33 ^ 2 + t30;
t58 = mrSges(6,3) * t49;
t57 = m(5) * pkin(3);
t34 = sin(qJ(3));
t35 = sin(qJ(2));
t37 = cos(qJ(3));
t38 = cos(qJ(2));
t18 = -t34 * t35 + t37 * t38;
t19 = t34 * t38 + t37 * t35;
t31 = sin(pkin(9));
t32 = cos(pkin(9));
t6 = -t32 * t18 + t31 * t19;
t56 = t6 ^ 2;
t20 = -t36 * mrSges(6,1) + t33 * mrSges(6,2);
t55 = 0.2e1 * t20;
t53 = pkin(2) * t34;
t26 = t37 * pkin(2) + pkin(3);
t12 = t32 * t26 - t31 * t53;
t52 = t12 * mrSges(5,1);
t13 = t31 * t26 + t32 * t53;
t51 = t13 * mrSges(5,2);
t50 = Ifges(6,5) * t33 + Ifges(6,6) * t36;
t11 = pkin(7) + t13;
t48 = t49 * t11;
t24 = t31 * pkin(3) + pkin(7);
t47 = t49 * t24;
t46 = Ifges(6,2) * t30 + Ifges(4,3) + Ifges(5,3) + (Ifges(6,1) * t33 + 0.2e1 * Ifges(6,4) * t36) * t33;
t45 = t32 * mrSges(5,1) - t31 * mrSges(5,2);
t44 = -mrSges(6,1) * t33 - mrSges(6,2) * t36;
t43 = 0.2e1 * t58;
t42 = (t37 * mrSges(4,1) - t34 * mrSges(4,2)) * pkin(2);
t8 = t31 * t18 + t32 * t19;
t41 = t18 * mrSges(4,1) - t19 * mrSges(4,2) + (-mrSges(5,1) + t20) * t6 + (-mrSges(5,2) + t58) * t8;
t25 = -t32 * pkin(3) - pkin(4);
t10 = -pkin(4) - t12;
t5 = t8 ^ 2;
t1 = [m(2) + m(6) * (t49 * t5 + t56) + m(4) * (t18 ^ 2 + t19 ^ 2) + m(5) * (t5 + t56) + m(3) * (t35 ^ 2 + t38 ^ 2); t38 * mrSges(3,1) - t35 * mrSges(3,2) + m(6) * (t10 * t6 + t8 * t48) + m(5) * (-t12 * t6 + t13 * t8) + m(4) * (t18 * t37 + t19 * t34) * pkin(2) + t41; 0.2e1 * t52 - 0.2e1 * t51 + t10 * t55 + Ifges(3,3) + 0.2e1 * t42 + t11 * t43 + m(6) * (t49 * t11 ^ 2 + t10 ^ 2) + m(5) * (t12 ^ 2 + t13 ^ 2) + m(4) * (t34 ^ 2 + t37 ^ 2) * pkin(2) ^ 2 + t46; m(6) * (t25 * t6 + t8 * t47) + (t31 * t8 - t32 * t6) * t57 + t41; m(6) * (t25 * t10 + t11 * t47) - t51 + t52 + (t10 + t25) * t20 + t42 + (m(5) * (t12 * t32 + t13 * t31) + t45) * pkin(3) + (t47 + t48) * mrSges(6,3) + t46; t25 * t55 + t24 * t43 + m(6) * (t49 * t24 ^ 2 + t25 ^ 2) + t46 + (0.2e1 * t45 + (t31 ^ 2 + t32 ^ 2) * t57) * pkin(3); 0; 0; 0; m(6) * t49 + m(5); t44 * t8; t44 * t11 + t50; t44 * t24 + t50; -t20; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
