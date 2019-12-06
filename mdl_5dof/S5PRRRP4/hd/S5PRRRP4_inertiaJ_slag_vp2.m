% Calculate joint inertia matrix for
% S5PRRRP4
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
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP4_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP4_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP4_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:45:37
% EndTime: 2019-12-05 16:45:39
% DurationCPUTime: 0.33s
% Computational Cost: add. (262->89), mult. (576->113), div. (0->0), fcn. (392->6), ass. (0->40)
t45 = cos(qJ(4));
t41 = t45 ^ 2;
t42 = sin(qJ(4));
t59 = t42 ^ 2 + t41;
t77 = mrSges(5,3) + mrSges(6,2);
t43 = sin(qJ(3));
t31 = t43 * pkin(2) + pkin(7);
t76 = t59 * t31;
t44 = sin(qJ(2));
t46 = cos(qJ(3));
t47 = cos(qJ(2));
t20 = t43 * t47 + t46 * t44;
t75 = t59 * t20;
t74 = -m(6) * pkin(4) - mrSges(6,1);
t73 = t77 * t59;
t18 = t43 * t44 - t46 * t47;
t72 = t18 ^ 2;
t24 = -t45 * mrSges(6,1) - t42 * mrSges(6,3);
t71 = 0.2e1 * t24;
t70 = t76 * t20;
t69 = m(6) * t42;
t68 = t46 * pkin(2);
t34 = t42 * mrSges(6,2);
t63 = t75 * pkin(7);
t62 = t76 * pkin(7);
t61 = t59 * t31 ^ 2;
t60 = t59 * pkin(7) ^ 2;
t58 = qJ(5) * t45;
t21 = -t45 * pkin(4) - t42 * qJ(5) - pkin(3);
t56 = Ifges(4,3) + (Ifges(6,3) + Ifges(5,2)) * t41 + ((Ifges(6,1) + Ifges(5,1)) * t42 + 0.2e1 * (Ifges(5,4) - Ifges(6,5)) * t45) * t42;
t55 = (t46 * mrSges(4,1) - t43 * mrSges(4,2)) * pkin(2);
t25 = -t45 * mrSges(5,1) + t42 * mrSges(5,2);
t53 = -t20 * mrSges(4,2) + (-mrSges(4,1) + t24 + t25) * t18 + t77 * t75;
t52 = mrSges(6,2) * t58 - pkin(4) * t34 + (Ifges(5,6) - Ifges(6,6)) * t45 + (Ifges(6,4) + Ifges(5,5)) * t42;
t51 = 0.2e1 * t73;
t50 = m(6) * t58 + (-mrSges(5,2) + mrSges(6,3)) * t45 + (-mrSges(5,1) + t74) * t42;
t32 = -pkin(3) - t68;
t14 = t20 ^ 2;
t12 = t21 - t68;
t1 = [m(2) + m(3) * (t44 ^ 2 + t47 ^ 2) + m(4) * (t14 + t72) + 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (t59 * t14 + t72); t47 * mrSges(3,1) - t44 * mrSges(3,2) + m(5) * (t32 * t18 + t70) + m(6) * (t12 * t18 + t70) + m(4) * (-t18 * t46 + t20 * t43) * pkin(2) + t53; t12 * t71 + 0.2e1 * t32 * t25 + Ifges(3,3) + 0.2e1 * t55 + m(6) * (t12 ^ 2 + t61) + m(5) * (t32 ^ 2 + t61) + m(4) * (t43 ^ 2 + t46 ^ 2) * pkin(2) ^ 2 + t51 * t31 + t56; m(5) * (-pkin(3) * t18 + t63) + m(6) * (t21 * t18 + t63) + t53; (t32 - pkin(3)) * t25 + (t12 + t21) * t24 + t55 + m(6) * (t21 * t12 + t62) + m(5) * (-pkin(3) * t32 + t62) + t56 + (pkin(7) + t31) * t73; -0.2e1 * pkin(3) * t25 + t21 * t71 + m(6) * (t21 ^ 2 + t60) + m(5) * (pkin(3) ^ 2 + t60) + t51 * pkin(7) + t56; t50 * t20; t50 * t31 + t52; t50 * pkin(7) + t52; Ifges(6,2) + Ifges(5,3) + 0.2e1 * pkin(4) * mrSges(6,1) + 0.2e1 * qJ(5) * mrSges(6,3) + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2); t20 * t69; t31 * t69 + t34; pkin(7) * t69 + t34; t74; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
