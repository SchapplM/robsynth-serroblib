% Calculate joint inertia matrix for
% S5RPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR8_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR8_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR8_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:05:43
% EndTime: 2019-12-31 19:05:44
% DurationCPUTime: 0.47s
% Computational Cost: add. (537->126), mult. (925->174), div. (0->0), fcn. (778->6), ass. (0->54)
t45 = sin(qJ(5));
t46 = sin(qJ(4));
t48 = cos(qJ(5));
t49 = cos(qJ(4));
t24 = t45 * t49 + t48 * t46;
t47 = sin(qJ(3));
t16 = t24 * t47;
t23 = t45 * t46 - t48 * t49;
t17 = t23 * t47;
t30 = -t49 * mrSges(5,1) + t46 * mrSges(5,2);
t4 = t23 * mrSges(6,1) + t24 * mrSges(6,2);
t50 = cos(qJ(3));
t87 = t49 ^ 2;
t68 = t46 ^ 2 + t87;
t79 = t68 * mrSges(5,3);
t88 = -(mrSges(4,2) - t79) * t47 - (t4 + t30 - mrSges(4,1)) * t50 - (-t16 * t24 - t17 * t23) * mrSges(6,3);
t51 = -pkin(1) - pkin(2);
t28 = -t47 * qJ(2) + t50 * t51;
t26 = pkin(3) - t28;
t75 = t49 * pkin(4);
t18 = t26 + t75;
t86 = t18 * t4;
t36 = -pkin(3) - t75;
t85 = t36 * t4;
t82 = -Ifges(6,1) * t24 ^ 2 - Ifges(5,2) * t87 - Ifges(4,3) - (Ifges(5,1) * t46 + 0.2e1 * Ifges(5,4) * t49) * t46 + (0.2e1 * Ifges(6,4) * t24 - Ifges(6,2) * t23) * t23;
t81 = (t45 * t23 + t48 * t24) * mrSges(6,3);
t78 = Ifges(6,5) * t24 - Ifges(6,6) * t23;
t77 = -0.2e1 * t30;
t76 = -pkin(8) - pkin(7);
t29 = t50 * qJ(2) + t47 * t51;
t27 = -pkin(7) + t29;
t74 = pkin(8) - t27;
t71 = t28 * mrSges(4,1);
t70 = t29 * mrSges(4,2);
t67 = 0.2e1 * mrSges(6,3);
t64 = t68 * t27;
t63 = t68 * t47;
t62 = -t16 * mrSges(6,1) + t17 * mrSges(6,2);
t61 = -t46 * mrSges(5,1) - t49 * mrSges(5,2);
t60 = Ifges(5,5) * t46 + Ifges(5,6) * t49;
t13 = t74 * t46;
t14 = t74 * t49;
t2 = t48 * t13 + t45 * t14;
t3 = t45 * t13 - t48 * t14;
t57 = t2 * mrSges(6,1) - t3 * mrSges(6,2) - t78;
t33 = t76 * t46;
t34 = t76 * t49;
t11 = t48 * t33 + t45 * t34;
t12 = t45 * t33 - t48 * t34;
t56 = t11 * mrSges(6,1) - t12 * mrSges(6,2) + t78;
t54 = (t48 * mrSges(6,1) - t45 * mrSges(6,2)) * pkin(4);
t44 = t50 ^ 2;
t42 = t47 ^ 2;
t1 = [(2 * pkin(1) * mrSges(3,1)) - 0.2e1 * t71 + 0.2e1 * t70 + 0.2e1 * qJ(2) * mrSges(3,3) - 0.2e1 * t86 + t26 * t77 + Ifges(3,2) + Ifges(2,3) + (t2 * t24 + t3 * t23) * t67 - 0.2e1 * t27 * t79 + m(6) * (t18 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(5) * (t68 * t27 ^ 2 + t26 ^ 2) + m(4) * (t28 ^ 2 + t29 ^ 2) + m(3) * ((pkin(1) ^ 2) + qJ(2) ^ 2) - t82; -m(3) * pkin(1) - mrSges(3,1) + m(6) * (-t16 * t2 - t17 * t3 - t50 * t18) + m(5) * (-t50 * t26 + t27 * t63) + m(4) * (t50 * t28 + t47 * t29) - t88; m(3) + m(4) * (t42 + t44) + m(5) * (t68 * t42 + t44) + m(6) * (t16 ^ 2 + t17 ^ 2 + t44); t71 - t70 + t86 - t85 + (pkin(3) + t26) * t30 + m(6) * (t11 * t2 + t12 * t3 + t36 * t18) + m(5) * (-pkin(3) * t26 + pkin(7) * t64) + ((t11 - t2) * t24 + (t12 - t3) * t23) * mrSges(6,3) + (-t68 * pkin(7) + t64) * mrSges(5,3) + t82; m(5) * (pkin(3) * t50 + pkin(7) * t63) + m(6) * (-t11 * t16 - t12 * t17 - t36 * t50) + t88; pkin(3) * t77 + 0.2e1 * t85 + (-t11 * t24 - t12 * t23) * t67 + 0.2e1 * pkin(7) * t79 + m(6) * (t11 ^ 2 + t12 ^ 2 + t36 ^ 2) + m(5) * (t68 * pkin(7) ^ 2 + pkin(3) ^ 2) - t82; t61 * t27 + (m(6) * (t2 * t48 + t3 * t45) + t81) * pkin(4) + t57 - t60; t61 * t47 + m(6) * (-t16 * t48 - t17 * t45) * pkin(4) + t62; t61 * pkin(7) + (m(6) * (t11 * t48 + t12 * t45) - t81) * pkin(4) + t56 + t60; Ifges(5,3) + Ifges(6,3) + m(6) * (t45 ^ 2 + t48 ^ 2) * pkin(4) ^ 2 + 0.2e1 * t54; t57; t62; t56; Ifges(6,3) + t54; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
