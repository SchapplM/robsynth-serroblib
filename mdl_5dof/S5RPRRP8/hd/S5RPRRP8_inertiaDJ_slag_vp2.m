% Calculate time derivative of joint inertia matrix for
% S5RPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP8_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP8_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP8_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP8_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:08
% EndTime: 2019-12-31 18:47:10
% DurationCPUTime: 0.57s
% Computational Cost: add. (485->113), mult. (1040->156), div. (0->0), fcn. (547->4), ass. (0->51)
t45 = sin(qJ(4));
t47 = cos(qJ(4));
t71 = t45 ^ 2 + t47 ^ 2;
t46 = sin(qJ(3));
t48 = cos(qJ(3));
t49 = -pkin(1) - pkin(2);
t23 = -t46 * qJ(2) + t48 * t49;
t12 = t48 * qJD(2) + qJD(3) * t23;
t92 = mrSges(4,2) - (mrSges(6,2) + mrSges(5,3)) * t71;
t97 = t92 * t12;
t69 = qJD(3) * t48;
t96 = t71 * t69;
t95 = t71 * t12;
t67 = qJD(4) * t47;
t68 = qJD(4) * t45;
t84 = Ifges(6,5) * t47;
t85 = Ifges(6,5) * t45;
t50 = ((-Ifges(5,4) * t45 + t85 + (-Ifges(5,2) - Ifges(6,3)) * t47) * t45 + (Ifges(5,4) * t47 - t84 + (Ifges(5,1) + Ifges(6,1)) * t45) * t47) * qJD(4) + ((Ifges(6,1) * t47 + t85) * qJD(4) + Ifges(5,1) * t67 - Ifges(5,4) * t68) * t45 - ((Ifges(6,3) * t45 + t84) * qJD(4) - Ifges(5,4) * t67 + Ifges(5,2) * t68) * t47;
t24 = t48 * qJ(2) + t46 * t49;
t13 = t46 * qJD(2) + qJD(3) * t24;
t27 = -t47 * mrSges(5,1) + t45 * mrSges(5,2);
t78 = mrSges(4,1) - t27;
t94 = t78 * t13;
t54 = pkin(4) * t45 - qJ(5) * t47;
t14 = qJD(4) * t54 - t45 * qJD(5);
t91 = -(qJ(5) * mrSges(6,2) + Ifges(5,6) - Ifges(6,6)) * t45 - (pkin(4) * mrSges(6,2) - Ifges(6,4) - Ifges(5,5)) * t47;
t56 = t45 * mrSges(5,1) + t47 * mrSges(5,2);
t16 = t56 * qJD(4);
t88 = -0.2e1 * t16;
t22 = -pkin(7) + t24;
t87 = t95 * t22;
t86 = t95 * pkin(7);
t81 = t46 * t12;
t80 = t48 * t13;
t15 = -mrSges(6,1) * t68 + mrSges(6,3) * t67;
t77 = -t15 + t16;
t72 = t96 * pkin(7);
t70 = qJD(3) * t46;
t66 = qJD(5) * t47;
t26 = t47 * mrSges(6,1) + t45 * mrSges(6,3);
t64 = t26 + t78;
t61 = t22 * t96 + t71 * t81;
t57 = (m(6) * pkin(7) + mrSges(6,2)) * t47;
t55 = -pkin(4) * t47 - qJ(5) * t45;
t25 = -pkin(3) + t55;
t52 = m(6) * t55 - t26 + t27;
t51 = -m(6) * t54 - t45 * mrSges(6,1) + t47 * mrSges(6,3) - t56;
t21 = pkin(3) - t23;
t9 = -t23 - t25;
t1 = t13 - t14;
t2 = [0.2e1 * t1 * t26 + 0.2e1 * t9 * t15 + t21 * t88 + 0.2e1 * t97 + 0.2e1 * m(5) * (t13 * t21 + t87) + 0.2e1 * m(6) * (t1 * t9 + t87) + 0.2e1 * m(4) * (t12 * t24 - t13 * t23) + t50 + 0.2e1 * t94 + 0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * qJD(2); t77 * t48 + (t64 * t46 + t92 * t48) * qJD(3) + m(5) * (t21 * t70 + t61 - t80) + m(6) * (-t1 * t48 + t70 * t9 + t61) + m(4) * (t81 - t80 + (-t23 * t46 + t24 * t48) * qJD(3)); 0.4e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (-0.1e1 + t71) * t46 * t69; (t14 - t1) * t26 + (pkin(3) + t21) * t16 + (t25 - t9) * t15 - t94 + m(6) * (t1 * t25 + t14 * t9 + t86) + m(5) * (-pkin(3) * t13 + t86) - t97 - t50; -t64 * t70 + m(5) * (-pkin(3) * t70 + t72) + m(6) * (t25 * t70 + t72) + (-m(6) * t14 - qJD(3) * t92 - t77) * t48; pkin(3) * t88 - 0.2e1 * t25 * t15 + 0.2e1 * (m(6) * t25 - t26) * t14 + t50; (m(6) * t22 - mrSges(6,2)) * t66 + t51 * t12 + (t22 * t52 - t91) * qJD(4); t51 * t69 + (m(6) * t66 + qJD(4) * t52) * t46; qJD(5) * t57 + (pkin(7) * t52 + t91) * qJD(4); 0.2e1 * (m(6) * qJ(5) + mrSges(6,3)) * qJD(5); -mrSges(6,2) * t67 + (t12 * t45 + t22 * t67) * m(6); (t45 * t69 + t46 * t67) * m(6); qJD(4) * t57; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
