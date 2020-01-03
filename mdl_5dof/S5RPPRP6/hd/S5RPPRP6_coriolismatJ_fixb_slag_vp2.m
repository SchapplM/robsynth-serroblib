% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRP6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP6_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP6_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP6_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP6_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:54:59
% EndTime: 2019-12-31 17:55:00
% DurationCPUTime: 0.38s
% Computational Cost: add. (1391->84), mult. (2521->105), div. (0->0), fcn. (2424->4), ass. (0->52)
t58 = sin(pkin(7));
t61 = sin(qJ(4));
t59 = cos(pkin(7));
t89 = cos(qJ(4));
t73 = t89 * t59;
t46 = t61 * t58 - t73;
t88 = t61 * t59;
t47 = t89 * t58 + t88;
t29 = t47 * mrSges(6,1) + t46 * mrSges(6,3);
t52 = t58 * pkin(3) + qJ(2);
t68 = -pkin(4) * t47 - t46 * qJ(5);
t20 = -t68 + t52;
t92 = m(6) * t20;
t99 = t29 + t92;
t98 = Ifges(5,4) - Ifges(6,5);
t97 = t46 ^ 2 + t47 ^ 2;
t54 = m(6) * qJ(5) + mrSges(6,3);
t96 = 0.2e1 * m(6);
t91 = pkin(4) * t46;
t60 = -pkin(1) - qJ(3);
t90 = -pkin(6) + t60;
t86 = t58 ^ 2 + t59 ^ 2;
t40 = t46 * mrSges(5,1);
t41 = t47 * mrSges(6,3);
t83 = t47 * qJ(5);
t67 = -t83 + t91;
t1 = t20 * t41 - t52 * t40 + (-t52 * mrSges(5,2) + t98 * t47) * t47 + (-t20 * mrSges(6,1) + (Ifges(5,1) + Ifges(6,1) - Ifges(5,2) - Ifges(6,3)) * t47 - t98 * t46) * t46 - t99 * t67;
t85 = t1 * qJD(1);
t51 = t90 * t58;
t26 = t61 * t51 - t90 * t73;
t27 = t89 * t51 + t90 * t88;
t74 = m(4) * t86;
t4 = t86 * mrSges(4,3) - t60 * t74 + t97 * (mrSges(5,3) + mrSges(6,2)) + (m(5) + m(6)) * (-t26 * t46 - t27 * t47);
t84 = t4 * qJD(1);
t65 = -t74 / 0.2e1 - 0.2e1 * (m(5) / 0.4e1 + m(6) / 0.4e1) * t97;
t69 = -m(4) / 0.2e1 - m(5) / 0.2e1 - m(6) / 0.2e1;
t8 = t65 + t69;
t82 = t8 * qJD(1);
t9 = -t46 * mrSges(6,1) - t47 * mrSges(5,2) - t40 + t41 + (-t67 / 0.4e1 - t91 / 0.4e1 + t83 / 0.4e1) * t96;
t81 = t9 * qJD(1);
t80 = qJD(4) * t47;
t66 = t47 * mrSges(5,1) - t46 * mrSges(5,2) + t29;
t10 = t58 * mrSges(4,1) + t59 * mrSges(4,2) + mrSges(3,3) + (m(4) + m(3)) * qJ(2) + m(5) * t52 + t92 + t66;
t79 = t10 * qJD(1);
t12 = t99 * t46;
t77 = t12 * qJD(1);
t22 = m(6) * t46;
t76 = t22 * qJD(1);
t75 = t54 * qJD(4);
t14 = m(6) * t27 - t47 * mrSges(6,2);
t7 = t65 - t69;
t2 = [t10 * qJD(2) + t4 * qJD(3) + t1 * qJD(4) + t12 * qJD(5), t7 * qJD(3) + t79, t7 * qJD(2) + t84, t85 + t14 * qJD(5) + (pkin(4) * mrSges(6,2) - Ifges(6,4) - Ifges(5,5)) * t80 + ((-m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1)) * t27 + (mrSges(5,2) - t54) * t26 + (qJ(5) * mrSges(6,2) + Ifges(5,6) - Ifges(6,6)) * t46) * qJD(4), t14 * qJD(4) + t77; t8 * qJD(3) - t79, 0, t82, -t66 * qJD(4) + (t68 * qJD(4) / 0.2e1 + t47 * qJD(5) / 0.2e1) * t96, m(6) * t80; -t8 * qJD(2) + t9 * qJD(4) + t22 * qJD(5) - t84, -t82, 0, t81, t76; -t9 * qJD(3) - t85, 0, -t81, t54 * qJD(5), t75; -t22 * qJD(3) - t77, 0, -t76, -t75, 0;];
Cq = t2;
