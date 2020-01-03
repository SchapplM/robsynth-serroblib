% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR5_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR5_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:30
% EndTime: 2019-12-31 17:56:31
% DurationCPUTime: 0.39s
% Computational Cost: add. (1050->89), mult. (1733->120), div. (0->0), fcn. (1201->6), ass. (0->56)
t59 = sin(qJ(5));
t61 = cos(qJ(5));
t69 = -Ifges(6,4) * t59 + (Ifges(6,1) - Ifges(6,2)) * t61;
t104 = t69 * t59;
t62 = cos(qJ(4));
t58 = t61 ^ 2;
t85 = t59 ^ 2 + t58;
t78 = t85 * t62;
t102 = -t58 * Ifges(6,4) - t104;
t49 = -t61 * mrSges(6,1) + t59 * mrSges(6,2);
t60 = sin(qJ(4));
t46 = t60 * t49;
t77 = t62 * mrSges(5,2) - mrSges(6,3) * t78;
t98 = -t60 * mrSges(5,1) + t46 - t77;
t53 = -cos(pkin(8)) * pkin(1) - pkin(2) - pkin(3);
t56 = sin(pkin(8)) * pkin(1) + qJ(3);
t28 = t62 * t53 - t60 * t56;
t97 = -t28 / 0.2e1;
t96 = m(6) * (-0.1e1 + t85) * t62 * t60;
t89 = t61 * mrSges(6,2);
t90 = t59 * mrSges(6,1);
t50 = t89 + t90;
t95 = pkin(4) * t50;
t92 = Ifges(6,4) * t61;
t26 = pkin(4) - t28;
t91 = t26 * t50;
t87 = t62 * t50;
t29 = t60 * t53 + t62 * t56;
t27 = -pkin(7) + t29;
t79 = t85 * t28;
t64 = mrSges(6,3) * t79 - t28 * mrSges(5,2) + (-mrSges(5,1) + t49) * t29;
t1 = -m(6) * (t26 * t29 + t27 * t79) + t64;
t84 = t1 * qJD(1);
t6 = -mrSges(4,3) - m(6) * (t26 * t60 + t27 * t78) - m(5) * (-t28 * t60 + t29 * t62) - m(4) * t56 + t98;
t83 = t6 * qJD(1);
t10 = -t91 - t102;
t82 = t10 * qJD(1);
t71 = -t89 / 0.2e1 - t90 / 0.2e1;
t18 = (-t50 / 0.2e1 + t71) * t62;
t81 = t18 * qJD(1);
t80 = qJD(3) * t96;
t72 = -Ifges(6,5) * t61 + Ifges(6,6) * t59;
t63 = -t46 / 0.2e1 + m(6) * ((t85 * t27 - t29) * t62 + (t26 + t79) * t60) / 0.2e1;
t65 = m(6) * (-pkin(4) * t60 + pkin(7) * t78);
t2 = (mrSges(5,1) - t49 / 0.2e1) * t60 - t65 / 0.2e1 + t63 + t77;
t70 = t2 * qJD(1) + t80;
t67 = t71 * t62;
t12 = t95 + t102;
t17 = (t50 / 0.2e1 + t71) * t62;
t4 = (mrSges(6,2) * t97 + t92) * t61 + (-t26 / 0.2e1 - pkin(4) / 0.2e1) * t50 + (mrSges(6,1) * t97 + t69) * t59;
t66 = t4 * qJD(1) + t17 * qJD(3) + t12 * qJD(4);
t20 = t87 / 0.2e1 + t67;
t19 = -t87 / 0.2e1 + t67;
t5 = t91 / 0.2e1 + t95 / 0.2e1 + t71 * t28 - t92 * t61 - t104;
t3 = t65 / 0.2e1 + t46 / 0.2e1 + t63;
t7 = [-t6 * qJD(3) - t1 * qJD(4) + t10 * qJD(5), 0, t3 * qJD(4) + t20 * qJD(5) + t80 - t83, -t84 + t3 * qJD(3) + (m(6) * (-pkin(4) * t29 + pkin(7) * t79) + t64) * qJD(4) + t5 * qJD(5), t82 + t20 * qJD(3) + t5 * qJD(4) + (t49 * t27 + t72) * qJD(5); 0, 0, 0, 0, t50 * qJD(5); t2 * qJD(4) - t18 * qJD(5) + t83, 0, qJD(4) * t96, (t65 + t98) * qJD(4) + t19 * qJD(5) + t70, t19 * qJD(4) + qJD(5) * t46 - t81; -t2 * qJD(3) - t4 * qJD(5) + t84, 0, -t17 * qJD(5) - t70, -t12 * qJD(5), (t49 * pkin(7) - t72) * qJD(5) - t66; t18 * qJD(3) + t4 * qJD(4) - t82, 0, t17 * qJD(4) + t81, t66, 0;];
Cq = t7;
