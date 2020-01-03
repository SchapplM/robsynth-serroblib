% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRR8_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR8_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR8_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR8_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR8_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:06
% EndTime: 2019-12-31 18:01:07
% DurationCPUTime: 0.48s
% Computational Cost: add. (1663->101), mult. (2599->139), div. (0->0), fcn. (2266->6), ass. (0->70)
t68 = sin(qJ(5));
t69 = cos(qJ(5));
t97 = t68 ^ 2 + t69 ^ 2;
t115 = -t69 / 0.2e1;
t105 = sin(qJ(4));
t106 = cos(qJ(4));
t93 = sin(pkin(8));
t94 = cos(pkin(8));
t50 = -t105 * t94 - t106 * t93;
t84 = t69 * mrSges(6,1) - t68 * mrSges(6,2);
t37 = t50 * t84;
t89 = Ifges(6,1) * t69 / 0.2e1 - Ifges(6,4) * t68 + Ifges(6,2) * t115;
t114 = t89 * t68;
t70 = -pkin(1) - pkin(2);
t53 = t94 * qJ(2) + t93 * t70;
t87 = t105 * t53;
t76 = -t93 * qJ(2) + t94 * t70;
t51 = -pkin(3) + t76;
t88 = t106 * t51;
t31 = -t87 + t88;
t85 = t97 * t31;
t113 = t84 + mrSges(5,1);
t49 = t105 * t93 - t106 * t94;
t10 = m(6) * (-0.1e1 + t97) * t50 * t49;
t112 = t10 * qJD(2);
t111 = t97 * mrSges(6,3);
t110 = -t31 / 0.2e1;
t86 = t97 * t49;
t108 = m(6) * (pkin(4) * t50 - pkin(7) * t86);
t98 = t69 * mrSges(6,2);
t99 = t68 * mrSges(6,1);
t55 = t98 + t99;
t107 = pkin(4) * t55;
t103 = Ifges(6,4) * t69;
t102 = t49 * t55;
t32 = t105 * t51 + t106 * t53;
t30 = -pkin(7) + t32;
t75 = pkin(4) - t31;
t81 = -t31 * mrSges(5,2) - t113 * t32;
t3 = -m(6) * (t30 * t85 + t75 * t32) + t81 + t111 * t31;
t96 = t3 * qJD(1);
t40 = t49 * mrSges(5,2);
t73 = t50 * mrSges(5,1) - t111 * t49 + t40;
t4 = -m(6) * (-t30 * t86 - t75 * t50) + t37 - m(5) * (t31 * t50 - t49 * t32) - m(4) * (t53 * t94 - t76 * t93) - m(3) * qJ(2) - mrSges(3,3) - t94 * mrSges(4,2) - t93 * mrSges(4,1) + t73;
t95 = t4 * qJD(1);
t57 = -Ifges(6,2) * t68 + t103;
t58 = Ifges(6,1) * t68 + t103;
t90 = t58 / 0.2e1 + t57 / 0.2e1;
t72 = t90 * t69 + t114;
t74 = t75 * t55;
t13 = -t74 + t72;
t92 = t13 * qJD(1);
t79 = t99 / 0.2e1 + t98 / 0.2e1;
t18 = (t55 / 0.2e1 + t79) * t49;
t91 = t18 * qJD(1);
t83 = -Ifges(6,5) * t69 + Ifges(6,6) * t68;
t71 = -t37 / 0.2e1 + m(6) * ((-t75 - t85) * t50 + (-t97 * t30 + t32) * t49) / 0.2e1;
t80 = t108 / 0.2e1 + t37 / 0.2e1;
t1 = -t71 + t73 + t80;
t82 = -t1 * qJD(1) + t112;
t78 = t79 * t49;
t16 = t72 - t107;
t17 = (-t55 / 0.2e1 + t79) * t49;
t6 = (-t87 / 0.2e1 + t88 / 0.2e1 - pkin(4)) * t55 + (mrSges(6,2) * t110 + t90) * t69 + (mrSges(6,1) * t110 + t89) * t68;
t77 = t6 * qJD(1) + t17 * qJD(2) - t16 * qJD(4);
t20 = -t102 / 0.2e1 + t78;
t19 = t102 / 0.2e1 + t78;
t7 = t74 / 0.2e1 + t107 / 0.2e1 - t79 * t31 - t114 + (t58 + t57) * t115;
t2 = t71 + t80;
t5 = [-t4 * qJD(2) - t3 * qJD(4) + t13 * qJD(5), t2 * qJD(4) + t20 * qJD(5) + t112 - t95, 0, -t96 + t2 * qJD(2) + (m(6) * (-pkin(4) * t32 + pkin(7) * t85) + mrSges(6,3) * t85 + t81) * qJD(4) + t7 * qJD(5), t92 + t20 * qJD(2) + t7 * qJD(4) + (-t84 * t30 + t83) * qJD(5); -t1 * qJD(4) - t18 * qJD(5) + t95, t10 * qJD(4), 0, (-mrSges(6,3) * t86 + t113 * t50 + t108 + t40) * qJD(4) + t19 * qJD(5) + t82, t19 * qJD(4) + qJD(5) * t37 - t91; 0, 0, 0, 0, -t55 * qJD(5); t1 * qJD(2) - t6 * qJD(5) + t96, -t17 * qJD(5) - t82, 0, t16 * qJD(5), (-t84 * pkin(7) - t83) * qJD(5) - t77; t18 * qJD(2) + t6 * qJD(4) - t92, t17 * qJD(4) + t91, 0, t77, 0;];
Cq = t5;
