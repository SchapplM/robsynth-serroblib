% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% m [6x1]
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
% Datum: 2022-01-23 09:23
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:22:27
% EndTime: 2022-01-23 09:22:30
% DurationCPUTime: 1.24s
% Computational Cost: add. (4621->107), mult. (8453->149), div. (0->0), fcn. (9237->8), ass. (0->61)
t115 = cos(pkin(9));
t88 = sin(pkin(9));
t91 = sin(qJ(3));
t93 = cos(qJ(3));
t79 = t115 * t93 - t88 * t91;
t80 = -t115 * t91 - t88 * t93;
t90 = sin(qJ(5));
t92 = cos(qJ(5));
t106 = t92 * t79 + t80 * t90;
t63 = t79 * t90 - t80 * t92;
t164 = t63 * mrSges(6,1) + t106 * mrSges(6,2);
t162 = qJD(5) * t164;
t85 = sin(pkin(8)) * pkin(1) + pkin(6);
t116 = qJ(4) + t85;
t77 = t116 * t91;
t78 = t116 * t93;
t147 = -t115 * t77 - t88 * t78;
t149 = pkin(7) * t80 + t147;
t44 = -t115 * t78 + t88 * t77;
t42 = pkin(7) * t79 - t44;
t153 = t149 * t92 - t42 * t90;
t31 = t149 * t90 + t42 * t92;
t5 = -t31 * mrSges(6,1) - t153 * mrSges(6,2) + Ifges(6,5) * t106 - Ifges(6,6) * t63;
t165 = t5 * qJD(5);
t135 = t91 * pkin(3);
t159 = m(5) * t135;
t156 = -t80 * mrSges(5,1) + t79 * mrSges(5,2);
t134 = Ifges(6,4) * t63;
t138 = t63 / 0.2e1;
t139 = -t63 / 0.2e1;
t140 = t106 / 0.2e1;
t110 = -cos(pkin(8)) * pkin(1) - pkin(2);
t81 = -pkin(3) * t93 + t110;
t64 = -t79 * pkin(4) + t81;
t3 = t64 * t164 + (0.2e1 * Ifges(6,4) * t106 + (Ifges(6,1) - Ifges(6,2)) * t63) * t140 + (Ifges(6,2) * t106 + t134) * t139 + (Ifges(6,1) * t106 - t134) * t138;
t141 = m(5) * pkin(3);
t137 = pkin(3) * t88;
t117 = t3 * qJD(1);
t95 = (t115 * t80 + t79 * t88) * t141;
t108 = t115 * pkin(3);
t86 = t108 + pkin(4);
t72 = -t90 * t137 + t86 * t92;
t73 = t92 * t137 + t86 * t90;
t99 = m(6) * (t106 * t73 - t63 * t72);
t94 = t99 / 0.2e1 + t95 / 0.2e1;
t97 = t164 + t156;
t66 = -pkin(4) * t80 + t135;
t98 = t159 / 0.2e1 + m(6) * t66 / 0.2e1;
t13 = -t94 + t97 + t98;
t114 = t13 * qJD(1);
t112 = t164 * qJD(1);
t27 = -mrSges(6,1) * t73 - mrSges(6,2) * t72;
t111 = t27 * qJD(5);
t1 = (t110 * mrSges(4,1) - Ifges(4,4) * t91) * t91 + (-mrSges(5,2) * t135 - Ifges(5,4) * t80) * t80 + (Ifges(4,4) * t93 + t110 * mrSges(4,2) + (Ifges(4,1) - Ifges(4,2)) * t91) * t93 + (-mrSges(5,1) * t135 + Ifges(5,4) * t79 + (-Ifges(5,1) + Ifges(5,2)) * t80) * t79 + (t156 + t159) * t81 + t3 + (m(6) * t64 - mrSges(6,1) * t106 + t63 * mrSges(6,2)) * t66;
t102 = t1 * qJD(1);
t8 = (t106 ^ 2 + t63 ^ 2) * mrSges(6,3) + (t79 ^ 2 + t80 ^ 2) * mrSges(5,3) + m(6) * (t106 * t31 - t153 * t63) + m(5) * (t147 * t80 - t44 * t79);
t100 = qJD(1) * t8;
t4 = (t139 + t138) * Ifges(6,6) + (t140 - t106 / 0.2e1) * Ifges(6,5);
t96 = t4 * qJD(1) + t27 * qJD(3);
t18 = t94 + t98;
t2 = [qJD(3) * t1 + qJD(4) * t8 + qJD(5) * t3, 0, (m(6) * (t153 * t73 - t31 * t72) - t147 * mrSges(5,2) + t44 * mrSges(5,1) + (t115 * t44 + t147 * t88) * t141 + Ifges(5,5) * t79 + Ifges(5,6) * t80 + Ifges(4,5) * t93 - Ifges(4,6) * t91 + (-t93 * mrSges(4,1) + t91 * mrSges(4,2)) * t85 + (-t106 * t72 - t63 * t73) * mrSges(6,3) + (-t79 * t108 + t80 * t137) * mrSges(5,3) + t5) * qJD(3) + t18 * qJD(4) + t165 + t102, qJD(3) * t18 + t100, t5 * qJD(3) + t117 + t165; 0, 0, (-t91 * mrSges(4,1) - t93 * mrSges(4,2) + t95 - t97 + t99) * qJD(3) - t162, 0, -qJD(3) * t164 - t162; -qJD(4) * t13 + qJD(5) * t4 - t102, 0, t111, -t114, t96 + t111; qJD(3) * t13 - t100 + t162, 0, t114, 0, t112; -t4 * qJD(3) - qJD(4) * t164 - t117, 0, -t96, -t112, 0;];
Cq = t2;
