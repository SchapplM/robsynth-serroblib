% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PPRRR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR2_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR2_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:14
% EndTime: 2019-12-05 15:14:16
% DurationCPUTime: 0.81s
% Computational Cost: add. (2229->92), mult. (5342->150), div. (0->0), fcn. (5716->8), ass. (0->66)
t91 = sin(qJ(5));
t92 = sin(qJ(4));
t94 = cos(qJ(5));
t95 = cos(qJ(4));
t104 = t91 * t92 - t94 * t95;
t105 = t91 * t95 + t92 * t94;
t119 = cos(pkin(9));
t143 = cos(qJ(3));
t90 = sin(pkin(9));
t93 = sin(qJ(3));
t74 = -t119 * t143 + t90 * t93;
t46 = t105 * t74;
t48 = t104 * t74;
t26 = -t104 * t46 + t105 * t48;
t163 = m(6) * t26;
t147 = -pkin(7) - pkin(6);
t81 = t147 * t92;
t82 = t147 * t95;
t109 = t81 * t94 + t82 * t91;
t62 = t81 * t91 - t82 * t94;
t11 = -t62 * mrSges(6,1) - t109 * mrSges(6,2) - Ifges(6,5) * t104 - Ifges(6,6) * t105;
t162 = t11 * qJD(5);
t145 = pkin(4) * t92;
t161 = m(6) * t145;
t52 = mrSges(6,1) * t105 - mrSges(6,2) * t104;
t158 = t52 * qJD(5);
t75 = t119 * t93 + t143 * t90;
t45 = t104 * t75;
t47 = t105 * t75;
t110 = mrSges(6,1) * t45 + t47 * mrSges(6,2);
t156 = qJD(5) * t110;
t86 = -pkin(4) * t95 - pkin(3);
t7 = -Ifges(6,4) * t105 ^ 2 - (-Ifges(6,4) * t104 - (-Ifges(6,1) + Ifges(6,2)) * t105) * t104 + t86 * t52;
t107 = -t95 * mrSges(5,1) + t92 * mrSges(5,2);
t138 = t74 * t52;
t152 = t138 / 0.2e1;
t88 = t92 ^ 2;
t149 = m(6) / 0.2e1;
t146 = pkin(4) * t91;
t144 = pkin(4) * t94;
t135 = t92 * mrSges(5,1);
t131 = t95 * mrSges(5,2);
t121 = -t95 ^ 2 - t88;
t116 = qJD(1) * t163;
t115 = -t138 / 0.2e1;
t111 = t121 * t74;
t106 = t131 + t135;
t80 = (mrSges(6,1) * t91 + mrSges(6,2) * t94) * pkin(4);
t103 = t80 * qJD(4);
t102 = t46 * mrSges(6,1) / 0.2e1 - t48 * mrSges(6,2) / 0.2e1;
t58 = t74 * t75;
t8 = m(6) * (-t45 * t48 - t46 * t47 + t58) + m(5) * (t111 * t75 + t58);
t100 = -t8 * qJD(1) - qJD(2) * t163 / 0.2e1;
t96 = (t46 * t94 + t48 * t91) * pkin(4) * t149 + t102;
t97 = t74 * t161;
t1 = t115 - t97 / 0.2e1 + t96;
t53 = mrSges(6,1) * t104 + mrSges(6,2) * t105;
t3 = t86 * t161 + t53 * t145 - t88 * Ifges(5,4) - pkin(3) * t135 + (-pkin(3) * mrSges(5,2) + Ifges(5,4) * t95 + (Ifges(5,1) - Ifges(5,2)) * t92) * t95 + t7;
t99 = t1 * qJD(1) - t3 * qJD(3);
t5 = t115 + t102;
t98 = qJD(1) * t5 - qJD(3) * t7;
t76 = t80 * qJD(5);
t6 = t102 + t152;
t4 = t26 * qJD(3) * t149;
t2 = t97 / 0.2e1 + t96 + (t106 / 0.2e1 + t131 / 0.2e1 + t135 / 0.2e1) * t74 + t152;
t9 = [t8 * qJD(3), t4, t2 * qJD(4) + t6 * qJD(5) - t100 + ((-t104 * t48 - t105 * t46) * mrSges(6,3) + (-mrSges(4,1) + t107 + t53) * t75 + (mrSges(5,3) * t121 + mrSges(4,2)) * t74 + 0.2e1 * (t109 * t46 + t62 * t48 + t86 * t75) * t149 + m(5) * (-pkin(3) * t75 + pkin(6) * t111)) * qJD(3), t2 * qJD(3) + (m(6) * (t144 * t45 - t146 * t47) + t75 * t107 + t110) * qJD(4) + t156, qJD(3) * t6 + qJD(4) * t110 + t156; t4, 0, t116 / 0.2e1, -t158 + (-t106 - t52 + (-t104 * t146 - t105 * t144) * m(6)) * qJD(4), -qJD(4) * t52 - t158; -qJD(4) * t1 - qJD(5) * t5 + t100, -t116 / 0.2e1, qJD(4) * t3 + qJD(5) * t7, t162 - t99 + (Ifges(5,5) * t95 - Ifges(5,6) * t92 + (m(6) * (t109 * t91 - t62 * t94) + (t104 * t94 - t105 * t91) * mrSges(6,3)) * pkin(4) + t107 * pkin(6) + t11) * qJD(4), t11 * qJD(4) + t162 - t98; t1 * qJD(3), 0, t99, -t76, -t103 - t76; t5 * qJD(3), 0, t98, t103, 0;];
Cq = t9;
