% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPR3
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
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:36:03
% EndTime: 2020-01-03 11:36:06
% DurationCPUTime: 0.87s
% Computational Cost: add. (2862->162), mult. (5978->233), div. (0->0), fcn. (4842->8), ass. (0->110)
t166 = qJD(1) + qJD(3);
t101 = cos(qJ(5));
t135 = t101 * mrSges(6,2);
t99 = sin(qJ(5));
t141 = t99 * mrSges(6,1);
t124 = t135 + t141;
t96 = sin(pkin(9));
t94 = t96 ^ 2;
t98 = cos(pkin(9));
t95 = t98 ^ 2;
t165 = t94 * t124 + (t94 + t95) * mrSges(5,3);
t100 = sin(qJ(3));
t102 = cos(qJ(3));
t150 = pkin(1) * sin(pkin(8));
t90 = cos(pkin(8)) * pkin(1) + pkin(2);
t79 = t100 * t90 + t102 * t150;
t75 = qJ(4) + t79;
t66 = t94 * t75;
t91 = t94 * qJ(4);
t137 = t66 + t91;
t142 = t98 * t99;
t78 = -t100 * t150 + t102 * t90;
t83 = -t98 * pkin(4) - t96 * pkin(7) - pkin(3);
t58 = -t78 + t83;
t37 = t101 * t58 - t75 * t142;
t132 = t98 * t101;
t38 = t75 * t132 + t99 * t58;
t62 = -qJ(4) * t142 + t101 * t83;
t63 = qJ(4) * t132 + t99 * t83;
t164 = ((-t37 - t62) * t99 + (t38 + t63) * t101) * t98 + t137;
t41 = t101 * t79 - t78 * t142;
t42 = t78 * t132 + t99 * t79;
t139 = t99 * t96;
t130 = mrSges(6,3) * t139;
t143 = t98 * mrSges(6,2);
t81 = -t130 + t143;
t134 = t101 * t96;
t129 = mrSges(6,3) * t134;
t144 = t98 * mrSges(6,1);
t82 = -t129 - t144;
t163 = (-mrSges(4,2) + t165) * t78 + t41 * t82 + t42 * t81;
t123 = mrSges(6,1) * t101 - mrSges(6,2) * t99;
t149 = t101 / 0.2e1;
t156 = t101 ^ 2;
t157 = t99 ^ 2;
t104 = (t156 / 0.2e1 + t157 / 0.2e1) * t94 * mrSges(6,3) + (t98 * t123 / 0.2e1 + t99 * t81 / 0.2e1 + t82 * t149) * t96;
t162 = qJD(2) * t104;
t106 = (-t144 / 0.2e1 + t82 / 0.2e1) * t99 + (-t143 / 0.2e1 - t81 / 0.2e1) * t101;
t161 = qJD(4) * t106;
t140 = t99 * t82;
t109 = t81 * t149 - t140 / 0.2e1 + (-t135 / 0.2e1 - t141 / 0.2e1) * t98;
t160 = t109 * qJD(4);
t59 = (t156 - 0.1e1 + t157) * t98 * t96;
t151 = m(6) * t59;
t127 = t151 / 0.2e1;
t138 = qJD(4) * t127 - qJD(5) * t104;
t154 = -m(6) / 0.2e1;
t126 = qJD(2) * t154;
t159 = -qJD(5) * t106 + t59 * t126;
t158 = qJD(2) * t127 + qJD(5) * t109;
t155 = m(5) / 0.2e1;
t153 = m(6) / 0.2e1;
t20 = (t101 * t42 - t41 * t99 - t78 * t98) * t96;
t152 = m(6) * t20;
t148 = Ifges(6,4) * t99;
t147 = Ifges(6,5) * t98;
t146 = Ifges(6,6) * t98;
t145 = t95 * t75;
t136 = Ifges(6,4) * t101;
t133 = t95 * qJ(4);
t114 = t81 * t132 - t98 * t140 + t165;
t10 = m(5) * (t66 + t145) + m(6) * (t66 + (t101 * t38 - t37 * t99) * t98) + t114;
t131 = qJD(1) * t10;
t128 = t152 / 0.2e1;
t107 = (t136 + (Ifges(6,1) - Ifges(6,2)) * t99) * t96 - t146;
t113 = (-Ifges(6,4) * t139 - t147) * t99;
t28 = t37 * t81;
t29 = t38 * t82;
t33 = t37 * t130;
t117 = t94 * t123;
t45 = t75 * t117;
t4 = -t45 - t28 + t29 - t33 + (t113 + (t38 * mrSges(6,3) + t107) * t101) * t96;
t122 = -t4 * qJD(1) - t162;
t16 = m(6) * (t91 + (t101 * t63 - t62 * t99) * t98) + m(5) * (t91 + t133) + t114;
t105 = ((qJ(4) + t75) * t95 + t137) * t155 + t114;
t108 = t79 * t155 + (t101 * t41 + t99 * t42) * t153;
t5 = t164 * t154 - t105 + t108;
t121 = qJD(1) * t5 - qJD(3) * t16;
t120 = t166 * t104;
t119 = t166 * t106;
t118 = t41 * mrSges(6,1) / 0.2e1 - t42 * mrSges(6,2) / 0.2e1;
t116 = (-Ifges(6,5) * t99 - Ifges(6,6) * t101) * t96;
t43 = t78 * t66;
t84 = -t98 * mrSges(5,1) + t96 * mrSges(5,2);
t3 = (-mrSges(4,1) + t84) * t79 + m(5) * (t78 * t145 + t43 + (-pkin(3) - t78) * t79) + m(6) * (t37 * t41 + t38 * t42 + t43) + t163;
t115 = -t3 * qJD(1) + t20 * t126;
t48 = t62 * t81;
t49 = t63 * t82;
t55 = t62 * t130;
t72 = qJ(4) * t117;
t103 = -(-t147 + (Ifges(6,1) * t101 - t148) * t96) * t139 / 0.2e1 - (-t146 + (-Ifges(6,2) * t99 + t136) * t96) * t134 / 0.2e1 - t98 * t116 / 0.2e1 + (-t38 / 0.2e1 - t63 / 0.2e1) * t129 + t28 / 0.2e1 - t29 / 0.2e1 + t33 / 0.2e1 + t45 / 0.2e1 + t48 / 0.2e1 - t49 / 0.2e1 + t55 / 0.2e1 + t72 / 0.2e1 + (-t99 * (-Ifges(6,2) * t101 - t148) / 0.2e1 + (-Ifges(6,1) * t99 - t136) * t149) * t94;
t1 = t103 - t118;
t8 = -t48 + t49 - t55 - t72 + (t113 + (t63 * mrSges(6,3) + t107) * t101) * t96;
t112 = t1 * qJD(1) - t8 * qJD(3) - t162;
t65 = t78 * t91;
t32 = 0.2e1 * (qJD(1) / 0.4e1 + qJD(3) / 0.4e1) * t151;
t7 = qJD(3) * t128 + t138;
t6 = t164 * t153 + t105 + t108;
t2 = t103 + t118;
t9 = [qJD(3) * t3 + qJD(4) * t10 - qJD(5) * t4, t7, t6 * qJD(4) + t2 * qJD(5) - t115 + (-t79 * mrSges(4,1) + t79 * t84 + 0.2e1 * (-pkin(3) * t79 + t78 * t133 + t65) * t155 + 0.2e1 * (t62 * t41 + t63 * t42 + t65) * t153 + t163) * qJD(3), qJD(3) * t6 + t131 + t158, t2 * qJD(3) + t160 + (-t38 * mrSges(6,1) - t37 * mrSges(6,2) + t116) * qJD(5) + t122; t7, 0, qJD(1) * t128 + t138, t32, -t123 * qJD(5) * t96 - t120; -t5 * qJD(4) + t1 * qJD(5) + t115, -qJD(1) * t152 / 0.2e1 + t138, qJD(4) * t16 - qJD(5) * t8, -t121 + t158, t160 + (-t63 * mrSges(6,1) - t62 * mrSges(6,2) + t116) * qJD(5) + t112; qJD(3) * t5 - t131 + t159, -t32, t121 + t159, 0, -t124 * qJD(5) - t119; -qJD(3) * t1 - t122 + t161, t120, -t112 + t161, t119, 0;];
Cq = t9;
