% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2020-01-03 11:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRP3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP3_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP3_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:46:58
% EndTime: 2020-01-03 11:47:00
% DurationCPUTime: 1.01s
% Computational Cost: add. (3137->119), mult. (5826->155), div. (0->0), fcn. (5622->6), ass. (0->77)
t144 = cos(qJ(4));
t145 = cos(qJ(3));
t95 = sin(qJ(4));
t96 = sin(qJ(3));
t86 = t144 * t145 - t95 * t96;
t115 = -cos(pkin(8)) * pkin(1) - pkin(2);
t88 = -t145 * pkin(3) + t115;
t58 = -t86 * pkin(4) + t88;
t87 = -t144 * t96 - t95 * t145;
t174 = m(6) * t58 - t86 * mrSges(6,1) - t87 * mrSges(6,2);
t108 = (Ifges(5,6) + Ifges(6,6)) * t87 + (Ifges(5,5) + Ifges(6,5)) * t86;
t89 = sin(pkin(8)) * pkin(1) + pkin(6);
t76 = (-pkin(7) - t89) * t96;
t111 = t145 * t89;
t77 = t145 * pkin(7) + t111;
t156 = t144 * t76 - t95 * t77;
t160 = t156 * mrSges(5,2);
t159 = t87 * qJ(5) + t156;
t164 = t159 * mrSges(6,2);
t54 = t144 * t77 + t95 * t76;
t166 = t54 * mrSges(5,1);
t33 = t86 * qJ(5) + t54;
t170 = t33 * mrSges(6,1);
t173 = t108 - t160 - t164 - t166 - t170;
t172 = -t160 / 0.2e1 - t164 / 0.2e1 - t166 / 0.2e1 - t170 / 0.2e1;
t117 = t144 * pkin(3);
t92 = t117 + pkin(4);
t167 = m(6) * (t117 - t92);
t165 = mrSges(5,2) + mrSges(6,2);
t161 = Ifges(6,4) + Ifges(5,4);
t162 = mrSges(6,3) * t159 + t161 * t86;
t158 = -mrSges(5,1) - mrSges(6,1);
t155 = t165 * t144;
t143 = mrSges(6,3) * t86;
t152 = m(6) / 0.2e1;
t154 = (-t143 / 0.2e1 - t33 * t152) * pkin(4) + t172;
t147 = t87 * pkin(4);
t74 = m(6) * t147;
t151 = t74 / 0.2e1;
t148 = pkin(3) * t95;
t146 = t96 * pkin(3);
t134 = t92 * t87;
t13 = m(6) * (t159 * t87 + t33 * t86) + (t86 ^ 2 + t87 ^ 2) * mrSges(6,3);
t124 = t13 * qJD(1);
t78 = t87 * mrSges(6,1);
t55 = t86 * mrSges(6,2) - t78;
t67 = t146 - t147;
t71 = t86 * t148;
t19 = 0.2e1 * (t67 / 0.4e1 - t134 / 0.4e1 - t71 / 0.4e1) * m(6) + t55;
t122 = t19 * qJD(1);
t34 = t55 - t74;
t121 = t34 * qJD(1);
t120 = t87 * t148;
t119 = t92 * t143;
t112 = t145 * mrSges(4,2);
t110 = t161 * t87;
t79 = t87 * mrSges(5,1);
t56 = t86 * mrSges(5,2) - t79;
t109 = -Ifges(5,1) - Ifges(6,1) + Ifges(5,2) + Ifges(6,2);
t107 = t86 * t117;
t102 = -t143 * t159 + t58 * t55 + t88 * t56;
t1 = t145 ^ 2 * Ifges(4,4) + t115 * t112 + t162 * t86 + (t115 * mrSges(4,1) - Ifges(4,4) * t96 + (m(5) * t88 - t86 * mrSges(5,1)) * pkin(3) + (Ifges(4,1) - Ifges(4,2)) * t145) * t96 + (-mrSges(5,2) * t146 + t109 * t86 - t110) * t87 + t102 + t174 * t67;
t105 = t1 * qJD(1);
t2 = (-t174 * pkin(4) - t110) * t87 + (t109 * t87 + t162) * t86 + t102;
t104 = t2 * qJD(1);
t103 = -t165 * t86 + t78 + t79;
t101 = (t71 + t134) * t152;
t100 = t87 * t167;
t18 = t151 + t100 / 0.2e1;
t43 = t155 * pkin(3) + (-t158 - t167) * t148;
t97 = t119 / 0.2e1 - mrSges(6,3) * t107 / 0.2e1 - t33 * t167 / 0.2e1 - t172;
t5 = t97 + t154;
t99 = t5 * qJD(1) + t18 * qJD(2) + t43 * qJD(3);
t22 = t67 * t152 + t101;
t14 = -t100 / 0.2e1 + t151 + t103;
t4 = -t97 + t108 + t154;
t3 = [t1 * qJD(3) + t2 * qJD(4) + t13 * qJD(5), 0, (-mrSges(4,1) * t111 + mrSges(6,3) * t120 - t119 + m(6) * (t148 * t159 - t92 * t33) + m(5) * (-t144 * t54 + t156 * t95) * pkin(3) + Ifges(4,5) * t145 + (t89 * mrSges(4,2) - Ifges(4,6)) * t96 + (t120 - t107) * mrSges(5,3) + t173) * qJD(3) + t4 * qJD(4) + t22 * qJD(5) + t105, t4 * qJD(3) + ((-m(6) * t33 - t143) * pkin(4) + t173) * qJD(4) + t104, t22 * qJD(3) + t124; 0, 0, t14 * qJD(4) + (-t96 * mrSges(4,1) - t112 - t55 - t56 + m(5) * (t87 * t117 + t71) + 0.2e1 * t101) * qJD(3), t14 * qJD(3) + (t74 + t103) * qJD(4), 0; -t5 * qJD(4) - t19 * qJD(5) - t105, -t18 * qJD(4), -t43 * qJD(4), ((-m(6) * pkin(4) + t158) * t95 - t155) * qJD(4) * pkin(3) - t99, -t122; t5 * qJD(3) - t34 * qJD(5) - t104, t18 * qJD(3), t99, 0, -t121; t19 * qJD(3) + t34 * qJD(4) - t124, 0, t122, t121, 0;];
Cq = t3;
