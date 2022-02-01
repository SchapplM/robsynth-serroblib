% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2022-01-20 09:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:48:38
% EndTime: 2022-01-20 09:48:42
% DurationCPUTime: 1.37s
% Computational Cost: add. (5014->127), mult. (9606->177), div. (0->0), fcn. (8868->8), ass. (0->93)
t121 = cos(qJ(4));
t229 = t121 ^ 2;
t119 = sin(qJ(4));
t190 = sin(qJ(5));
t191 = cos(qJ(5));
t100 = -t191 * t119 - t190 * t121;
t110 = cos(pkin(9)) * pkin(1) + pkin(2);
t120 = sin(qJ(3));
t122 = cos(qJ(3));
t189 = pkin(1) * sin(pkin(9));
t87 = t110 * t122 - t120 * t189;
t49 = t100 * t87;
t99 = -t190 * t119 + t191 * t121;
t50 = t99 * t87;
t24 = -t100 * t50 + t49 * t99;
t228 = m(6) * t24;
t106 = t119 * mrSges(5,1) + t121 * mrSges(5,2);
t197 = m(6) * pkin(4);
t208 = t49 * mrSges(6,1) / 0.2e1 - t50 * mrSges(6,2) / 0.2e1;
t227 = (t190 * t50 + t191 * t49) * t197 / 0.2e1 - t106 * t87 / 0.2e1 + t208;
t194 = -pkin(8) - pkin(7);
t107 = t194 * t119;
t108 = t194 * t121;
t143 = t191 * t107 + t190 * t108;
t75 = -t190 * t107 + t191 * t108;
t226 = t75 * mrSges(6,1) - t143 * mrSges(6,2);
t88 = t120 * t110 + t122 * t189;
t85 = pkin(7) + t88;
t192 = pkin(8) + t85;
t73 = t192 * t119;
t74 = t192 * t121;
t144 = -t190 * t74 - t191 * t73;
t39 = t190 * t73 - t191 * t74;
t225 = t39 * mrSges(6,1) - t144 * mrSges(6,2);
t169 = Ifges(6,5) * t99 + Ifges(6,6) * t100;
t10 = t169 + t225;
t224 = t10 * qJD(5);
t223 = t229 * Ifges(5,4);
t17 = t169 + t226;
t222 = t17 * qJD(5);
t165 = Ifges(6,4) * t100;
t216 = Ifges(6,1) - Ifges(6,2);
t93 = Ifges(6,4) * t99;
t129 = t99 * t93 + (-t216 * t99 - t165) * t100;
t221 = -Ifges(5,4) * t119 + (Ifges(5,1) - Ifges(5,2)) * t121;
t66 = -t100 * mrSges(6,1) + t99 * mrSges(6,2);
t214 = t66 * qJD(5);
t67 = -mrSges(6,1) * t99 - mrSges(6,2) * t100;
t211 = t223 + (pkin(4) * t67 + t221) * t119 + t129;
t204 = -mrSges(5,1) * t121 + mrSges(5,2) * t119;
t159 = t119 ^ 2 + t229;
t201 = (t159 * mrSges(5,3) - mrSges(4,2)) * t87 + (t100 * t49 + t50 * t99) * mrSges(6,3) + (-mrSges(4,1) + t67 + t204) * t88;
t199 = 0.2e1 * m(6);
t198 = m(6) / 0.2e1;
t195 = t99 / 0.2e1;
t188 = pkin(3) * t106;
t187 = pkin(4) * t119;
t186 = t121 * pkin(4);
t84 = -pkin(3) - t87;
t81 = t84 - t186;
t175 = t81 * t66;
t113 = -pkin(3) - t186;
t163 = t113 * t66;
t160 = t84 * t106;
t158 = qJD(1) * t228;
t157 = pkin(4) * t191;
t156 = pkin(4) * t190;
t147 = t159 * t87;
t145 = (t81 / 0.2e1 + t113 / 0.2e1) * t66;
t72 = t81 * t187;
t5 = m(6) * t72 + t160 + t175 + t211;
t136 = t5 * qJD(1);
t6 = m(6) * (t144 * t49 - t39 * t50 + t81 * t88) + m(5) * (t85 * t147 + t84 * t88) + t201;
t135 = -qJD(2) * t228 / 0.2e1 - t6 * qJD(1);
t8 = t129 + t175;
t134 = t8 * qJD(1);
t104 = t113 * t187;
t125 = (t104 + t72) * t198;
t2 = t125 + t67 * t187 + t175 / 0.2e1 + 0.2e1 * t93 * t195 - t188 / 0.2e1 + t160 / 0.2e1 + t163 / 0.2e1 + t223 + t221 * t119 + (-t165 + t216 * (-t99 / 0.2e1 - t195)) * t100 - t227;
t9 = m(6) * t104 + t163 - t188 + t211;
t132 = -t2 * qJD(1) - t9 * qJD(3);
t12 = t129 + t163;
t124 = t145 + t129;
t3 = t124 - t208;
t131 = -t3 * qJD(1) - t12 * qJD(3);
t128 = Ifges(5,5) * t121 - Ifges(5,6) * t119 + t169 + (t100 * t156 - t157 * t99) * mrSges(6,3);
t103 = (mrSges(6,1) * t190 + mrSges(6,2) * t191) * pkin(4);
t127 = t103 * qJD(4);
t98 = t103 * qJD(5);
t7 = t24 * qJD(3) * t199 / 0.4e1;
t4 = t124 + t208;
t1 = t125 + t145 + (-pkin(3) / 0.2e1 + t84 / 0.2e1) * t106 + t211 + t227;
t11 = [qJD(3) * t6 + qJD(4) * t5 + qJD(5) * t8, t7, t1 * qJD(4) + t4 * qJD(5) - t135 + (0.2e1 * (t113 * t88 + t143 * t49 - t50 * t75) * t198 + m(5) * (-pkin(3) * t88 + pkin(7) * t147) + t201) * qJD(3), t1 * qJD(3) + ((t144 * t190 + t191 * t39) * t197 + t128 + t204 * t85 + t225) * qJD(4) + t224 + t136, t4 * qJD(3) + t10 * qJD(4) + t134 + t224; t7, 0, t158 / 0.2e1, -t214 + (-t106 - t66 + (t100 * t157 + t99 * t156) * t199 / 0.2e1) * qJD(4), -qJD(4) * t66 - t214; t2 * qJD(4) + t3 * qJD(5) + t135, -t158 / 0.2e1, qJD(4) * t9 + qJD(5) * t12, ((t143 * t190 + t191 * t75) * t197 + t128 + t204 * pkin(7) + t226) * qJD(4) + t222 - t132, t17 * qJD(4) - t131 + t222; -t2 * qJD(3) - t136, 0, t132, -t98, -t98 - t127; -t3 * qJD(3) - t134, 0, t131, t127, 0;];
Cq = t11;
