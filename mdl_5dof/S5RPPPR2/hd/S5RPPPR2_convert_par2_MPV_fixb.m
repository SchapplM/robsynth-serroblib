% Return the minimum parameter vector for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% MPV [22x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RPPPR2_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR2_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR2_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t127 = 2 * pkin(6) * mrSges(6,3) + Ifges(6,2);
t116 = sin(pkin(9));
t119 = cos(pkin(9));
t126 = t119 * t116;
t125 = pkin(6) * m(6) + mrSges(6,3);
t111 = t125 * pkin(4) + Ifges(5,4);
t124 = t111 * t126;
t123 = pkin(4) ^ 2;
t122 = pkin(6) ^ 2;
t121 = cos(pkin(7));
t120 = cos(pkin(8));
t118 = sin(pkin(7));
t117 = sin(pkin(8));
t114 = t119 ^ 2;
t113 = t116 ^ 2;
t112 = t123 * m(6) + Ifges(5,2);
t110 = m(6) * t122 + Ifges(5,1) + t127;
t1 = [Ifges(2,3) + t121 ^ 2 * (t113 * t110 + t114 * t112 + Ifges(3,2) + Ifges(4,3) + 0.2e1 * t124) + (0.2e1 * t121 * (Ifges(3,4) - t120 * (Ifges(4,5) + (t114 - t113) * t111 + (t110 - t112) * t126) + t117 * (-t116 * Ifges(5,5) - t119 * Ifges(5,6) + Ifges(4,6))) + (Ifges(3,1) + t120 ^ 2 * (t114 * t110 + t113 * t112 + Ifges(4,1) - 0.2e1 * t124) + (-0.2e1 * t120 * (-t119 * Ifges(5,5) + t116 * Ifges(5,6) + Ifges(4,4)) + (Ifges(4,2) + Ifges(5,3) + (t122 + t123) * m(6) + t127) * t117) * t117) * t118) * t118; mrSges(2,1); mrSges(2,2); mrSges(3,1); mrSges(3,2); mrSges(3,3); m(3); mrSges(4,1); mrSges(4,2); mrSges(4,3); m(4); m(6) * pkin(4) + mrSges(5,1); mrSges(5,2) - t125; mrSges(5,3); m(5) + m(6); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
