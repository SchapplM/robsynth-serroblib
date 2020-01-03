% Return the minimum parameter vector for
% S5RRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% MPV [23x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RRPPR7_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_convert_par2_MPV_fixb: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR7_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR7_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR7_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t105 = (pkin(7) ^ 2);
t106 = (pkin(4) ^ 2);
t111 = 2 * pkin(7) * mrSges(6,3) + Ifges(6,2);
t95 = Ifges(4,2) + Ifges(5,3) + (t105 + t106) * m(6) + t111;
t97 = t106 * m(6) + Ifges(4,1) + Ifges(5,2);
t113 = -t95 + t97;
t101 = sin(pkin(8));
t98 = t101 ^ 2;
t102 = cos(pkin(8));
t99 = t102 ^ 2;
t112 = t99 - t98;
t110 = t101 * t102;
t109 = -pkin(7) * m(6) - mrSges(6,3);
t104 = Ifges(4,4) + Ifges(5,6);
t108 = t104 * t110;
t107 = t102 * mrSges(4,1) - t101 * mrSges(4,2);
t103 = Ifges(4,6) - Ifges(5,5);
t96 = t109 * pkin(4) - Ifges(5,4) + Ifges(4,5);
t1 = [Ifges(2,3) + Ifges(3,2) + t98 * t97 + 0.2e1 * t108 + t99 * t95 + (2 * pkin(6) * mrSges(3,3)) + ((pkin(1) ^ 2 + pkin(6) ^ 2) * m(3)); m(3) * pkin(1) + mrSges(2,1); -pkin(6) * m(3) + mrSges(2,2) - mrSges(3,3); t113 * t112 + Ifges(3,1) - Ifges(3,2) - 0.4e1 * t108; t112 * t104 + t113 * t110 + Ifges(3,4); -t101 * t103 + t102 * t96 + Ifges(3,5); t101 * t96 + t102 * t103 + Ifges(3,6); (m(6) * t105) + 0.2e1 * pkin(2) * t107 + Ifges(5,1) + Ifges(3,3) + Ifges(4,3) + t111; mrSges(3,1) + t107; t101 * mrSges(4,1) + t102 * mrSges(4,2) + mrSges(3,2); mrSges(4,3); m(4); m(6) * pkin(4) + mrSges(5,1); mrSges(5,2) + t109; mrSges(5,3); m(5) + m(6); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
