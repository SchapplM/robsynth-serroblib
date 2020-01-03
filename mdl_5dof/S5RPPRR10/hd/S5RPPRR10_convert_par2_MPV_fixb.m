% Return the minimum parameter vector for
% S5RPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
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
% MPV [25x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RPPRR10_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_convert_par2_MPV_fixb: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR10_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR10_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR10_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t99 = (m(5) + m(6));
t100 = (pkin(7) ^ 2);
t102 = (pkin(4) ^ 2);
t106 = (Ifges(5,2) + (t100 + t102) * m(6));
t105 = -pkin(7) * m(6) - mrSges(6,3);
t104 = (mrSges(5,3) - t105);
t103 = 2 * pkin(7) * mrSges(6,3) + 2 * pkin(6) * t104 + Ifges(6,2) + t106;
t101 = pkin(6) ^ 2;
t98 = cos(pkin(8));
t97 = sin(pkin(8));
t1 = [Ifges(2,3) + t98 ^ 2 * (Ifges(3,2) + Ifges(4,3) + (pkin(3) ^ 2 + t101) * t99 + t103) + (0.2e1 * t98 * (Ifges(3,4) - Ifges(4,5)) + (t101 * t99 + Ifges(3,1) + Ifges(4,1) + t103) * t97) * t97; mrSges(2,1); mrSges(2,2); mrSges(3,1); mrSges(3,2); mrSges(3,3); m(3); pkin(3) * t99 + mrSges(4,1); -pkin(6) * t99 + mrSges(4,2) - t104; mrSges(4,3); m(4) + t99; m(6) * t100 + Ifges(5,1) - t106; Ifges(5,4); t105 * pkin(4) + Ifges(5,5); Ifges(5,6); t102 * m(6) + Ifges(5,3); m(6) * pkin(4) + mrSges(5,1); mrSges(5,2); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
