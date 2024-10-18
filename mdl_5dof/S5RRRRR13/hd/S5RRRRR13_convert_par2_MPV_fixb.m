% Return the minimum parameter vector for
% S5RRRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
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
% MPV [23x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RRRRR13_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR13_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR13_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR13_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t99 = m(5) + m(6);
t103 = (t99 * pkin(9));
t97 = m(4) + t99;
t102 = -pkin(10) * m(6) - mrSges(6,3);
t101 = (pkin(4) ^ 2);
t100 = (pkin(10) ^ 2);
t98 = sin(pkin(5));
t95 = (t100 + t101);
t94 = m(3) + t97;
t93 = (mrSges(5,3) - t102);
t1 = [pkin(1) ^ 2 * t94 + Ifges(2,3); pkin(1) * t94 + mrSges(2,1); mrSges(2,2); pkin(2) ^ 2 * t97 + Ifges(3,3); pkin(2) * t97 + mrSges(3,1); mrSges(3,2); pkin(3) ^ 2 * t99 + Ifges(4,3) + (t95 * m(6) + 2 * pkin(10) * mrSges(6,3) + Ifges(5,2) + Ifges(6,2) + (2 * t93 + t103) * pkin(9)) * t98 ^ 2; pkin(3) * t99 + mrSges(4,1); mrSges(4,2) + (-t93 - t103) * t98; Ifges(5,1) - Ifges(5,2) + (t100 - t95) * m(6); Ifges(5,4); t102 * pkin(4) + Ifges(5,5); Ifges(5,6); t101 * m(6) + Ifges(5,3); m(6) * pkin(4) + mrSges(5,1); mrSges(5,2); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
