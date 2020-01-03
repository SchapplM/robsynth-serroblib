% Return the minimum parameter vector for
% S5RRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% MPV [31x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RRRRR7_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR7_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR7_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR7_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t100 = (m(5) + m(6));
t112 = -pkin(8) * t100 - mrSges(5,3);
t93 = (mrSges(4,3) - t112);
t98 = (m(4) + t100);
t111 = -pkin(7) * t98 - t93;
t109 = 2 * pkin(9) * mrSges(6,3) + Ifges(6,2);
t108 = pkin(9) * m(6) + mrSges(6,3);
t107 = (mrSges(3,3) - t111);
t106 = (pkin(2) ^ 2);
t105 = (pkin(3) ^ 2);
t104 = (pkin(4) ^ 2);
t103 = (pkin(7) ^ 2);
t102 = (pkin(8) ^ 2);
t101 = pkin(9) ^ 2;
t97 = (t103 + t106);
t96 = (t102 + t105);
t94 = (m(3) + t98);
t1 = [Ifges(2,3) + Ifges(3,2) + Ifges(4,2) + t104 * m(6) + Ifges(5,2) + 2 * pkin(8) * mrSges(5,3) + t96 * t100 + 2 * pkin(7) * t93 + t97 * t98 + 2 * pkin(6) * t107 + (pkin(1) ^ 2 + pkin(6) ^ 2) * t94; pkin(1) * t94 + mrSges(2,1); -pkin(6) * t94 + mrSges(2,2) - t107; Ifges(3,1) - Ifges(3,2) + (t103 - t97) * t98; Ifges(3,4); t111 * pkin(2) + Ifges(3,5); Ifges(3,6); t106 * t98 + Ifges(3,3); pkin(2) * t98 + mrSges(3,1); mrSges(3,2); Ifges(4,1) - Ifges(4,2) + (t102 - t96) * t100; Ifges(4,4); t112 * pkin(3) + Ifges(4,5); Ifges(4,6); t105 * t100 + Ifges(4,3); pkin(3) * t100 + mrSges(4,1); mrSges(4,2); Ifges(5,1) - Ifges(5,2) + (t101 - t104) * m(6) + t109; t108 * pkin(4) + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t101 + t104) * m(6) + t109; m(6) * pkin(4) + mrSges(5,1); mrSges(5,2) - t108; Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
