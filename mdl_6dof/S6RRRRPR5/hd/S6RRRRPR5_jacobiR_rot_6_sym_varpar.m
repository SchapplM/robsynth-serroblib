% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR5
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPR5_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:33:03
% EndTime: 2019-02-26 22:33:03
% DurationCPUTime: 0.08s
% Computational Cost: add. (123->32), mult. (182->32), div. (0->0), fcn. (270->8), ass. (0->34)
t140 = qJ(2) + qJ(3);
t138 = sin(t140);
t141 = sin(qJ(6));
t142 = sin(qJ(4));
t144 = cos(qJ(6));
t145 = cos(qJ(4));
t148 = t141 * t142 + t144 * t145;
t160 = t148 * t138;
t139 = cos(t140);
t146 = cos(qJ(1));
t152 = t146 * t145;
t143 = sin(qJ(1));
t156 = t143 * t142;
t132 = t139 * t156 + t152;
t153 = t146 * t142;
t155 = t143 * t145;
t133 = t139 * t155 - t153;
t159 = t132 * t144 - t133 * t141;
t157 = t143 * t139;
t154 = t146 * t139;
t150 = t132 * t141 + t133 * t144;
t134 = t139 * t153 - t155;
t135 = t139 * t152 + t156;
t121 = t134 * t144 - t135 * t141;
t122 = t134 * t141 + t135 * t144;
t149 = t141 * t145 - t142 * t144;
t127 = t149 * t138;
t130 = t148 * t139;
t129 = t149 * t139;
t126 = t146 * t160;
t125 = t146 * t127;
t124 = t143 * t160;
t123 = t143 * t127;
t1 = [-t150, -t126, -t126, -t121, 0, t121; t122, -t124, -t124, -t159, 0, t159; 0, t130, t130, t127, 0, -t127; -t159, t125, t125, t122, 0, -t122; t121, t123, t123, t150, 0, -t150; 0, -t129, -t129, t160, 0, -t160; t143 * t138, -t154, -t154, 0, 0, 0; -t146 * t138, -t157, -t157, 0, 0, 0; 0, -t138, -t138, 0, 0, 0;];
JR_rot  = t1;
