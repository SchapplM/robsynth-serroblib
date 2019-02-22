% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR12
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:10
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR12_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:10:16
% EndTime: 2019-02-22 12:10:16
% DurationCPUTime: 0.10s
% Computational Cost: add. (112->31), mult. (219->62), div. (0->0), fcn. (320->10), ass. (0->37)
t136 = cos(pkin(6));
t138 = sin(qJ(2));
t142 = cos(qJ(1));
t144 = t142 * t138;
t139 = sin(qJ(1));
t141 = cos(qJ(2));
t146 = t139 * t141;
t127 = t136 * t144 + t146;
t137 = sin(qJ(3));
t140 = cos(qJ(3));
t135 = sin(pkin(6));
t148 = t135 * t142;
t121 = -t127 * t140 + t137 * t148;
t143 = t142 * t141;
t147 = t139 * t138;
t126 = -t136 * t143 + t147;
t134 = pkin(12) + qJ(5);
t132 = sin(t134);
t133 = cos(t134);
t157 = t121 * t132 + t126 * t133;
t156 = t121 * t133 - t126 * t132;
t153 = t132 * t140;
t152 = t133 * t140;
t151 = t135 * t137;
t150 = t135 * t140;
t149 = t135 * t141;
t145 = t140 * t141;
t119 = -t127 * t137 - t140 * t148;
t129 = -t136 * t147 + t143;
t128 = t136 * t146 + t144;
t125 = t136 * t137 + t138 * t150;
t124 = t136 * t140 - t138 * t151;
t123 = t129 * t140 + t139 * t151;
t122 = t129 * t137 - t139 * t150;
t118 = t123 * t133 + t128 * t132;
t117 = -t123 * t132 + t128 * t133;
t1 = [t156, -t128 * t152 + t129 * t132, -t122 * t133, 0, t117, 0; t118, -t126 * t152 + t127 * t132, t119 * t133, 0, t157, 0; 0 (t132 * t138 + t133 * t145) * t135, t124 * t133, 0, -t125 * t132 - t133 * t149, 0; -t157, t128 * t153 + t129 * t133, t122 * t132, 0, -t118, 0; t117, t126 * t153 + t127 * t133, -t119 * t132, 0, t156, 0; 0 (-t132 * t145 + t133 * t138) * t135, -t124 * t132, 0, -t125 * t133 + t132 * t149, 0; t119, -t128 * t137, t123, 0, 0, 0; t122, -t126 * t137, -t121, 0, 0, 0; 0, t137 * t149, t125, 0, 0, 0;];
JR_rot  = t1;
