% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR13
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:47
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR13_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:47:13
% EndTime: 2019-02-22 11:47:13
% DurationCPUTime: 0.12s
% Computational Cost: add. (147->33), mult. (275->63), div. (0->0), fcn. (402->10), ass. (0->39)
t156 = cos(pkin(6));
t161 = cos(qJ(2));
t162 = cos(qJ(1));
t163 = t162 * t161;
t158 = sin(qJ(2));
t159 = sin(qJ(1));
t166 = t159 * t158;
t146 = -t156 * t163 + t166;
t157 = sin(qJ(4));
t160 = cos(qJ(4));
t155 = sin(pkin(6));
t168 = t155 * t162;
t141 = -t146 * t157 + t160 * t168;
t164 = t162 * t158;
t165 = t159 * t161;
t147 = t156 * t164 + t165;
t154 = qJ(5) + qJ(6);
t152 = sin(t154);
t153 = cos(t154);
t134 = t141 * t152 + t147 * t153;
t135 = t141 * t153 - t147 * t152;
t173 = t152 * t157;
t172 = t153 * t157;
t171 = t155 * t158;
t170 = t155 * t160;
t169 = t155 * t161;
t167 = t157 * t158;
t140 = t146 * t160 + t157 * t168;
t149 = -t156 * t166 + t163;
t148 = t156 * t165 + t164;
t145 = t156 * t160 - t157 * t169;
t144 = -t156 * t157 - t160 * t169;
t139 = t148 * t157 + t159 * t170;
t138 = t159 * t155 * t157 - t148 * t160;
t137 = -t145 * t153 - t152 * t171;
t136 = -t145 * t152 + t153 * t171;
t133 = t139 * t153 + t149 * t152;
t132 = -t139 * t152 + t149 * t153;
t1 = [t135, -t148 * t152 + t149 * t172, 0, -t138 * t153, t132, t132; t133, -t146 * t152 + t147 * t172, 0, t140 * t153, t134, t134; 0 (t152 * t161 + t153 * t167) * t155, 0, t144 * t153, t136, t136; -t134, -t148 * t153 - t149 * t173, 0, t138 * t152, -t133, -t133; t132, -t146 * t153 - t147 * t173, 0, -t140 * t152, t135, t135; 0 (-t152 * t167 + t153 * t161) * t155, 0, -t144 * t152, t137, t137; t140, -t149 * t160, 0, t139, 0, 0; t138, -t147 * t160, 0, -t141, 0, 0; 0, -t158 * t170, 0, t145, 0, 0;];
JR_rot  = t1;
