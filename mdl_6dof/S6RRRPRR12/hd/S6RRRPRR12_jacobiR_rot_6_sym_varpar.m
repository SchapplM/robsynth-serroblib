% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function JR_rot = S6RRRPRR12_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:10:06
% EndTime: 2019-02-22 12:10:06
% DurationCPUTime: 0.11s
% Computational Cost: add. (194->32), mult. (275->62), div. (0->0), fcn. (402->10), ass. (0->39)
t156 = cos(pkin(6));
t158 = sin(qJ(2));
t162 = cos(qJ(1));
t164 = t162 * t158;
t159 = sin(qJ(1));
t161 = cos(qJ(2));
t166 = t159 * t161;
t147 = t156 * t164 + t166;
t157 = sin(qJ(3));
t160 = cos(qJ(3));
t155 = sin(pkin(6));
t168 = t155 * t162;
t141 = -t147 * t160 + t157 * t168;
t163 = t162 * t161;
t167 = t159 * t158;
t146 = -t156 * t163 + t167;
t154 = pkin(12) + qJ(5) + qJ(6);
t152 = sin(t154);
t153 = cos(t154);
t133 = t141 * t152 + t146 * t153;
t134 = t141 * t153 - t146 * t152;
t173 = t152 * t160;
t172 = t153 * t160;
t171 = t155 * t157;
t170 = t155 * t160;
t169 = t155 * t161;
t165 = t160 * t161;
t139 = -t147 * t157 - t160 * t168;
t149 = -t156 * t167 + t163;
t148 = t156 * t166 + t164;
t145 = t156 * t157 + t158 * t170;
t144 = t156 * t160 - t158 * t171;
t143 = t149 * t160 + t159 * t171;
t142 = t149 * t157 - t159 * t170;
t138 = -t145 * t153 + t152 * t169;
t137 = -t145 * t152 - t153 * t169;
t136 = t143 * t153 + t148 * t152;
t135 = -t143 * t152 + t148 * t153;
t1 = [t134, -t148 * t172 + t149 * t152, -t142 * t153, 0, t135, t135; t136, -t146 * t172 + t147 * t152, t139 * t153, 0, t133, t133; 0 (t152 * t158 + t153 * t165) * t155, t144 * t153, 0, t137, t137; -t133, t148 * t173 + t149 * t153, t142 * t152, 0, -t136, -t136; t135, t146 * t173 + t147 * t153, -t139 * t152, 0, t134, t134; 0 (-t152 * t165 + t153 * t158) * t155, -t144 * t152, 0, t138, t138; t139, -t148 * t157, t143, 0, 0, 0; t142, -t146 * t157, -t141, 0, 0, 0; 0, t157 * t169, t145, 0, 0, 0;];
JR_rot  = t1;
