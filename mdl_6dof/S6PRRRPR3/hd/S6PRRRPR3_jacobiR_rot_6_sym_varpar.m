% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:55
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRPR3_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:55:04
% EndTime: 2019-02-22 09:55:04
% DurationCPUTime: 0.08s
% Computational Cost: add. (123->31), mult. (214->65), div. (0->0), fcn. (313->10), ass. (0->37)
t153 = qJ(3) + qJ(4);
t151 = sin(t153);
t158 = sin(qJ(6));
t170 = t151 * t158;
t160 = cos(qJ(6));
t169 = t151 * t160;
t154 = sin(pkin(11));
t155 = sin(pkin(6));
t168 = t154 * t155;
t156 = cos(pkin(11));
t167 = t155 * t156;
t159 = sin(qJ(2));
t166 = t155 * t159;
t157 = cos(pkin(6));
t165 = t157 * t159;
t161 = cos(qJ(2));
t164 = t157 * t161;
t163 = t158 * t161;
t162 = t160 * t161;
t152 = cos(t153);
t147 = -t154 * t165 + t156 * t161;
t146 = t154 * t164 + t156 * t159;
t145 = t154 * t161 + t156 * t165;
t144 = t154 * t159 - t156 * t164;
t143 = t157 * t151 + t152 * t166;
t142 = t151 * t166 - t157 * t152;
t141 = t143 * t160;
t140 = t143 * t158;
t139 = t147 * t152 + t151 * t168;
t138 = t147 * t151 - t152 * t168;
t137 = t145 * t152 - t151 * t167;
t136 = t145 * t151 + t152 * t167;
t135 = t139 * t160;
t134 = t139 * t158;
t133 = t137 * t160;
t132 = t137 * t158;
t1 = [0, -t146 * t170 + t147 * t160, t134, t134, 0, t138 * t160 - t146 * t158; 0, -t144 * t170 + t145 * t160, t132, t132, 0, t136 * t160 - t144 * t158; 0 (t151 * t163 + t159 * t160) * t155, t140, t140, 0, t142 * t160 + t155 * t163; 0, -t146 * t169 - t147 * t158, t135, t135, 0, -t138 * t158 - t146 * t160; 0, -t144 * t169 - t145 * t158, t133, t133, 0, -t136 * t158 - t144 * t160; 0 (t151 * t162 - t158 * t159) * t155, t141, t141, 0, -t142 * t158 + t155 * t162; 0, -t146 * t152, -t138, -t138, 0, 0; 0, -t144 * t152, -t136, -t136, 0, 0; 0, t155 * t161 * t152, -t142, -t142, 0, 0;];
JR_rot  = t1;
