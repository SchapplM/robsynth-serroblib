% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PPRPRR1_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_jacobiR_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:39:48
% EndTime: 2019-02-26 19:39:48
% DurationCPUTime: 0.10s
% Computational Cost: add. (114->30), mult. (323->66), div. (0->0), fcn. (449->14), ass. (0->39)
t145 = sin(pkin(11));
t147 = sin(pkin(6));
t162 = t145 * t147;
t152 = cos(pkin(6));
t161 = t145 * t152;
t150 = cos(pkin(11));
t160 = t147 * t150;
t151 = cos(pkin(7));
t159 = t147 * t151;
t158 = t150 * t152;
t143 = sin(pkin(13));
t148 = cos(pkin(13));
t154 = sin(qJ(3));
t156 = cos(qJ(3));
t157 = t156 * t143 + t154 * t148;
t140 = t154 * t143 - t156 * t148;
t146 = sin(pkin(7));
t132 = t157 * t146;
t134 = t157 * t151;
t144 = sin(pkin(12));
t149 = cos(pkin(12));
t136 = -t145 * t144 + t149 * t158;
t137 = t144 * t158 + t145 * t149;
t124 = -t132 * t160 + t136 * t134 - t137 * t140;
t138 = -t150 * t144 - t149 * t161;
t139 = -t144 * t161 + t150 * t149;
t126 = t132 * t162 + t138 * t134 - t139 * t140;
t128 = t152 * t132 + (t134 * t149 - t140 * t144) * t147;
t155 = cos(qJ(5));
t153 = sin(qJ(5));
t135 = -t147 * t149 * t146 + t152 * t151;
t133 = t140 * t151;
t131 = t140 * t146;
t130 = -t138 * t146 + t145 * t159;
t129 = -t136 * t146 - t150 * t159;
t127 = -t152 * t131 + (-t133 * t149 - t144 * t157) * t147;
t125 = -t131 * t162 - t138 * t133 - t139 * t157;
t123 = t131 * t160 - t136 * t133 - t137 * t157;
t1 = [0, 0, t125 * t155, 0, -t126 * t153 + t130 * t155, 0; 0, 0, t123 * t155, 0, -t124 * t153 + t129 * t155, 0; 0, 0, t127 * t155, 0, -t128 * t153 + t135 * t155, 0; 0, 0, -t125 * t153, 0, -t126 * t155 - t130 * t153, 0; 0, 0, -t123 * t153, 0, -t124 * t155 - t129 * t153, 0; 0, 0, -t127 * t153, 0, -t128 * t155 - t135 * t153, 0; 0, 0, t126, 0, 0, 0; 0, 0, t124, 0, 0, 0; 0, 0, t128, 0, 0, 0;];
JR_rot  = t1;
