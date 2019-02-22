% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:22
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPR9_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:22:03
% EndTime: 2019-02-22 12:22:03
% DurationCPUTime: 0.08s
% Computational Cost: add. (129->30), mult. (214->58), div. (0->0), fcn. (313->10), ass. (0->38)
t149 = cos(pkin(6));
t150 = sin(qJ(2));
t153 = cos(qJ(1));
t155 = t153 * t150;
t151 = sin(qJ(1));
t152 = cos(qJ(2));
t156 = t151 * t152;
t139 = t149 * t155 + t156;
t145 = qJ(3) + qJ(4);
t143 = sin(t145);
t144 = cos(t145);
t147 = sin(pkin(6));
t158 = t147 * t153;
t130 = -t139 * t143 - t144 * t158;
t146 = sin(pkin(12));
t166 = t130 * t146;
t154 = t153 * t152;
t157 = t151 * t150;
t141 = -t149 * t157 + t154;
t159 = t147 * t151;
t133 = t141 * t143 - t144 * t159;
t165 = t133 * t146;
t160 = t147 * t150;
t136 = -t143 * t160 + t149 * t144;
t164 = t136 * t146;
t163 = t144 * t146;
t148 = cos(pkin(12));
t162 = t144 * t148;
t161 = t144 * t152;
t132 = -t139 * t144 + t143 * t158;
t140 = t149 * t156 + t155;
t138 = t149 * t154 - t157;
t137 = t149 * t143 + t144 * t160;
t135 = t136 * t148;
t134 = t141 * t144 + t143 * t159;
t129 = t133 * t148;
t128 = t130 * t148;
t1 = [t132 * t148 + t138 * t146, -t140 * t162 + t141 * t146, -t129, -t129, 0, 0; t134 * t148 + t140 * t146, t138 * t162 + t139 * t146, t128, t128, 0, 0; 0 (t146 * t150 + t148 * t161) * t147, t135, t135, 0, 0; -t132 * t146 + t138 * t148, t140 * t163 + t141 * t148, t165, t165, 0, 0; -t134 * t146 + t140 * t148, -t138 * t163 + t139 * t148, -t166, -t166, 0, 0; 0 (-t146 * t161 + t148 * t150) * t147, -t164, -t164, 0, 0; t130, -t140 * t143, t134, t134, 0, 0; t133, t138 * t143, -t132, -t132, 0, 0; 0, t147 * t152 * t143, t137, t137, 0, 0;];
JR_rot  = t1;
